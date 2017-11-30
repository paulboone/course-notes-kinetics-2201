
import math

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

pressure0 = 10 # atm
temp0 = 360.0
const_r = 0.082057 # L atm K−1 mol−1
const_r_J = 8.314 # J / mol K
c_h = 0.01692573879495781 / 2.0

# prefactors
k_f = 1e1
k_b = 1e-1

dH = -118780 # J / mol
dG = -76200 # J / mol
e_a_diff = dG + const_r_J * 298.15* math.log(k_f / k_b)
e_a_diff
e_a_f = 10000 # J / mol
e_a_b = e_a_f - e_a_diff
e_a_b

desired_conversion = 0.1

temps = np.linspace(360, 1000, 641)
y = np.exp(-(-118780 / (const_r_J * temps) + 142.94 / const_r_J))

# fig1 = plt.figure(figsize=(7,7))
# ax = fig1.add_subplot(1, 1, 1)
# ax.plot(temps, y, zorder=2)
# # ax.plot(temps, k1, zorder=2)
# # ax.plot(temps, k2, zorder=2)
# ax.set_xlabel('temp [K]')
# ax.set_ylabel('K_eq')
# ax.grid(linestyle='-', color='0.7', zorder=0)
# ax.set_ylim(0.1,1e9)
# ax.set_yscale('log')
#
# fig1.savefig("Keq graph.png", dpi=288)

#
# t_opt = (e_a_b - e_a_f)  / (const_r_J * math.log(((e_a_b * k_b) / (e_a_f * k_f)) / (desired_conversion / (1 - desired_conversion))))
# t_opt


ua_v = 50 # J / L s K rough low estimate for gas inside tube, liquid outside from https://www.engineersedge.com/thermodynamics/overall_heat_transfer-table.htm
cp_f = 148.64 # J/mol*K cyclohexane at 400K
f_hx = 1
cp_hx = 70



def reactor_func(vars, t):
    """
        Define the right-hand side of equation dy/dt = a*y
    """
    f0, f1, f2, temp, temp_hx, pressure = vars


    k1 = k_f * np.exp(-e_a_f / (const_r_J * temp))
    k2 = k_b * np.exp(-e_a_b / (const_r_J * temp))

    if (t == 0.0):
        print(t, k1,k2)

    total_flow = f0 + f1 + f2  # mol / s
    ig_volumetric_flow = pressure / (const_r * temp) # mol / L
    # r1 = (k1 * f0 * f1  - k2 * f2) * (ig_volumetric_flow / total_flow) ** 2  # mol / L s

    ka = 2 # L / mole
    kb = 2 # L / mole
    pa = pressure * f0 / total_flow # mole / L
    pb = pressure * f1 / total_flow # mole / L
    krxn = 5e-2

    beta0=1 #       %Ergun parameter - calculated separately
    # Ar=1.14e-3 #      %Tube cross section in m2

    r1 = krxn * ka*pa*kb*pb / (1 + ka*pa + kb*pb)

    dTdV = ((-dH * r1) - ua_v * (temp - temp_hx)) / (total_flow * cp_f)
    dTedVs = ua_v * (temp - temp_hx) / (f_hx * cp_hx)
    dPdz = -beta0 * (pressure0 / pressure) * (temp / temp0) * (total_flow / (c_h + c_h))
    # print(dPdz, pressure0, pressure, temp0, temp, c_h * 2, total_flow)
    return np.array([-r1, -r1, r1, dTdV, dTedVs, dPdz])

# Initial concentrations
f0 = [c_h, c_h, 0, temp0,  400.00, pressure0]

t_end = 1
num_points = 1000
t = np.linspace(0, t_end, num_points)

# Solve the equation.
y = odeint(reactor_func, f0, t)

fig = plt.figure(figsize=(11,5))
ax = fig.add_subplot(1, 2, 1)
ax.plot(t, y[:,0], zorder=2, label="F Cyclohexene")
ax.plot(t, y[:,1], zorder=2, label="H2")
ax.plot(t, y[:,2], zorder=2, label="Cyclohexane")


ax.set_xlim(0, 0.010)
ax.set_xlim(0,t_end)
ax.grid(linestyle='-', color='0.7', zorder=0)
ax.set_xlabel('V [liters]')
ax.set_ylabel('Flow rate [mol / s]')
ax.legend(loc='upper right')


ax2 = ax.twinx()
ax2.plot(t, y[:,3], zorder=2, label='Temp', color='black')
ax2.set_ylim(300, 900)
ax2.set_ylabel('Temp')
ax2.legend(loc='right')

# ax3 = ax.twinx()
# ax2.plot(t, y[:,5], zorder=2, label='Pressure', color='grey')

conversion_xx = [i for i,yval in enumerate(y) if yval[0] < (1 - desired_conversion) * f0[0]]
subheader = ""
temp_subheader = ""
if conversion_xx:
    volume_at_conversion_xx = conversion_xx[0] * (t_end / num_points)
    temp_at_conversion_xx = y[conversion_xx[0],3]
    subheader = "\nVolume at %2.0f%% conversion:  %2.3f L" % (desired_conversion * 100, volume_at_conversion_xx)
    temp_subheader = "\nTemp at %2.0f%% conversion:  %2.3f K" % (desired_conversion * 100, temp_at_conversion_xx)
    ax.axvspan(0,volume_at_conversion_xx, color='0.9')

ax.set_title("Flowrate vs Cumulative Volume (PFR)" + subheader)


ax = fig.add_subplot(1, 2, 2)
fig.subplots_adjust(wspace=0.4, hspace=0.4)
ax.plot(t, y[:,5], zorder=2)
ax.set_xlabel('V [liters]')
ax.set_ylabel('Pressure')

fig.savefig("pfr-flow-vs-vol.png", dpi=288)
