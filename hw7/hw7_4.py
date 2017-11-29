
import math

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt


ua_v = 50
cp_f = 60
f_hx = 10
cp_hx = 7

pressure = 1 # atm
const_r = 0.082057 # L atm K−1 mol−1
const_r_J = 8.314 # J / mol K
c_h = 0.01692573879495781 / 2.0

k_f = 1e3
k_b = 1e3

dG = -76200 # J / mol
e_a_diff = dG + const_r_J * 298.15 * math.log(k_f / k_b)
e_a_diff
e_a_f = 5000 # J / mol
e_a_b = e_a_f - e_a_diff
e_a_b

desired_conversion = 0.5

temps = np.linspace(360, 1000, 641)
y = np.exp(-(-118.78 / (const_r_J * temps) + 0.14294  / const_r_J))

fig1 = plt.figure(figsize=(7,7))
ax = fig1.add_subplot(1, 1, 1)
ax.plot(temps, y, zorder=2)
ax.set_xlabel('temp [K]')
ax.set_ylabel('K_eq')
ax.grid(linestyle='-', color='0.7', zorder=0)
# fig.savefig("conc-vs-time.png", dpi=288)


t_opt = (e_a_b - e_a_f)  / (const_r_J * math.log(((e_a_b * k_b) / (e_a_f * k_f)) * (desired_conversion / (1 - desired_conversion))))
t_opt


def reactor_func(vars, t):
    """
        Define the right-hand side of equation dy/dt = a*y
    """
    f0, f1, f2, temp, temp_hx = vars
    k1 = 1
    # k1 = 1 * np.exp(-50000 / (8.31451 * temp))

    total_flow = f0 + f1 + f2  # mol / s
    ig_volumetric_flow = pressure / (const_r * temp) # mol / L
    r1 = k1 * f0 * f1 * (ig_volumetric_flow / total_flow) ** 2 # mol / L s

    dTdV = 0 #(30000 * r1  + 34000 * r2 - ua_v * (temp - temp_hx)) / ((f0 + f1 + f2 + f3 + f4 + f5) * cp_f)
    dTedVs = 0 # ua_v * (temp - temp_hx) / (f_hx * cp_hx)
    return np.array([-r1, -r1, r1, dTdV, dTedVs])

# Initial concentrations
f0 = [c_h, c_h, 0, 360.00, 463.15]

t_end = 20
num_points = 100
t = np.linspace(0, t_end, num_points)

# Solve the equation.
y = odeint(reactor_func, f0, t)

fig = plt.figure(figsize=(10,5))
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

conversion_85 = [i for i,yval in enumerate(y) if yval[0] < 0.15*f0[0]]
subheader = ""
temp_subheader = ""
if conversion_85:
    volume_at_conversion_85 = conversion_85[0] * (t_end / num_points)
    temp_at_conversion_85 = y[conversion_85[0],3]
    subheader = "\nVolume at 85%% conversion:  %2.1f L" % volume_at_conversion_85
    temp_subheader = "\nTemp at 85%% conversion:  %2.1f K" % temp_at_conversion_85
    ax.axvspan(0,volume_at_conversion_85, color='0.9')

ax.set_title("Flowrate vs Cumulative Volume (PFR)" + subheader)


# ax = fig.add_subplot(1, 2, 2)
# ax.plot(t, y[:,6], zorder=2)
# ax.plot(t, y[:,7], zorder=2)
# ax.legend(["T reactor", "T HX"])
# ax.set_xlabel('V [liters]')
# ax.set_ylabel('Temp [K]')
# ax.grid(linestyle='-', color='0.7', zorder=0)
# ax.set_xlim(0,t_end)
# ax.set_ylim(450, 700)
# temp_subheader = "\nUA/V = %2.1f" % ua_v
# ax.set_title('Temp vs Cumulative Volume (PFR)' + temp_subheader)
# ax.axhspan(450,673.15, color='0.9')
# fig.savefig("pfr-flow-vs-vol.png", dpi=288)
