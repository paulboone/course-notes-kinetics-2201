
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt


ua_v = 50
cp_f = 60
f_hx = 10
cp_hx = 70


def reactor_func(vars, t):
    """
        Define the right-hand side of equation dy/dt = a*y
    """
    f0, f1, f2, f3, f4, f5, temp, temp_hx = vars
    k1 = 1.8e6 * np.exp(-26000 / (8.31451 * temp))
    k2 = 2.4e6 * np.exp(-32000 / (8.31451 * temp))

    ig_volumetric_flow = (1 / (82.0578 * 473.15)) * 1000 # mol / L
    total_flow = f0 + f1 + f2 + f3 + f4 + f5 # mol / s
    r1 = k1 * f0 * f1 * (ig_volumetric_flow / total_flow) ** 2 # mol / L s
    r2 = k2 * f2 * f4 * (ig_volumetric_flow / total_flow) ** 2 # mol / L s

    dTdV = (30000 * r1  + 34000 * r2 - ua_v * (temp - temp_hx)) / ((f0 + f1 + f2 + f3 + f4 + f5) * cp_f)
    dTedVs =  ua_v * (temp - temp_hx) / (f_hx * cp_hx)
    return np.array([-r1, -r1, r1 - r2, r1, -r2, r2, dTdV, dTedVs])

# Initial concentrations
f0 = [0.5, 0.5, 0.0, 0.0, 0.8, 0, 473.15, 463.15]

t_end = 20
num_points = 100
t = np.linspace(0, t_end, num_points)

# Solve the equation.
y = odeint(reactor_func, f0, t)

fig = plt.figure(figsize=(10,5))
ax = fig.add_subplot(1, 2, 1)
ax.plot(t, y[:,0], zorder=2)
ax.plot(t, y[:,1], zorder=2)
ax.plot(t, y[:,2], zorder=2)
ax.plot(t, y[:,3], zorder=2)
ax.plot(t, y[:,4], zorder=2)
ax.plot(t, y[:,5], zorder=2)
ax.set_xlim(0,t_end)
ax.grid(linestyle='-', color='0.7', zorder=0)
ax.set_xlabel('V [liters]')
ax.set_ylabel('Flow rate [mol / s]')
ax.legend(["F_a","F_b", "F_c", "F_d", "F_e", "F_p"])

conversion_85 = [i for i,yval in enumerate(y) if yval[0] < 0.15*f0[0]]
subheader = ""
temp_subheader = ""
if conversion_85:
    volume_at_conversion_85 = conversion_85[0] * (t_end / num_points)
    temp_at_conversion_85 = y[conversion_85[0],6]
    subheader = "\nVolume at 85%% conversion:  %2.1f L" % volume_at_conversion_85
    temp_subheader = "\nTemp at 85%% conversion:  %2.1f K" % temp_at_conversion_85
    ax.axvspan(0,volume_at_conversion_85, color='0.9')

ax.set_title("Flowrate vs Cumulative Volume (PFR)" + subheader)


ax = fig.add_subplot(1, 2, 2)
ax.plot(t, y[:,6], zorder=2)
ax.plot(t, y[:,7], zorder=2)
ax.legend(["T reactor", "T HX"])
ax.set_xlabel('V [liters]')
ax.set_ylabel('Temp [K]')
ax.grid(linestyle='-', color='0.7', zorder=0)
ax.set_xlim(0,t_end)
ax.set_ylim(450, 700)
temp_subheader = "\nUA/V = %2.1f" % ua_v
ax.set_title('Temp vs Cumulative Volume (PFR)' + temp_subheader)
ax.axhspan(450,673.15, color='0.9')
fig.savefig("pfr-flow-vs-vol.png", dpi=288)
