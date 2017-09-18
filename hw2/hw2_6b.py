
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def reactor_func(vars, t):
    """
        Define the right-hand side of equation dy/dt = a*y
    """
    f0, f1, f2, f3, f4, f5 = vars

    k1 = 12000.0 # L/(mole*s)
    k2 = 18000.0 # L/(mole*s)

    ig_volumetric_flow = (1 / (82.0578 * 473.15)) * 1000 # mol / L
    # print(ig_volumetric_flow)
    total_flow = f0 + f1 + f2 + f3 + f4 + f5 # mol / s
    r1 = k1 * f0 * f1 * (ig_volumetric_flow / total_flow) ** 2 # mol / L s
    r2 = k2 * f2 * f4 * (ig_volumetric_flow / total_flow) ** 2 # mol / L s

    return np.array([-r1, -r1, r1 - r2, r1, -r2, r2])

# Initial concentrations
y0 = [0.5, 0.5, 0, 0, 0.8, 0]

t_end = 5
num_points = 100
t = np.linspace(0, t_end, num_points)

# Solve the equation.
y = odeint(reactor_func, y0, t)

# print(y)

fig = plt.figure(figsize=(7,7))
ax = fig.add_subplot(1, 1, 1)
ax.plot(t, y[:,0], zorder=2)
ax.plot(t, y[:,1], zorder=2)
ax.plot(t, y[:,2], zorder=2)
ax.plot(t, y[:,3], zorder=2)
ax.plot(t, y[:,4], zorder=2)
ax.plot(t, y[:,5], zorder=2)
ax.set_xlim(0,t_end)
ax.grid(linestyle='-', color='0.7', zorder=0)

ax.set_xlabel('cumulative volume [liters]')
ax.set_ylabel('Flow rate [mol / s]')

ax.legend(["F_a","F_b", "F_c", "F_d", "F_e", "F_p"])
conversion_85 = [i for i,yval in enumerate(y) if yval[0] < 0.15*y0[0]]
subheader = ""
if conversion_85:
    time_at_conversion_85 = conversion_85[0] * (t_end / num_points)
    subheader = "\nV at 85%% conversion:  %s" % time_at_conversion_85
    ax.axvspan(0,time_at_conversion_85, color='0.9')

ax.set_title("Concentration vs Cumulative Volume\nTotal PFR volume=%.2f liters" % time_at_conversion_85)



# fig.savefig("conc-vs-time.png", dpi=144)

## flow rate vs cumulative volume
#
# fig = plt.figure(figsize=(7,7))
# ax = fig.add_subplot(1, 1, 1)
# # y = y / 69.8
# t = t * 69.8
# ax.plot(t, y[:,0], zorder=2)
# ax.plot(t, y[:,1], zorder=2)
# ax.plot(t, y[:,2], zorder=2)
# ax.plot(t, y[:,3], zorder=2)
# ax.plot(t, y[:,4], zorder=2)
# ax.plot(t, y[:,5], zorder=2)
# ax.set_xlim(0,0.07 * 69.8)
# ax.grid(linestyle='-', color='0.7', zorder=0)
# ax.set_xlabel('cumulative volume [liters]')
# # ax.set_ylabel('flow rate [moles / second]')
# ax.set_ylabel('concentration [M]')
# ax.legend(["C_a","C_b", "C_c", "C_d", "C_e", "C_p"])
# ax.axvspan(0,time_at_conversion_85*69.8, color='0.9')
# ax.set_title("Concentration vs Cumulative Volume\nTotal PFR volume=%.2f liters" % (time_at_conversion_85*69.8))
# fig.savefig("conc-vs-cum_volume.png", dpi=144)

# ax.legend(["F_a","F_b", "F_c", "F_d", "F_e", "F_p"])


# fig.savefig("conc-vs-time.png", dpi=144)
