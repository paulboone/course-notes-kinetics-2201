
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def reactor_func(concs, t):
    """
        Define the right-hand side of equation dy/dt = a*y
    """

    k1 = 4.64e-4 # 1/(M*s)
    r1 = -k1 * concs[0] * concs[1]

    return np.array([r1, r1, -r1, -r1])


# Initial concentrations
y0 = [1.0, 1.1, 0, 0]

t_end = 80000
num_points = 100
t = np.linspace(0, t_end, num_points)

# Solve the equation.
y = odeint(reactor_func, y0, t)

fig = plt.figure(figsize=(7,7))
ax = fig.add_subplot(1, 1, 1)
ax.plot(t, y[:,0], zorder=2)
ax.plot(t, y[:,1], zorder=2)
ax.plot(t, y[:,2], zorder=2)
ax.plot(t, y[:,3], zorder=2)
ax.set_xlim(0,t_end)
ax.grid(linestyle='-', color='0.7', zorder=0)
ax.set_xlabel('t [seconds]')
ax.set_ylabel('concentration [M]')
ax.legend(["C_a","C_b", "C_c", "C_d", "C_e", "C_p"])
conversion_limit = [i for i,yval in enumerate(y) if yval[0] < 0.05*y0[0]]
subheader = ""
if conversion_limit:
    time_at_conversion_limit = conversion_limit[0] * (t_end / num_points)
    subheader = "\nTime at 95%% conversion:  %s seconds" % time_at_conversion_limit
    ax.axvspan(0,time_at_conversion_limit, color='0.9')

ax.set_title("Concentration vs Time (w/o catalyst)" + subheader)
fig.savefig("conc-vs-time.png", dpi=144)
