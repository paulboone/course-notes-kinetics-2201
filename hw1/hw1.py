
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt



def reactor_func(y, t, k):
    """
        Three functions:
            dCa/dt = -r
            dCb/dt = -r
            dCp/dt = 2r

        r = k[A][B] = kCaCb
        Define the right-hand side of equation dy/dt = a*y
    """

    c_a = y[0]
    c_b = y[1]
    c_p = y[2]
    r = k * c_a * c_b
    f = r * np.array([-1, -1, 2])
    # print("time: ", t,y, f)
    return f


# Initial concentrations
y0 = [1.0, 1.3, 0.0]

t_end = 50000
num_points = 1000
# Times at which the solution if to be computed.
t = np.linspace(0, t_end, num_points)

# Parameter value to use in `fun`.
# r = [-1, -1, 2]
k = 1.4e-4

# Solve the equation.
y = odeint(reactor_func, y0, t, args=(k,))


fig = plt.figure(figsize=(7,7))
ax = fig.add_subplot(1, 1, 1)
ax.plot(t, y[:,0], zorder=2)
ax.plot(t, y[:,1], zorder=2)
ax.plot(t, y[:,2], zorder=2)
ax.set_xlim(0,t_end)
ax.grid(linestyle='-', color='0.7', zorder=0)
ax.set_xlabel('t [seconds]')
ax.set_ylabel('concentration [M]')
ax.legend(["C_a","C_b", "C_p"])
conversion_90 = [i for i,yval in enumerate(y) if yval[0] < 0.1*y0[0]]
subheader = ""
if conversion_90:
    time_at_conversion_90 = conversion_90[0] * (t_end / num_points)
    subheader = "\nTime at 90%% conversion:  %s" % time_at_conversion_90
    ax.axvspan(0,time_at_conversion_90, color='0.9')

ax.set_title("Concentration vs Time" + subheader)

fig.savefig("conc-vs-time.png", dpi=144)
