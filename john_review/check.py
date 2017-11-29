
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def reactor_func(concs, t):
    """
        Define the right-hand side of equation dy/dt = a*y
    """
    x = concs[0]

    k1 = 3.18e-5

    # print(t, concs)
    r1 = -k1 * ((y0[0] / 2 - x)**2) / x

    return np.array([r1])


# Initial concentrations
y0 = [2]

t_end = 1e4*60
num_points = 100
t = np.linspace(0, t_end, num_points)

# Solve the equation.
y = odeint(reactor_func, y0, t)

fig = plt.figure(figsize=(7,7))
ax = fig.add_subplot(1, 1, 1)
ax.plot(t, y[:,0], zorder=2)
ax.set_xlim(0,t_end)
ax.set_ylim(0,y0[0])
ax.grid(linestyle='-', color='0.7', zorder=0)
ax.set_xlabel('t [seconds]')
ax.set_ylabel('concentration [M]')

ax.set_title("Concentration vs Time")
fig.savefig("conc-vs-time.png", dpi=144)
