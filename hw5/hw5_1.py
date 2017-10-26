
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt


def reactor_func(vars, t):
    """
        Define the right-hand side of equation dy/dt = a*y
    """
    f0, f1, f2, f3, f4, temp = vars
    k1 = 5e4
    k1r = k1 / 10
    k2 = 5e3

    ig_volumetric_flow = (1 / (82.0578 * 473.15)) * 1000 # mol / L
    total_flow = f0 + f1 + f2 + f3 + f4 # mol / s
    r1 = (k1 * f0**2  - k1r * f1 * f2) * (ig_volumetric_flow / total_flow) ** 2 # mol / L s
    r2 = k2 * f2 * f3 * (ig_volumetric_flow / total_flow) ** 2 # mol / L s

    return np.array([-2*r1, r1, r1 - r2, -r2, 2*r2, 0])

def reactor_func_esa(vars, t):
    """
        Define the right-hand side of equation dy/dt = a*y
    """
    f0, f1, f2, f3, f4, temp = vars
    k1 = 5e4
    k1r = k1 / 10
    k2 = 5e3

    f0start = 0.5
    k_eq = k1 / k1r

    ig_volumetric_flow = (1 / (82.0578 * 473.15)) * 1000 # mol / L
    total_flow = f0 + f1 + f2 + f3 + f4 # mol / s
    c_b = k_eq * f0start / 2 # 2 for stoich
    r2 =k2 * (f0**2 * f3) * (ig_volumetric_flow / total_flow)**2 # mol / L s

    return np.array([-2*r2, 0, 0, -r2, 2*r2, 0])

# Initial concentrations
f0 = [0.5, 0.0, 0.0, 0.8, 0.0, 473.15]

t_end = 20
num_points = 100
t = np.linspace(0, t_end, num_points)

# Solve the equation.
y = odeint(reactor_func, f0, t)

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(2, 1, 1)
ax.plot(t, y[:,0], zorder=2)
ax.plot(t, y[:,1], zorder=2)
ax.plot(t, y[:,2], zorder=2)
ax.plot(t, y[:,3], zorder=2)
ax.plot(t, y[:,4], zorder=2)

ax.set_xlim(0,t_end)
ax.set_ylim(0,1.0)
ax.grid(linestyle='-', color='0.7', zorder=0)
ax.set_xlabel('V [liters]')
ax.set_ylabel('Flow rate [mol / s]')
ax.legend(["F_a","F_b", "F_c", "F_d", "F_e"])
ax.set_title("Flowrate vs Cumulative Volume (PFR): Full Simulation")


# Solve the equation.
y = odeint(reactor_func_esa, f0, t)

ax = fig.add_subplot(2, 1, 2)
ax.plot(t, y[:,0], zorder=2)
ax.plot(t, y[:,1], zorder=2)
ax.plot(t, y[:,2], zorder=2)
ax.plot(t, y[:,3], zorder=2)
ax.plot(t, y[:,4], zorder=2)

ax.set_xlim(0,t_end)
ax.set_ylim(0,1.0)
ax.grid(linestyle='-', color='0.7', zorder=0)
ax.set_xlabel('V [liters]')
ax.set_ylabel('Flow rate [mol / s]')
ax.legend(["F_a","F_b", "F_c", "F_d", "F_e"])
ax.set_title("Flowrate vs Cumulative Volume (PFR): ESA")

fig.savefig("pfr-flow-vs-vol.png", dpi=288)
