import math

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

## protonation equilibrium constants (in M-1)
keq_co3_p = 1.95e10
keq_hco3_p = 2.46e3
keq_rnh2_p = 2.24e9
keq_oh_p = 6.76e13
keq_rnhcoo = 3.09e6#3.09e7

## protonation equilibrium rate constant
## (arbitrarily high to dwarf other rate constants and create near-instant equilibrium)
keq = 1e5

def reactor_func(c, t):
    """
        Define the right-hand side of equation dy/dt = a*y
    """

    co2g, co2l, h2o, oh, h, h2co3, hco3, co3, rnh2, rnh3, rnhcooh, rnhcoo = c

    r0 = r1 = r2 = r3 = r4 = r5 = r6 = r7 = r8 = r9 = r10 = 0

    r0 = 0
    r1 = 4.68e-2 * co2l - 40.65 * h2co3
    r2 = 1.24e4 * co2l * oh - 3.86e-4 * hco3
    r3 = keq * (keq_co3_p * co3 * h - hco3)
    r4 = keq * (keq_hco3_p * hco3 * h - h2co3)
    r5 = keq * (keq_rnh2_p * rnh2 * h - rnh3)
    r6 = keq * (keq_oh_p * oh * h - h2o)
    r7 = 9.16e+2 * h2co3 * rnh2 - 5.14e-3 * rnhcooh
    r8 = 1.05e-3 * hco3 * rnh2 - 7.43e-5 * rnhcoo
    r9 = 6.11e3 * co2l * rnh2 - 3*29.8 * rnhcooh
    r10 = keq * (keq_rnhcoo * rnhcoo * h - rnhcooh)

    # if h > 0.0 and oh > 0.0:
    #     print(h, oh, -math.log10(h))
    # print(co2l, h2co3, r1)
    # print(oh, h, h2o, keq_oh_p*oh*h, h2o)

    return np.array([r0, r0 - r1 - r2 - r9,
                    -r1 + r6 + r7 + r8, -r2 - r6, -r3 - r4 - r5 - r6 - r10,
                    r1 + r4 - r7, r2 + r3 - r4 - r8, -r3,
                    -r5 - r7 - r8 - r9, r5, r7 + r9 + r10, r8 - r10])


## Set Initial Concentrations
## based on intial concentrations of h2o, rnh2, co2l, we can calculate some initial equilibrium
## concentrations. See reaction specification spreadsheet for which index applies to which species.
# h2o0 = 55.5
# rnh20 = 0.004
# co2l0 = 0.006
# y0 = [0.0, co2l0,
#       h2o0,1*math.sqrt(h2o0 / keq_oh_p), 1*math.sqrt(h2o0 / keq_oh_p), 0.0, 0.0, 0.0,
#       rnh20, rnh20 / keq_rnh2_p, 0.0, 0.0]


h2o0 = 0
rnh20 = 0.004
co2l0 = 0.006
y0 = [0.0, co2l0,
      h2o0,1*math.sqrt(h2o0 / keq_oh_p), 1*math.sqrt(h2o0 / keq_oh_p), 0.0, 0.0, 0.0,
      rnh20, rnh20 / keq_rnh2_p, 0.0, 0.0]


# h2o0 = 52.5
# rnh20 = 3.01
# co2l0 = 0.0
# co30 = 1.2
# h0 = 2.83
#
# y0 = [0.0, co2l0,
#       h2o0,0, h0, 0.0, 0.0, co30,
#       rnh20, 0.0, 0.0, 0.0]

t_end = 0.2
num_points = 100
t = np.linspace(0, t_end, num_points)

# Solve the equation.
y = odeint(reactor_func, y0, t)
print(y[-1])

print(matplotlib.style.available)
# matplotlib.style.use("v2.0")
fig = plt.figure(figsize=(7,7))
ax = fig.add_subplot(1, 1, 1)
ax.plot(t, y[:,0], zorder=2, label="CO2(g)")
ax.plot(t, y[:,1], zorder=10, label="CO2(l)")
# ax.plot(t, y[:,2], zorder=2, label="H2O")
ax.plot(t, y[:,3], zorder=2, label="OH-")
ax.plot(t, y[:,4], zorder=2, label="H+")
ax.plot(t, y[:,5], zorder=2, label="H2CO3")
ax.plot(t, y[:,6], zorder=2, label="HCO3-")
ax.plot(t, y[:,7], zorder=2, label="CO3(2-)")
ax.plot(t, y[:,8], zorder=9, label="RNH2")
ax.plot(t, y[:,9], zorder=2, label="RNH3+")
ax.plot(t, y[:,10], zorder=2, label="RNHCOOH")
ax.plot(t, y[:,11], zorder=2, label="RNHCOO-")
ax.set_xlim(0,t_end)
ax.grid(linestyle='-', color='0.7', zorder=0)
ax.set_xlabel('t [seconds]')
ax.set_ylabel('concentration [M]')
ax.legend(loc='upper right')

ax2 = ax.twinx()
ax2.set_ylim(6, 11)
ax2.set_ylabel('pH')
ph = -np.log10(y[:,4])
ax2.plot(t, ph, '--', color='black', zorder=5, label="pH", linewidth=3)
ax2.legend(loc='center right')

fig.show()

# ax.set_title("CTime with Catalyst" + subheader)
# fig.savefig("conc-vs-time.png", dpi=144)
