import matplotlib.pyplot as plt
import numpy as np


deltaH = -14
deltaS = -27
r =
temps = np.linspace(300,800,501,)
k_eq = np.exp(-(deltaH - temps * deltaS) / (r * temps))
k_eq

fig = plt.figure(figsize=(7,7))
ax = fig.add_subplot(1, 1, 1)
ax.plot(t, y[:,0], zorder=2)
ax.set_xlabel('t [seconds]')
ax.set_ylabel('concentration [M]')

ax.grid(linestyle='-', color='0.7', zorder=0)
# ax.set_xlim(0,t_end)
# ax.legend(['C_a'])
# ax.set_title('Concentration vs Time' + subheader)
# fig.savefig("conc-vs-time.png", dpi=144)
