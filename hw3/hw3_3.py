import matplotlib.pyplot as plt
import numpy as np


deltaH = -14000
deltaS = -27
r = 8.31451
temps = np.linspace(300,800,501,)
k_eq = np.exp(-(deltaH - temps * deltaS) / (r * temps))


fig = plt.figure(figsize=(7,7))
fig.subplots_adjust(hspace=0.4)
ax = fig.add_subplot(2, 1, 1)
ax.plot(temps, k_eq, zorder=2)
ax.set_xlabel('temp [K]')
ax.set_ylabel('k_eq')
ax.set_title("k_eq vs Temp")
ax.grid(linestyle='-', color='0.7', zorder=0)

ax = fig.add_subplot(2, 1, 2)
ax.plot(temps, k_eq / (k_eq + 1), zorder=2)
ax.set_xlabel('temp [K]')
ax.set_ylabel('X')
ax.set_title("X vs Temp")
ax.grid(linestyle='-', color='0.7', zorder=0)
fig.savefig("X vs Temp", dpi=288)
