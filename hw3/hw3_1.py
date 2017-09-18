import matplotlib.pyplot as plt
import numpy as np

temp = np.array([2, 8, 14, 20, 25, 32, 40, 45, 50]) + 273.15 # in Kelvins
rate = np.array([0.0011, 0.0015, 0.0023, 0.0033, 0.0055, 0.0076, 0.0127, 0.0171, 0.0219])
ln_rate = np.log(rate)
fc = np.polyfit(1/temp, ln_rate, 1)
inv_temp_bounds = np.array([0.0, 0.004])
fcf = lambda invt: fc[1] + fc[0]*invt #I'm sure there is a better way to do this with linear algebra...
rate_fit = [ fcf(t) for t in inv_temp_bounds]

fig = plt.figure(figsize=(7,7))
fig.subplots_adjust(top=0.85)
ax = fig.add_subplot(1, 1, 1)
ax.plot(1/temp, ln_rate, '+', zorder=4)
ax.plot(inv_temp_bounds, rate_fit, zorder=3)
ax.plot([1/373.15], [fcf(1/373.15)], 'r^', zorder=5)
ax.grid(linestyle='-', color='0.7', zorder=0)
ax.set_xlim(inv_temp_bounds)
ax.set_xlabel('1 / temp [K]')
ax.set_ylabel('ln(rate)')
ax.legend(['Measured values (CSTR)', 'Linear Fit', 'Extrapolated k @ 100°C'])
ax.set_title("Arrhennius Plot\ny-intercept = %2.2f\nslope (EA / R) = %2.2f K\nln(rate) @ 100°C = %2.2f" % (fc[1], fc[0], fcf(1/373.15)))
fig.savefig("arrhenius plot", dpi=144)
