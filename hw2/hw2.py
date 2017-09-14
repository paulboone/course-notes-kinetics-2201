
import numpy as np
import matplotlib.pyplot as plt
# from matplotlib import rc
# plt.rcParams["font.family"] = "sans-serif"
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# rc('text', usetex=True)

c_a = np.array([2, 1.8, 1.6, 1.4, 1.2, 1, 0.8, 0.6, 0.4, 0.2])
c_b = np.array([2.2, 2.0, 1.8, 1.4, 1.2, 1.0, 0.8, 0.6, 0.4])
c_c = np.array([0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8])
c_d = np.array([0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8])

x_a = (c_a[0] - c_a) / c_a[0]


feed_rate = 10000 / (24 * 60 * 60)
print(feed_rate, " moles / s")

r = np.array([0.023, 0.061, 0.12, 0.081, 0.077, 0.069, 0.045, 0.029, 0.021, 0.017])

delta_X = 0.1
volume_pfr = 0.0
fa0_r = feed_rate / r
print("feed_rate / r: ", fa0_r)
for i in range(0,len(fa0_r) - 1):
    print(((fa0_r[i] + fa0_r[i+1])/2) * delta_X)
    volume_pfr += ((fa0_r[i] + fa0_r[i+1])/2) * delta_X

print("volume PFR = %f" % volume_pfr)
print("volume CSTR = %f" % (fa0_r[-1] * 0.9))

delta_X = 0.1
r_inv = 0.0

print("1 / r: ", 1/r)
for i in range(0,len(r) - 1):
    print((((1/r)[i] + (1/r)[i+1])/2) * delta_X)
    r_inv += (((1/r)[i] + (1/r)[i+1])/2) * delta_X

print("1/r = %f" % r_inv)



fig = plt.figure(figsize=(7,7))
ax = fig.add_subplot(1, 1, 1)
ax.plot(x_a, 1/r, zorder=2)
ax.grid(linestyle='-', color='0.7', zorder=0)
ax.set_xlabel('$X_A$')
ax.set_ylabel("1/-r_A [-] $")


fig = plt.figure(figsize=(7,7))
ax = fig.add_subplot(1, 1, 1)
ax.plot(x_a, feed_rate/r, zorder=2)
ax.grid(linestyle='-', color='0.7', zorder=0)
ax.set_xlabel('$X_A$')
ax.set_ylabel("$ \\frac{F_{A0}}{-r_A} [liters] $")


fig.savefig("levonspiel plot.png", dpi=144)


 # The target conversion is 90% with respect to the reactant A and at a feed rate of 10 kmol/day of
 # the same reactant A.


# 10 kmol / 1 day (24 hours * 60 minutes * 60 second)

# volumetric_flow_rate =
