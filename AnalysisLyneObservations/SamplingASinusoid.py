""" Script to investigate sampling from a sin wave """

import numpy as np
import matplotlib.pyplot as plt
import TNtools as TN
import sys


TN.PlotDefaults()

# Data
N = 10000
y0 = 0.0
A = 1.0
x = np.linspace(0, 100, N)
y = y0 + A * np.sin(x)

if "noise" in sys.argv:
    noise = np.random.normal(0, 0.1, N)
    label = "with_noise"
else:
    noise = 0
    label = "without_noise"

y = y + noise

# Randomly sample points
SampleSize = 4000
SampleIndexes = np.random.randint(0, N, SampleSize)
SampleIndexes = np.sort(SampleIndexes)

y = y[SampleIndexes]
x = x[SampleIndexes]

fig, (ax1, ax2) = plt.subplots(nrows=2)

ax1.plot(x, y, "o", markersize=2.0)
ax1.set_ylabel("$\sin(x)$ "+label.replace("_", " "))
ax1.set_xlabel("$x$", labelpad=0.01)
ax1.set_ylim(1.5*(y0-A), 1.5*(y0 + A))


ax2.hist(y, bins=100, color="b")
ax2.set_xlim(-2, 2)
ax2.set_xlabel("Measured values of $y=\sin(x)$")
ax2.set_ylabel("Count")
ax2.set_xlim(1.5*(y0-A), 1.5*(y0 + A))

fig.subplots_adjust(hspace=0.5)
plt.savefig("img/Sinx_"+label+".pdf")
plt.show()

