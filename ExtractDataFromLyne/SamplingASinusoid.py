""" Script to investigate sampling from a sin wave """

import numpy as np
import matplotlib.pyplot as plt
import TNtools as TN
TN.PlotDefaults()

# Data
N = 10000
x = np.linspace(0, 100, N)
y = np.sin(x)

# Add experimental noise
noise = np.random.normal(0, 0.1, N)
y += noise

# Randomly sample points
SampleSize = 10000
SampleIndexes = np.random.randint(0, N, SampleSize)
SampleIndexes = np.sort(SampleIndexes)

y = y[SampleIndexes]
x = x[SampleIndexes]

fig, (ax1, ax2) = plt.subplots(nrows=2)

ax1.plot(x, y, ".")
ax1.set_ylabel("$\sin(x)$ + noise")
ax1.set_xlabel("$x$", labelpad=0.01)


ax2.hist(y, bins=100, color="b")
ax2.set_xlim(-2, 2)
ax2.set_xlabel("Samples from $sin(x)$")
ax2.set_ylabel("Count")

fig.subplots_adjust(hspace=0.5)
plt.savefig("img/Sinx_example.pdf")
plt.show()

