import matplotlib.pyplot as plt
import numpy as np
from numpy import pi

#plt.xkcd()

N=5000
x = np.arange(N)

rand = np.random.normal(0, 0.9, N)
def phi(t, t0=0, f0=3, f1=-1e-1, f2=1e-2):
    dt = t - t0
    return 2 * pi * (f0 * dt + (f1 * dt **2.0)/2.0 + (f2 * dt **3.0)/6.0)

time = np.linspace(0, 1, N)
taylor = phi(time)

signal = 1e-24 * np.sin(taylor)

fig, ax = plt.subplots(figsize=(8, 4))
noise_signal =  (signal + np.random.normal(0, 1e-24, size=N) 
                 )  
ax.plot(x, noise_signal, "o", alpha=0.8, markersize=3.0, color="k", label="Data")
ax.plot(x, signal, lw=3.0, color="r", label="Template")
ax.set_xticks([])
ax.set_xlabel("Time", size=14)
ax.set_yticks([-1e-24, 0.0, 1e-24])
ax.set_ylim(-5.2e-24, 5.2e-24)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
#ax.set_ylabel("$h_{+}(t)$", rotation="horizontal", labelpad=20.5, size=20)
ax.set_ylabel("Strain", size=14); ax.set_yticks([])
fig.tight_layout()
ax.legend(loc=1, frameon=True)
plt.savefig("CGW_example.pdf")
plt.show()
