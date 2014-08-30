import matplotlib.pyplot as plt
import numpy as np
from numpy import pi

N=1000

rand = np.random.normal(0, 0.9, N)
def phi(t, t0=0, f0=3, f1=-1e-1, f2=1e-2):
    dt = t - t0
    return 2 * pi * (f0 * dt + (f1 * dt **2.0)/2.0 + (f2 * dt **3.0)/6.0)

time = np.linspace(0, 1, N)
taylor = phi(time)

fig, ax2 = plt.subplots(figsize=(8, 4))
ax2.plot(np.random.normal(0, 1e-24, size=N), alpha=0.8, color="k")
ax2.plot(1e-24 * np.sin(taylor), lw=3.0, color="r")
ax2.set_xticks([])
ax2.set_xlabel("Time", size=14)
ax2.set_yticks([-1e-24, 0.0, 1e-24])
ax2.set_ylim(-1.2e-24, 1.2e-24)
ax2.set_ylabel("$h_{+}(t)$", rotation="horizontal", labelpad=20.5, size=20)
fig.tight_layout()
plt.savefig("CGW_example.pdf")
plt.show()
