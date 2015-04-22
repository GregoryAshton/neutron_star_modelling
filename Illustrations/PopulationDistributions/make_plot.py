import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc_file
import pandas as pd
import sys

rc_file("/home/greg/Neutron_star_modelling/matplotlibrc")

DATA_FILE = sys.argv[1]

data = np.genfromtxt(DATA_FILE, skip_header=4, skip_footer=1, dtype=None)
name = data[:, 0]
F0 = np.genfromtxt(data[:, 1])
F1 = np.genfromtxt(data[:, 2])
F2 = np.genfromtxt(data[:, 3])

# Binary = data[:, 3] # If not "*" then it is a binary
# Type = data[:, 4]
# Age = np.genfromtxt(data[:, 5])

df = pd.DataFrame({'name': name,
                   'F0': F0,
                   'F1': F1,
                   'F2': F2,
                   })


fig = plt.figure(figsize=(6, 8))
ax = fig.add_subplot(111)
ax.hist(np.log10(df.F0[np.isfinite(df.F0.values)]), bins=100, normed=True)
ax.set_xlabel("$\log_{10}(f)$")
ax.set_ylabel("Normalised count")
fig.tight_layout()
plt.savefig("F0_distribution.pdf")
plt.show()


fig = plt.figure(figsize=(6, 8))
ax = fig.add_subplot(111)
F1 = np.abs(df.F1[np.isfinite(df.F1.values)]).values
F1 = F1[F1 > 0]
ax.hist(np.log10(F1), bins=100, normed=True)
ax.set_xlabel("$\log_{10}(\dot{f})$")
ax.set_ylabel("Normalised count")
fig.tight_layout()
plt.savefig("F1_distribution.pdf")
plt.show()

fig = plt.figure(figsize=(6, 8))
ax = fig.add_subplot(111)
F2 = np.abs(df.F2[np.isfinite(df.F2.values)]).values
ax.hist(np.log10(F2), bins=100, normed=True)
ax.set_xlabel("$\log_{10}(\ddot{f})$")
ax.set_ylabel("Normalised count")
fig.tight_layout()
plt.savefig("F2_distribution.pdf")
plt.show()


