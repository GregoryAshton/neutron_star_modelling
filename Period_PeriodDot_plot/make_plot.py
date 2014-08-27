import numpy as np
import matplotlib.pyplot as plt

DATA_FILE = "name_p0_p1.txt"

data = np.genfromtxt(DATA_FILE, usecols=(3, 6), skip_header=4, skip_footer=1,
                     autostrip=True)

p0 = data[:, 0]
p1 = data[:, 1]

print p1
fig = plt.figure()
ax = fig.add_subplot(111)

ax.scatter(p0, p1, marker=".")
ax.set_yscale("symlog")
ax.set_xscale("log")
plt.show()
