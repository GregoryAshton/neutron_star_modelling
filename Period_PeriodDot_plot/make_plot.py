import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc_file
rc_file("/home/greg/Neutron_star_modelling/matplotlibrc")

DATA_FILE = "name_p0_p1_AGE_BSURF.txt"

data = np.genfromtxt(DATA_FILE, usecols=(3, 6, 9, 10), skip_header=4, skip_footer=1,
                     autostrip=True)

p0 = data[:, 0]
p1 = data[:, 1]
AGE = data[:, 2]
BSURF = data[:, 3]

neg_idxs = p1 > 0
p0 = p0[neg_idxs]
p1 = p1[neg_idxs]

fig = plt.figure()
ax = fig.add_subplot(111)

# Millisecond pulsars
milliseconds_idxs = p0 < 1e-1
ax.loglog(p0[milliseconds_idxs], p1[milliseconds_idxs], 'o', 
          color='r', markersize=2.0, label="Millisecond pulsars")

# Normal pulsars
normal_idxs = p0 > 1e-1
ax.loglog(p0[normal_idxs], p1[normal_idxs], 'o', 
          color='k', markersize=2.0, label="Normal pulsars")

ax.legend(loc=2, frameon=False)
ax.set_xlabel("Period [s]")
ax.set_ylabel("Period deviative")

plt.savefig("Period_PeriodDot.pdf")
plt.show()
