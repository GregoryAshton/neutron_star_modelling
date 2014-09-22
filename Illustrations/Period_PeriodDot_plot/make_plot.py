import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc_file
import TNtools as TN

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

fig = plt.figure(figsize=(6,8))
ax = fig.add_subplot(111)

xlims = (1e-4, 100)

# Millisecond pulsars
milliseconds_idxs = p0 < 1e-1
ax.loglog(p0[milliseconds_idxs], p1[milliseconds_idxs], 'o', 
          color='k', markersize=3.0, label="Millisecond pulsars")

# Normal pulsars
normal_idxs = p0 > 1e-1
ax.loglog(p0[normal_idxs], p1[normal_idxs], 'o', 
          color='k', markersize=3.0, label="Normal pulsars")

P_range = np.logspace(-4, 2, 100)
# Add lines of constant magnetic field
def Pdot(P, B0):
    CONST = 3.2e19
    return (B0 / CONST) ** 2 / P  
B0s = np.logspace(10, 16, 7)
text_X = 1e1
for B0 in B0s:
    power = np.log10(B0)
    if power != int(power):
        print "pow = {} is not an int".format(power)
    ax.text(text_X, 0.8*Pdot(text_X, B0), 
            "$B_{{0}}=$" + "$10^{{{}}}$".format(int(power)),
            rotation=-35.3,            
            bbox=dict(facecolor='w', alpha=0.8, edgecolor="w")
            )
    ax.plot(P_range, Pdot(P_range, B0), ls="--", color="k")

# Add lines of constant characteristic age
def Pdot(P, tau):
    return 0.5 * P/tau 
tau_years = np.logspace(1, 9, 9)
text_X = 2e-4
for tau_year in tau_years:
    tau = tau_year * (365. * 86400)
    rotation=28
    power = np.log10(tau_year)
    if power != int(power):
        print "pow = {} is not an int".format(power)
    ax.text(text_X, 2.5*Pdot(text_X, tau), 
            r"$\tau=$" + "$10^{{{}}}$".format(int(power)),
            rotation=rotation,
            bbox=dict(facecolor='w', alpha=0.8, edgecolor="w")
            )
    ax.plot(P_range, Pdot(P_range, tau), ls="--", color="k")

# Add the DEATH LINE
#CONST = 1e-16
#ax.plot(P_range, CONST * P_range, ls="-", color="r")

#ax.legend(loc=2, frameon=False)
ax.set_xlabel("Period [s]")
ax.set_ylabel("Period deviative")

ax.set_xlim(*xlims)
ax.set_ylim(1e-21, 1e-7)

plt.savefig("Period_PeriodDot.pdf")
plt.show()
