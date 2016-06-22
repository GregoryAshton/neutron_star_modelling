import matplotlib.pyplot as plt
import numpy as np
import sys
from nsmod import File_Functions, Physics_Functions
from math import pi
import matplotlib.ticker as ticker

plt.style.use('thesis')

# We need to do this for two distinct sets of data
if "30" in sys.argv:
    f_average = "Averaged_ODEs/averaged_chi_30.0.txt"
    f_exact = "no_anom_chi_30.0_epsI_1.0e-9_epsA_5.0e-11_omega0_1.0e4_eta_1.0e-4.txt"
    f_exact = "../Pulsar_A_No_Anomalous_Torque/data/one-component-model_eta_0.00e+00_chi0_3.0000000000e+01_omega0_1.00e+04_epsI3_1.00e-09_epsA_5.00e-11_a0_5.0000000000e+01_T_1.00e+10_epsI1_0.00e+00_AnomTorque_0.hdf5"
    label_title = r"$\chi=30^{\circ}$"
    label_loc = 2
    maxn = 130000/4

if "75" in sys.argv:
    f_average = "Averaged_ODEs/averaged_chi_75.0.txt"
    f_exact = "../Pulsar_A_No_Anomalous_Torque/data/one-component-model_eta_0.00e+00_chi0_7.5000000000e+01_omega0_1.00e+04_epsI3_1.00e-09_epsA_5.00e-11_a0_5.0000000000e+01_T_1.00e+10_epsI1_0.00e+00_AnomTorque_0.hdf5"
    label_title = r"$\chi=75^{\circ}$"
    label_loc = 2
    maxn = -1

# Import the relevant data and define the axis
fig, ax = plt.subplots(figsize=(3.5, 3))

(time, omega_x, omega_y, omega_z) = File_Functions.One_Component_Import(f_exact, maxn)

# Transform to spherical polar coordiantes
(omega, a, phi) = Physics_Functions.Cartesian_2_Spherical(omega_x, omega_y, omega_z)

# Plot the data scaling the axis as we go
ax.plot(time, a, label="Exact soln.", color="k", alpha=1.0, ls="-", lw=1.0)
ax.set_xlabel(r"time [s]")
ax.set_ylabel(r"$a$ [deg]")

# Plot the averaged solution after solving the Goldreich averaged ODEs

# Manual import of the data
data = open(f_average,"r")
time_average = []
omega_average = []
a_average = []
params = []

for line in data:
    line = line.split()
    if len(params) == 0:
        params = (line)
    else:
        time_average.append(float(line[0]))
        omega_average.append(float(line[1]))
        a_average.append(float(line[2])*180/pi)

d = (3, 1.5)
ax.plot(time_average, a_average, label="Average soln.", lw=1, color="red",
        dashes=d)

if "close" in sys.argv:
    # Inset a plot of the fast oscillations this is not automated
    scale = 1e-9
    xmin = 1
    xmax = 1.3
    if "75" in sys.argv:
        ax2 = plt.axes([.45, .3, .45, .38])
        plt.xlim(xmin, xmax)
        plt.ylim(80.5, 81.5)
        plt.xticks(np.linspace(xmin, xmax, 3))
        ax2.yaxis.set_major_locator(ticker.MaxNLocator(5))
    if "30" in sys.argv:
        ax2 = plt.axes([.4, .35, .45, .38])
        plt.xlim(xmin, xmax)
        plt.ylim(0.87, 1.3)
        plt.xticks(np.linspace(xmin, xmax, 3))
        ax2.yaxis.set_major_locator(ticker.MaxNLocator(5))

    plt.plot(np.array(time)*scale, a, label="Exact soln.", color="k", alpha=1.0,
             ls="-", lw=1.0)
    plt.plot(np.array(time_average)*scale, a_average, color="red", dashes=d,
             lw=1.0)
    for tick in ax2.xaxis.get_major_ticks():
        tick.label.set_fontsize(4)
    for tick in ax2.yaxis.get_major_ticks():
        tick.label.set_fontsize(4)

ax.legend(loc=label_loc, frameon=False)
ax.grid(True, linestyle='-', linewidth=0.1)
ax.set_xlim(0.0, 2e9)


fig.tight_layout()
if '30' in sys.argv:
    plt.savefig('Plot_a_averaged_and_exact_chi_30.png')
else:
    plt.savefig('Plot_a_averaged_and_exact_chi_75.png')

