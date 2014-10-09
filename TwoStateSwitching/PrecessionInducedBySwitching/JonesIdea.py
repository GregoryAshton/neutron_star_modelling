import matplotlib.pyplot as plt
import numpy as np
from nsmod import Plot
import sys
from nsmod.manual_switching_torque_with_Euler import main
from nsmod import Plot
from nsmod.Physics_Functions import Beta_Function
from nsmod import File_Functions

R = 1e6
c = 3e10

def BetaPrime(epsI, epsA, upsilon, chi):
    chi = np.radians(chi)
    epsA = upsilon * epsA
    top = epsI - epsA * np.cos(2*chi) - np.sqrt(epsA**2 + epsI**2 - 2 * epsA * 
                                                       epsI * np.cos(2 * chi))
    bottom = 2 * epsA * np.sin(chi) * np.cos(chi)

    return np.arctan2(top, bottom)

if "BetaPlot" in sys.argv:

    ax = plt.subplot(111)
    chi = 65
    epsI = 1e-3
    epsA = np.logspace(-5, -1, 100)

    for upsilon, upsilon_label in zip([0.2, 0.5, 0.8], 
                                  [r"$\frac{1}{5}$", r"$\frac{1}{2}$", r"$\frac{4}{5}$"]):
        DeltaBeta = np.abs(BetaPrime(epsI, epsA, 1, chi) - 
                           BetaPrime(epsI, epsA, upsilon, chi))

        ax.semilogx(epsA/epsI, np.degrees(DeltaBeta), 
                    label="$\upsilon=$" + upsilon_label, rasterized=True)
        ax.set_xlabel(r"$\frac{\epsilon_{A}}{\epsilon_{I}}$")
        ax.set_ylabel(r"$\Delta \beta = |\beta' - \beta|$ [degrees]", rotation="vertical")

    ax.legend(loc=1, frameon=False, fontsize=20)
    plt.tight_layout()
    plt.savefig("img/DeltaBetaPlot.pdf", dpi=270)
    plt.show()

if "SingleSwitch" in sys.argv:
    epsI3 = 6e-6
    epsA = 1e-8
    chi0 = 40.0
    Beta = Beta_Function(epsI3, epsA, chi0)
    a0 = np.degrees(Beta)
    upsilon = 0.1
    omega0 = 1e4
    T = 2.0e3
    file_name = main(epsI1=0.0, epsI3=epsI3, epsA=epsA, omega0=omega0, chi0=chi0,
    a0=a0, T=T, AnomTorque=True , upsilon=upsilon, n=500000, error=1e-10,
    cleanup=False, SwitchTime=1.0e3
    )

    File_Functions.PrintParameterDictionary(file_name)

    Plot.Spherical_Plot(file_name)
    plt.show()
      
    Plot.nu_dot(file_name)
    plt.show()

    Plot.timing_residual(file_name, order=3)
    plt.show()

    print "Size of variations"
    theta = BetaPrime(epsI3, epsA, 1, chi0) - BetaPrime(epsI3, epsA, upsilon, chi0)
    print "dPhi_FP = {:1.2e}".format(theta * np.cos(np.radians(chi0)) / np.sin(np.radians(chi0)))
    F1 = -(2/3.0) * (R * omega0**3 / c) * epsA * np.sin(np.radians(chi-a0))**2
    print "dPhi_TS = {:1.2e}".format(upsilon * F1 * T**2 * np.pi/8.)

if "B1828" in sys.argv:
    epsI = np.logspace(-12, -10, 100)
    chi = 30
    upsilon = 0.71

    T = (1/0.73) * 86400 * 365.25
    nu0 = 2.47
    nudot0 = -365.e-15
    omega0 = 2 * np.pi * nu0 
    omegadot0 = 2 * np.pi * nudot0
    epsA = -omegadot0 * c /(((2/3.) * R * omega0**3))

    DeltaBeta = np.degrees(np.abs(BetaPrime(epsI, epsA, 1, chi) - BetaPrime(epsI, epsA, upsilon, chi)))
    plt.semilogx(epsA/epsI, DeltaBeta)
    plt.show()

    print "DeltaBeta = {:1.2e}".format(np.max(DeltaBeta))
    print "EpsilonA = {}".format(epsA)
    print "dPhi_TS = {:1.2f}".format(upsilon * nudot0 * T**2 * np.pi/16.)
 
