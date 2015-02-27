import matplotlib.pyplot as plt
import numpy as np
import sys

from nsmod.one_component_model import main
from nsmod.one_component_model_with_Euler import main as mainEuler
from nsmod import Plot, File_Functions, Physics_Functions

plt.rcParams['lines.linewidth'] = 2

# Parameters
omega0 = 2*np.pi*10
epsA = 3e-7
epsI3 = 1.05e-6
chi0 = 88.0
a0 = 3.0
n = 50000
error = 1e-13

tauP = 2 * np.pi/abs(epsI3 * omega0)
T = 4.4 * tauP

if "omegaz" in sys.argv:
    # Generate and plot the numerical result
    fig, ax = plt.subplots(figsize=(7, 4))

    FullSoln = main(chi0=chi0, epsI3=epsI3, epsA=epsA, omega0=omega0, T=T,
                     n=n, error=error, a0=a0, AnomTorque=True, cleanup=False)
    (time, wx, wy, wz) = File_Functions.One_Component_Import(FullSoln)
    ax.plot(time, wz, ls="-", color="k", label="Full numerical solution")

    NoAnom = main(chi0=chi0, epsI3=epsI3, epsA=epsA, omega0=omega0, T=T,
                     n=n, error=error, a0=a0, AnomTorque=False, cleanup=False)
    (time, wx, wy, wz) = File_Functions.One_Component_Import(NoAnom)
    ax.plot(time, wz, ls="--", color="k", 
            label="Numerical soln. without anom. torque")


    # Plot the analytic result
    def wz_analytic(t, epsA, omega0, a0):
        R=1e6
        c=3e10
        wz0 = omega0*np.cos(np.radians(a0))
        return 1.0 / np.sqrt(4*R/(3*c) * epsA * t + wz0**-2)
    ax.plot(time, wz_analytic(time, epsA, omega0, a0), ls="-", color="r",
            label="Analytic", zorder=-10)


    ax.set_ylabel(r"$\omega_{z} $", fontsize=20, labelpad=20, rotation="horizontal")
    ax.set_xlabel(r"time")
    ax.legend(framealpha=0.6, fontsize=14, frameon=False)

    fig.tight_layout()

    plt.savefig("img/omegaz.pdf")

    File_Functions.PropertiesTable(FullSoln, 'omegaz', 
                   include=['omega0', 'chi', 'a0', 'tauP', 'tauA', 'tauS'])
    plt.show()

if "psi" in sys.argv:
    fig, ax = plt.subplots(figsize=(7, 4))

    NoAnom = mainEuler(chi0=chi0, epsI3=epsI3, epsA=epsA, omega0=omega0, T=T,
                     n=n, error=error, a0=a0, AnomTorque=False, cleanup=False)
    (time, wx, wy, wz, theta, phi, psi) = File_Functions.Euler_Angles_Import(NoAnom)
    ax.plot(time, psi, ls="--", color="k", label="Numerical caln. without anom. torque")

    # Analytic plot
    def psi_analytic(t, epsI, epsA, a0, omega0):
        R = 1e6
        c = 3e10
        wz0 = omega0 * np.cos(np.radians(a0))
        return -epsI * wz0 * t + R/(3*c) * epsA * epsI * wz0**3 * t**2 + np.pi/2
    ax.plot(time, psi_analytic(time, epsI3, epsA, a0, omega0),
            ls="-", color="r", label="Analytic soln.", zorder=-10)

    ax.set_xlabel(r"time")
    ax.set_ylabel(r"$\psi(t)$", rotation="horizontal", labelpad=20)
    plt.tight_layout()
    plt.savefig("img/psi.pdf")
    plt.show()


