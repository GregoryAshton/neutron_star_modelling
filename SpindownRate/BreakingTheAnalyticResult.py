import numpy as np
from matplotlib import pyplot as plt
import sys 
import os
from numpy import cos, sin
import pandas as pd

from nsmod.manual_switching_torque_with_Euler import main
from nsmod import File_Functions, Physics_Functions

def SignalModel_EM2(theta, t):
    omega0, epsI, a0, chi, epsA = theta
    psidot = -epsI*omega0
    psi = psidot*t + np.pi/2.

    GEOMETRIC = (psidot**2*(2*(sin(chi)*cos(a0) - sin(a0)*sin(psi)*cos(chi))
                      *(sin(chi)*sin(a0)*sin(psi) + cos(chi)*cos(a0))*sin(chi)
                      - ((sin(chi)*sin(psi)*cos(a0) - sin(a0)*cos(chi))**2
                      + sin(chi)**2*cos(psi)**2)*cos(chi))*sin(chi)*sin(a0)
                      *cos(psi)/((sin(chi)*sin(psi)*cos(a0)
                    - sin(a0)*cos(chi))**2 
                      + sin(chi)**2*cos(psi)**2)**2
                        ) / (2*np.pi) 
    

    Sin2Theta = 1 - (sin(a0)*sin(psi)*sin(chi) + cos(a0)*cos(chi))**2
    Phidot = omega0*(1 - epsI * sin(chi) * (cos(a0) *sin(chi) - sin(psi) * sin(a0) * cos(chi)) / (
                        (sin(a0)*cos(chi) - cos(a0)*sin(psi)*sin(chi))**2 + (cos(psi) * sin(chi))**2))

    # Standard values for c , R in cgs
    c = 3e10
    R = 1e6
    k = 2/3.0 * R/c * epsA 
    PHIDDOT = -k * Phidot**3 * Sin2Theta / (2*np.pi)

    return PHIDDOT + GEOMETRIC 

results_file = "ResultsBreakingAnalytic.dat"

if "data" in sys.argv:
    # Parameters
    omega0 = 2*np.pi*1000
    epsI3 = 1.05e-3
    epsA = 3e-10
    chi0 = 87.0
    a0 = 3.0
    n = 50000
    error = 1e-11
    divisor=7
    tauP = 2.0 * np.pi/abs(epsI3 * omega0)
    T = 2.4 * tauP
    c = 3e10
    R = 1e6

    if not os.path.isfile(results_file):
        with open(results_file, "w+") as file:
            file.write("Aem epsI3 epsA omega0 a0 chi res\n")

    for chi0 in np.linspace(45, 89., 50):
        for a0 in np.linspace(1, 45, 50):
            for epsA in np.logspace(-6, -4, 50):
                file_name = main(chi0=chi0, epsI3=epsI3, epsA=epsA, omega0=omega0, T=T, 
                             n=n, error=error, a0=a0, cleanup=False, DryRun=False, 
                             AnomTorque=True)

                out_EA = File_Functions.Euler_Angles_Import(file_name, time_offset=None)
                [time, w1, w2, w3, theta, phi, psi] = out_EA

                time, nu_dot = Physics_Functions.nu_dot(time, w1, w2, w3, theta, phi, psi, 
                                                        np.radians(chi0), tauP, divisor=divisor)

                theta_EM = np.array([omega0, epsI3, np.radians(a0), np.radians(chi0), epsA])
                nu_dot_analytic = SignalModel_EM2(theta_EM, time)

                res = np.sum((nu_dot_analytic - nu_dot)**2)
                Aem = 2*R*epsA*(omega0/(2*np.pi))/(3 * c* epsI3**2) * (2*np.pi)**2
                results = "{:1.5e} {:1.5e} {:1.5e} {:1.5e} {:1.5e} {:1.5e} {:1.5e}\n".format(
                           Aem, epsI3, epsA, omega0, a0, chi0, res)
                with open("ResultsBreakingAnalytic.dat", "a") as file:
                    file.write(results)

                os.remove(file_name)

if "plot" in sys.argv:
    df = pd.read_csv(results_file, delim_whitespace=True)
    df.tauP = 1.0 / (df.epsI3 * df.omega0)
    df.tauS = 3.*3e10/(2.*1e6 * df.epsA * df.omega0**2) 

    z = np.log10(df.res.values)
    x = np.log10(np.unique(df.tauP.values/df.tauS.values))
    y = np.unique(df.a0.values)
    print len(x), len(y), len(z)

    X, Y = np.meshgrid(x, y)
    Z = z.reshape(X.shape)
    ax = plt.subplot(111)
    pcm = ax.pcolormesh(X, Y, Z)
    plt.colorbar(pcm, label="$\log_{10}(\mathrm{residual})$")


    plt.show()
