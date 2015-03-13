import numpy as np
from matplotlib import pyplot as plt
import sys 
import os
from numpy import cos, sin
import pandas as pd

from nsmod.manual_switching_torque_with_Euler import main
from nsmod import File_Functions, Physics_Functions, Plot


#def SignalModel_EM2(theta, t):
#    omega0, epsI, a0, chi, epsA = theta
#    psidot = -epsI*omega0
#    psi = psidot*t + np.pi/2.
#
#    GEOMETRIC = (psidot**2*(2*(sin(chi)*cos(a0) - sin(a0)*sin(psi)*cos(chi))
#                      *(sin(chi)*sin(a0)*sin(psi) + cos(chi)*cos(a0))*sin(chi)
#                      - ((sin(chi)*sin(psi)*cos(a0) - sin(a0)*cos(chi))**2
#                      + sin(chi)**2*cos(psi)**2)*cos(chi))*sin(chi)*sin(a0)
#                      *cos(psi)/((sin(chi)*sin(psi)*cos(a0)
#                    - sin(a0)*cos(chi))**2 
#                      + sin(chi)**2*cos(psi)**2)**2
#                        ) / (2*np.pi) 
#    
#
#    Sin2Theta = 1 - (sin(a0)*sin(psi)*sin(chi) + cos(a0)*cos(chi))**2
#    Phidot = omega0*(1 - epsI * sin(chi) * (cos(a0) *sin(chi) - sin(psi) * sin(a0) * cos(chi)) / (
#                        (sin(a0)*cos(chi) - cos(a0)*sin(psi)*sin(chi))**2 + (cos(psi) * sin(chi))**2))
#
#    # Standard values for c , R in cgs
#    c = 3e10
#    R = 1e6
#    k = 2/3.0 * R/c * epsA 
#    PHIDDOT = -k * Phidot**3 * Sin2Theta / (2*np.pi)
#
#    return PHIDDOT + GEOMETRIC 

#def SignalModel_EM20(theta1, t):
#    omega0, epsI, a0, chi, epsA = theta1
#    
#    c = 3e10
#    R = 1e6
#    k = 2/3.0 * R/c * epsA
#    
#    omega = 1/np.sqrt(2*k*t + 1/omega0**2)
#    
#    psi = np.pi/2 - epsI*omega*t + .5*k*epsI*omega**3*t**2
#    psidot = - epsI*omega + k*epsI*omega**3*t
#    
#    Sin2Theta = 1 - (sin(a0)*sin(psi)*sin(chi) + cos(a0)*cos(chi))**2
#    Phidot = omega - psidot * sin(chi) * (cos(a0) *sin(chi) - sin(psi) * sin(a0) * cos(chi)) / (
#                        (sin(a0)*cos(chi) - cos(a0)*sin(psi)*sin(chi))**2 + (cos(psi) * sin(chi))**2)
#    
#    T1 = -k * Phidot**3 * Sin2Theta / (2*np.pi)
#    
#    theta = a0
#    GEOMETRIC = psidot**2 * ((2*sin(chi)**3*sin(psi)*sin(theta)*cos(theta) - 
#                              sin(chi)**2*sin(psi)**2*sin(theta)**2*cos(chi) - 
#                              2*sin(chi)**2*sin(theta)**2*cos(chi) + 
#                              sin(chi)**2*cos(chi) - sin(theta)**2*cos(chi)**3
#                             )*sin(chi)*sin(theta)*cos(psi)/(
#                            (sin(chi)*sin(psi)*cos(theta) - sin(theta)*cos(chi))**2 + 
#                             sin(chi)**2*cos(psi)**2)**2
#                            )/ (2*np.pi)
#
#    return T1 + GEOMETRIC

def SignalModel_EM(theta, t):
    omega0, epsI, a0, chi, epsA = theta

    c = 3e10
    R = 1e6
    k = (2*R)/(3*c) * epsA

    #omega = omega0 - k*omega0**3*t*sin(chi)**2
    #psi = -epsI*omega0 + 0.5*k*epsI*omega0**3*sin(chi)**2*t**2 + np.pi/2

    omega = (2*k*sin(chi)**2*t + omega0**-2)**(-0.5)
    C =  np.pi/2 + epsI/(k*sin(chi)**2 * omega0)
    psi = -epsI/(k*sin(chi)**2) * (2*k*sin(chi)**2*t + omega0**-2)**(0.5) + C

    Sin2Theta = 1 - (sin(a0)*sin(psi)*sin(chi) + cos(a0)*cos(chi))**2
    Phidot = omega

    return -k * Phidot**3 * Sin2Theta / (2*np.pi)

if len(sys.argv) > 2:
    results_file = sys.argv[2]
else:
    results_file = "ResultsBreakingAnalytic.dat"
print results_file

if "data" in sys.argv:
    # Parameters
    omega0 = 2*np.pi*100
    epsI3 = -2e-6
    chi0 = 85.9
    n = 100000
    error = 1e-11
    divisor=7
    tauP = 2.0 * np.pi/abs(epsI3 * omega0)
    T = 2.4 * tauP
    c = 3e10
    R = 1e6

    if not os.path.isfile(results_file):
        with open(results_file, "w+") as file:
            file.write("Aem epsI3 epsA omega0 a0 chi res\n")

    D = 20
    for a0 in np.linspace(1.0, 10, D):
        for epsA in np.logspace(-10.5, -7.5, D):
            file_name = main(chi0=chi0, epsI3=epsI3, epsA=epsA, omega0=omega0, T=T, 
                         n=n, error=error, a0=a0, cleanup=False, DryRun=False, 
                         AnomTorque=True)

            File_Functions.PrintParameterDictionary(file_name)
            out_EA = File_Functions.Euler_Angles_Import(file_name, time_offset=None)
            [time, w1, w2, w3, theta, phi, psi] = out_EA

            time, nu_dot = Physics_Functions.nu_dot(time, w1, w2, w3, theta, phi, psi, 
                                                    np.radians(chi0), tauP, divisor=divisor)

            theta_EM = np.array([omega0, epsI3, np.radians(a0), np.radians(chi0), epsA])

            nu_dot_analytic = SignalModel_EM(theta_EM, time)

            res = np.sum((nu_dot_analytic - nu_dot)**2)
            Aem = 2*R*epsA*(omega0/(2*np.pi))/(3 * c* epsI3**2) * (2*np.pi)**2
            results = "{:1.5e} {:1.5e} {:1.5e} {:1.5e} {:1.5e} {:1.5e} {:1.5e}\n".format(
                       Aem, epsI3, epsA, omega0, a0, chi0, res)
            with open("ResultsBreakingAnalytic.dat", "a") as file:
                file.write(results)

                os.remove(file_name)

def ReadData(xlabel, ylabel, zlabel):
    df = pd.read_csv(results_file, delim_whitespace=True)
    df['tauP'] = np.abs(1.0 / (df.epsI3 * df.omega0))
    df['tauS'] = 3.*3e10/(2.*1e6 * df.epsA * df.omega0**2) 
    df['tauP_over_tauS'] = df.tauP/df.tauS

    z = np.log10(df[zlabel].values)
    y = np.log10(np.unique(df[ylabel].values))
    x = np.unique(df[xlabel].values)

    if len(x) != len(z):
        dx = x[1] - x[0]
        x = np.linspace(x[0], x[0] + len(y)*dx, len(y))
        z2 = np.zeros(len(x)**2)  + np.average(z)
        z2[:len(z)] = z
        z = z2

    X, Y = np.meshgrid(x, y)
    Z = z.reshape(X.shape)

    return X, Y, Z

if "plot" in sys.argv:
    X, Y, Z = ReadData("a0", "tauP_over_tauS", "res")

    ax = plt.subplot(111)
    ax.set_ylabel(r"$\log_{10}(\tau_P / \tau_S)$")
    ax.set_xlabel("$a_{0}$")
    pcm = ax.pcolormesh(X, Y, Z)
    plt.colorbar(pcm, label="$\log_{10}(\mathrm{residual})$")


    plt.show()

if "Aem" in sys.argv:
    X, Y, Z = ReadData("Aem", "tauS", "res")
    ax = plt.subplot(111)
    pcm = ax.pcolormesh(X, Y, Z)
    plt.colorbar(pcm, label="$\log_{10}(\mathrm{residual})$")
    ax.set_xlabel("$\log_{10}(A_{EM})$")
    ax.set_ylabel(r"$\log_{10}(\tau_P / \tau_S)$")

    plt.show()
