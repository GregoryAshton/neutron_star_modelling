import numpy as np
from matplotlib import pyplot as plt
import sys
import os
from numpy import cos, sin
import pandas as pd
import re

from nsmod.manual_switching_torque_with_Euler import main
from nsmod import File_Functions, Physics_Functions, Plot


def SignalModel_EM2(params, t):
    omega0, epsI, a0, chi, epsA = params

    theta = a0

    c = 3e10
    R = 1e6
    k = 2/3.0 * R/c * epsA

    psidot = -epsI*omega0
    psi = psidot*t + np.pi/2
    C = 1 - (cos(theta)*cos(chi))**2 - 0.5*(sin(theta)*sin(chi))**2
    tauP = 1./(epsI*omega0)
    tauS = 1./(k*epsA*omega0**2)

    I1 = tauP * cos(t/tauP) + .5*sin(chi)**2/(tauS*tauP) * (
                cos(t/tauP)*t**3/3.0 - tauP * sin(t/tauP)*t**2)
    I2 = -.5*tauP *sin(2*t/tauP) + sin(chi)**2/(tauS*tauP) * (
                sin(2*t/tauP)*t**3/3.0 + .5*tauP * cos(2*t/tauP)*t**2)

    Phidot3 = (2*k*(C*t + (.5*sin(2*theta)*sin(2*chi)*I1
                                      - .25*(sin(theta)*sin(chi))**2*I2))
               + omega0**-2
               )**(-3./2)
    Sin2Theta = 1 - (sin(theta)*sin(psi)*sin(chi) + cos(theta)*cos(chi))**2

    theta = a0
    GEOMETRIC = psidot**2 * ((2*sin(chi)**3*sin(psi)*sin(theta)*cos(theta) -
                              sin(chi)**2*sin(psi)**2*sin(theta)**2*cos(chi) -
                              2*sin(chi)**2*sin(theta)**2*cos(chi) +
                              sin(chi)**2*cos(chi) - sin(theta)**2*cos(chi)**3
                             )*sin(chi)*sin(theta)*cos(psi)/(
                            (sin(chi)*sin(psi)*cos(theta) - sin(theta)*cos(chi))**2 +
                             sin(chi)**2*cos(psi)**2)**2
                            )/ (2*np.pi)

    return  -k * Phidot3 * Sin2Theta / (2*np.pi) + 0*GEOMETRIC


def SignalModel(params, t):
    omega0, epsI, a0, chi, epsA = params

    theta = a0

    c = 3e10
    R = 1e6
    k = 2/3.0 * R/c * epsA
    tauP = 1./(epsI*omega0)

    psi = -epsI*omega0*t + np.pi/2 + 0.5*k*epsI*sin(chi)**2*omega0**3*t**2

    C = 1 - (cos(theta)*cos(chi))**2 - 0.5*(sin(theta)*sin(chi))**2

    T1 = -k * omega0**3 * C
    T2 =  3*k**2*omega0**5 * C**2 * t

    T3 = 0.5 * k * omega0**3 *(sin(2*theta)*sin(2*chi)*sin(psi) -
                               (sin(theta)*sin(chi))**2*cos(2*psi))

    #T4 = 1.5*k**2*omega0**5*C*(-sin(2*theta)*sin(2*chi)*(tauP*cos(psi) + t*sin(psi))
    #                           +(sin(theta)*sin(chi))**2*(.5*tauP*sin(2*psi) + t*cos(2*psi)))

    #T5 = 0.75*k**2*omega0**5*tauP*(
    #          (sin(2*theta)*sin(2*chi)*cos(psi)
    #           - 0.5*((sin(theta)*sin(chi))**2*sin(2*psi))) *
    #          (sin(2*theta)*sin(2*chi)*sin(psi)
    #           - ((sin(theta)*sin(chi))**2*cos(2*psi))))

    return  (T1 + T2 + T3)/(2*np.pi)

def SignalModelWithGeometric(params, t):
    omega0, epsI, a0, chi, epsA = params

    theta = a0

    c = 3e10
    R = 1e6
    k = 2/3.0 * R/c * epsA

    psi = -epsI*omega0*t + np.pi/2 + 0.5*k*epsI*sin(chi)**2*omega0**3*t**2

    C = 1 - (cos(theta)*cos(chi))**2 - 0.5*(sin(theta)*sin(chi))**2

    T1 = -k * omega0**3 * C
    T2 =  3*k**2*omega0**5 * C**2 * t

    T3 = 0.5 * k * omega0**3 *(sin(2*theta)*sin(2*chi)*sin(psi) -
                               (sin(theta)*sin(chi))**2*cos(2*psi))

    theta = a0
    psidot = -epsI*omega0
    GEOMETRIC = psidot**2 * ((2*sin(chi)**3*sin(psi)*sin(theta)*cos(theta) -
                              sin(chi)**2*sin(psi)**2*sin(theta)**2*cos(chi) -
                              2*sin(chi)**2*sin(theta)**2*cos(chi) +
                              sin(chi)**2*cos(chi) - sin(theta)**2*cos(chi)**3
                             )*sin(chi)*sin(theta)*cos(psi)/(
                            (sin(chi)*sin(psi)*cos(theta) - sin(theta)*cos(chi))**2 +
                             sin(chi)**2*cos(psi)**2)**2
                            )
    return  (T1 + T2 + T3 + GEOMETRIC)/(2*np.pi)

dat_files = [string for string in sys.argv if re.match("^.*\.dat", string)]
if len(dat_files) == 1:
    results_file = dat_files[0]
elif len(dat_files) > 1:
    print "Too many dat files.. using default"
elif len(dat_files) == 0:
    results_file = "ResultsBreakingAnalytic.dat"

print "Using results file {}".format(results_file)

# Parameters
omega0 = 100
epsI3 = 5e-6
chi0 = 88.0
n = 100000
error = 1e-16
divisor=9

tauP = 2.0*np.pi/abs(epsI3 * omega0)
T = 3.3 * tauP
c = 3e10
R = 1e6

if "data" in sys.argv:
    if not os.path.isfile(results_file):
        with open(results_file, "w+") as file:
            file.write("T Aem epsI3 epsA omega0 a0 chi res resWG\n")

    D = 25
    for a0 in np.linspace(0.1, 6, D):
        for epsA in np.logspace(-10, -5, D):
            file_name = main(chi0=chi0, epsI3=epsI3, epsA=epsA, omega0=omega0, T=T,
                         n=n, error=error, a0=a0, cleanup=False, DryRun=False,
                         AnomTorque=False)

            #File_Functions.PrintParameterDictionary(file_name)
            out_EA = File_Functions.Euler_Angles_Import(file_name, time_offset=None)
            [time, w1, w2, w3, theta, phi, psi] = out_EA

            time, nu_dot = Physics_Functions.nu_dot_Lyne(time, w1, w2, w3, theta, phi, psi,
                                                    np.radians(chi0), tauP, divisor=divisor)

            theta_EM = np.array([omega0, epsI3, np.radians(a0), np.radians(chi0), epsA])

            nu_dot_analytic = SignalModel(theta_EM, time)
            res = np.sum((nu_dot_analytic - nu_dot)**2) / np.mean(nu_dot_analytic)**2

            nu_dot_analyticWG = SignalModelWithGeometric(theta_EM, time)
            resWG = np.sum((nu_dot_analyticWG - nu_dot)**2) / np.mean(nu_dot_analyticWG)**2

            Aem = 2*R*epsA*(omega0/(2*np.pi))/(3 * c* epsI3**2) * (2*np.pi)**2

            results = "{:1.5e} {:1.5e} {:1.5e} {:1.5e} {:1.5e} {:1.5e} {:1.5e} {:1.5e} {:1.5e}\n".format(
                       T, Aem, epsI3, epsA, omega0, a0, chi0, res, resWG)
            with open("ResultsBreakingAnalytic.dat", "a") as file:
                file.write(results)

                #os.remove(file_name)

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
    Z = z.reshape(X.shape).T

    return X, Y, Z

if "plot" in sys.argv:
    if "WG" in sys.argv:
        zlabel = "resWG"
    else:
        zlabel = "res"

    print zlabel
    vmin, vmax = -11, 2

    X, Y, Z = ReadData("a0", "tauP_over_tauS", zlabel)

    if np.min(Z) < vmin or np.max(Z) > vmax:
        print("vmin {} and vmax {} not covering Z range {} -> {}".format(
              vmin, vmax, np.min(Z), np.max(Z)))


    ax = plt.subplot(111)
    #ax.set_ylabel(r"$\epsilon_A$", rotation='horizontal')
    ax.set_ylabel(r"$\log_{10}\left(\frac{\tau_{\mathrm{P}}}{\tau_{\mathrm{S}}}\right)$",
                  rotation='horizontal', labelpad=40)
    ax.set_xlabel("$a_{0}$ [degrees]")
    pcm = ax.pcolormesh(X, Y, Z, vmin=vmin, vmax=vmax)
    plt.colorbar(pcm, label="$\log_{10}(\mathrm{residual})$",
                 orientation='horizontal',
                 ticks=np.linspace(vmin, vmax, abs(vmax-vmin)+1))

    ax2 = ax.twinx()
    X, Y, Z = ReadData("a0", "Aem", zlabel)
    ax2.pcolormesh(X, Y, Z, vmin=vmin, vmax=vmax)
    ax2.set_ylabel(r"$\mathcal{A}_{\textrm{EM}}$", rotation="horizontal",
                   labelpad=30)

    ax2.set_xlim(np.min(X), np.max(X))
    ax2.set_ylim(np.min(Y), np.max(Y))


    plt.tight_layout()
    plt.savefig("Approximate_nu_dot_validation_{}.pdf".format(zlabel))
    plt.show()

if "corners" in sys.argv:
    if "WG" in sys.argv:
        zlabel = "resWG"
    else:
        zlabel = "res"

    plt.rcParams['lines.linewidth'] = 2.0

    df = pd.read_csv(results_file, delim_whitespace=True)
    a0 = df['a0']
    epsA = df['epsA']

    epsI3 = df['epsI3'][0]
    chi0 = df['chi'][0]
    omega0 = df['omega0'][0]
    #T = np.unique(df['T'])[0]
    #T = 2.07e4 # NEEDS FIXING

    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(14, 7))
    for i, a0 in enumerate([np.min(df['a0']), np.max(df['a0'])]):
        for j, epsA in enumerate([np.max(df['epsA']), np.min(df['epsA'])]):
            file_name = main(chi0=chi0, epsI3=epsI3, epsA=epsA, omega0=omega0,
                             T=T, n=n, error=error, a0=a0, cleanup=False,
                             DryRun=True, AnomTorque=False)
            try:
                ax = Plot.SpindownRate(file_name, ax=axes[j][i], 
                                       label="Numeric solution")
            except IOError:
                break

            out_EA = File_Functions.Euler_Angles_Import(file_name)
            [time, w1, w2, w3, theta, phi, psi] = out_EA


            theta_EM = np.array([omega0, epsI3, np.radians(a0),
                                 np.radians(chi0), epsA])

            if "WG" not in sys.argv:
                nu_dot_analytic = SignalModel(theta_EM, time)
                ax.plot(time, nu_dot_analytic, "-", color="b",
                        label="Analytic without geometric")

            else:
                nu_dot_analyticGeometric = SignalModelWithGeometric(theta_EM,
                                                                time)
                ax.plot(time, nu_dot_analyticGeometric, "--", color="r",
                    label="Analytic with geometric")

            if "verbose" in sys.argv:
                ax.annotate("{} {}".format(a0, epsA), (0.5, 0.5), 
                             xycoords="axes fraction")
                print j, i
                print df[(df.a0 == a0) & (df.epsA == epsA)]

    axes[1][1].legend(bbox_to_anchor=(0.8, 2.46), frameon=False,
                      fontsize=14, ncol=3)
    axes[0][0].set_xlabel("")
    axes[0][1].set_xlabel("")
    #fig.tight_layout()
    plt.savefig("CornersPlot_{}.pdf".format(zlabel))
    plt.show()
