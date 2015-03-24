import numpy as np
from matplotlib import pyplot as plt
import sys
import os
import pandas as pd
import re

from nsmod.manual_switching_torque_with_Euler import main
from nsmod import File_Functions, Physics_Functions, Plot
from nsmod.Spindown_SignalModels import SignalModelWithGeometric

def diff(x, y):
    #return np.abs(np.sqrt(np.sum((x-y)**2)) / np.mean(x))
    return np.abs(np.sqrt(np.mean((x-y)**2)) / np.mean(x))

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
n = 500000
error = 1e-17
divisor=7

tauP = 2.0*np.pi/abs(epsI3 * omega0)
T = 3.3 * tauP
c = 3e10
R = 1e6

print "Summary of fixed quantities"
print "omega0 : {}\ntauP : {:1.2e}\nchi0 : {}".format(omega0, tauP, chi0)

if "data" in sys.argv:
    if not os.path.isfile(results_file):
        with open(results_file, "w+") as file:
            file.write("T Aem epsI3 epsA omega0 a0 chi res resWG\n")

    D = 30
    for a0 in np.linspace(0.1, 6, D):
        for epsA in np.logspace(-10, -6, D):
            file_name = main(chi0=chi0, epsI3=epsI3, epsA=epsA, omega0=omega0, T=T,
                         n=n, error=error, a0=a0, cleanup=False, DryRun=False,
                         AnomTorque=False)

            #File_Functions.PrintParameterDictionary(file_name)
            out_EA = File_Functions.Euler_Angles_Import(file_name, time_offset=None)
            [time, w1, w2, w3, theta, phi, psi] = out_EA

            time, nu_dot = Physics_Functions.nu_dot_Lyne(time, w1, w2, w3, theta, phi, psi,
                                                    np.radians(chi0), tauP, divisor=divisor)

            theta_EM = np.array([omega0, epsI3, np.radians(a0), np.radians(chi0), epsA])

            nu_dot_analytic = SignalModelWithGeometric(theta_EM, time, geometric=0)
            res = diff(nu_dot, nu_dot_analytic)
            #res = np.sum((nu_dot_analytic - nu_dot)**2) / np.mean(nu_dot_analytic**2)

            nu_dot_analyticWG = SignalModelWithGeometric(theta_EM, time, geometric=1)
            resWG = diff(nu_dot, nu_dot_analyticWG)
            #resWG = np.sum((nu_dot_analyticWG - nu_dot)**2) / np.mean(nu_dot_analyticWG**2)

            Aem = 2*R*epsA*(omega0/(2*np.pi))/(3 * c* epsI3**2) * (2*np.pi)**2

            results = "{:1.5e} {:1.5e} {:2e} {:.2e} {:1.5e} {:1.10e} {:1.10e} {:1.5e} {:1.5e}\n".format(
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
    vmin, vmax = -6, 0

    X, Y, Z = ReadData("a0", "tauP_over_tauS", zlabel)

    print("Using limits {}, {} for data with range {}, {}".format(
          vmin, vmax, np.min(Z), np.max(Z)))


    ax = plt.subplot(111)
    #ax.set_ylabel(r"$\epsilon_A$", rotation='horizontal')
    ax.set_ylabel(r"$\log_{10}\left(\frac{\tau_{\mathrm{P}}}{\tau_{\mathrm{S}}}\right)$",
                  rotation='horizontal', labelpad=40)
    ax.set_xlabel(r"$\theta$ [degrees]")
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
                                       divisor=divisor,
                                       label="Numeric solution")
            except IOError:
                break
            print a0, epsA
            out_EA = File_Functions.Euler_Angles_Import(file_name)
            [time, w1, w2, w3, theta, phi, psi] = out_EA

            time, nu_dot = Physics_Functions.nu_dot_Lyne(time, w1, w2, w3, theta, phi, psi,
                                                    np.radians(chi0), tauP, divisor=divisor)

            time_fine = np.linspace(time[0], time[-1], 1000)

            theta_EM = np.array([omega0, epsI3, np.radians(a0),
                                 np.radians(chi0), epsA])

            if "WG" not in sys.argv:
                nu_dot_analytic = SignalModelWithGeometric(theta_EM, time_fine, geometric=0)
                line, = ax.plot(time_fine, nu_dot_analytic, "--", color="r",
                        label="Analytic without geometric")

            else:
                nu_dot_analyticWG = SignalModelWithGeometric(theta_EM, time_fine)
                line, = ax.plot(time_fine, nu_dot_analyticWG, "--", color="r",
                    label="Analytic with geometric")

            line.set_dashes((3,2))

            if "verbose" in sys.argv:
                print j, i
                line = df[(df.a0 == a0) & (df.epsA == epsA)]
                print line
                ax.annotate("{} \n{} \n{}".format(a0, epsA, line[zlabel].values[0]),
                            (0.2, 0.75), 
                             xycoords="axes fraction")



    axes[1][1].legend(bbox_to_anchor=(0.5, 2.46), frameon=False,
                      fontsize=14, ncol=3)
    axes[0][0].set_xlabel("")
    axes[0][1].set_xlabel("")

    for ax, label in zip(axes.flat, map(chr, range(65, 69))):
        ax.annotate(label, (0.008, 0.9), xycoords="axes fraction", size=22,
                    bbox=dict(facecolor='white', alpha=0.9, edgecolor='white'))
    #fig.tight_layout()
    plt.savefig("CornersPlot_{}.pdf".format(zlabel))
    plt.show()

def get_idxs(df, xlabel, ylabel, n):
    nx = len(df[xlabel])
    ny = len(df[ylabel])
    xidxs = np.arange(0, nx+1, nx/float(n-1) -1).astype(int)
    yidxs = np.arange(0, ny+1, ny/float(n-1) -1).astype(int)
    return xidxs, yidxs

if "points" in sys.argv:
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

    nrows = ncols = 3
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(14, 7))
    a0idxs, epsAidxs = get_idxs(df, 'a0', 'epsA', nrows)
    print "lookhere" 
    print a0idxs, df['a0'][a0idxs].values
    print epsAidxs, df['epsA'][epsAidxs].values
    print
    for i, a0 in enumerate(df['a0'][a0idxs].values):
        for j, epsA in enumerate(df['epsA'][epsAidxs].values):
            file_name = main(chi0=chi0, epsI3=epsI3, epsA=epsA, omega0=omega0,
                             T=T, n=n, error=error, a0=a0, cleanup=False,
                             DryRun=True, AnomTorque=False)
            try:
                ax = Plot.SpindownRate(file_name, ax=axes[nrows-1-j][i], 
                                       divisor=divisor,
                                       label="Numeric solution")
            except IOError:
                print file_name
                print "WARNING: file not found"
                continue
            out_EA = File_Functions.Euler_Angles_Import(file_name)
            [time, w1, w2, w3, theta, phi, psi] = out_EA

            time, nu_dot = Physics_Functions.nu_dot_Lyne(time, w1, w2, w3, theta, phi, psi,
                                                    np.radians(chi0), tauP, divisor=divisor)

            time_fine = np.linspace(time[0], time[-1], 1000)

            theta_EM = np.array([omega0, epsI3, np.radians(a0),
                                 np.radians(chi0), epsA])

            if "WG" not in sys.argv:
                nu_dot_analytic = SignalModelWithGeometric(theta_EM, time_fine, geometric=0)
                line, = ax.plot(time_fine, nu_dot_analytic, "--", color="r",
                        label="Analytic without geometric")

            else:
                nu_dot_analyticWG = SignalModelWithGeometric(theta_EM, time_fine)
                line, = ax.plot(time_fine, nu_dot_analyticWG, "--", color="r",
                    label="Analytic with geometric")

            line.set_dashes((3,2))

            if "verbose" in sys.argv:
                print j, i
                line = df[(df.a0 == a0) & (df.epsA == epsA)]
                print line
                ax.annotate("{} \n{} \n{}".format(a0, epsA, line[zlabel].values[0]),
                            (0.2, 0.75), 
                             xycoords="axes fraction")



    #axes[1][1].legend(bbox_to_anchor=(0.5, 2.46), frameon=False,
    #                  fontsize=14, ncol=3)
    axes[0][0].set_xlabel("")
    axes[0][1].set_xlabel("")

    for ax, label in zip(axes.flat, map(chr, range(65, 65+nrows*ncols))):
        ax.annotate(label, (0.008, 0.9), xycoords="axes fraction", size=22,
                    bbox=dict(facecolor='white', alpha=0.9, edgecolor='white'))
    #fig.tight_layout()
    plt.savefig("CornersPlot_{}.pdf".format(zlabel))
    plt.show()
