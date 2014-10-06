""" """

import numpy as np
import matplotlib.pyplot as plt
from nsmod.one_component_model_with_Euler import main
from nsmod import Plot
from nsmod import File_Functions
from nsmod import Physics_Functions
import matplotlib.pyplot as plt
import pandas as pd
import os, sys

def PlotNuDot(epsI1=0.0, epsI3=1.0e-6, epsA=1.0e-9 , omega0=10000,
                    error=1e-10, T=5.0e3 , chi0=10.0, 
                    a0=30.0, n=10000, cleanup=False,
                    divisor=10):
    """ Compare the timing residual plots """
    
    file_name = main(epsI1=epsI1, epsI3=epsI3, epsA=epsA, omega0=omega0,
                    error=error, T=T, chi0=chi0, a0=a0, n=n, cleanup=cleanup)

    for key, val in File_Functions.Parameter_Dictionary(file_name).iteritems():
        try:
            formatted_val = "{:1.4e}".format(float(val))
            print key, ":", formatted_val 
        except ValueError:
            print key, ":", val

    ax = Plot.nu_dot(file_name)
    
    out = File_Functions.Euler_Angles_Import(file_name)
    [time, w1, w2, w3, theta, phi, psi] = out
    chi0 = np.radians(File_Functions.Parameter_Dictionary(file_name)['chi0'])

    order=3
    Tres, coeffs = Physics_Functions.timing_residual(time, w1, w2, w3, 
                                             theta, phi, psi, chi0, 
                                             order=order, full=True)
    [F2, F1, F0, P0] = coeffs
    n = F2 * F0 / F1**2
    print "Braking index: {}".format(n)

    plt.show()

if "plot" in sys.argv:
    PlotNuDot()

def ReadResults(file_name, cleanup=False, add_Tobs=False, R55362=False):
    df = pd.read_csv(file_name, sep=" ", skipinitialspace=True)
    return df

def UpdateResults(file_name, results):
    """ Look for existing file and append new results, if it doesn't 
        exist then create it
    """
    if os.path.isfile(file_name):
        df = ReadResults(file_name)
        df = df.append(results, ignore_index=True)
    else:
        df = pd.DataFrame(results, index=[0])
    df.to_csv(file_name, sep=" ")

def GenerateData(epsI1=0.0, epsI3=5.0e-6, epsA=1.0e-8 , omega0=10000,
                    error=1e-10, T=1.5e3 , chi0=10.0, 
                    a0=30.0, n=10000, cleanup=False):
    """ """
    
    file_name = main(epsI1=epsI1, epsI3=epsI3, epsA=epsA, omega0=omega0,
                    error=error, T=T, chi0=chi0, 
                    a0=a0, n=n, cleanup=cleanup)

    out = File_Functions.Euler_Angles_Import(file_name)
    [time, w1, w2, w3, theta, phi, psi] = out
    chi0 = np.radians(File_Functions.Parameter_Dictionary(file_name)['chi0'])

    order=3
    Tres, coeffs = Physics_Functions.timing_residual(time, w1, w2, w3, 
                                             theta, phi, psi, chi0, 
                                             order=order, full=True)
    return coeffs

def Secs2Year(sec):
    return sec / (86400 * 365.25)

def MagneticField(epsA):
    I0 = 1e45
    c = 3e10
    R = 1e6
    return np.sqrt(4 * epsA * I0 * c**2 / R**5)

file_name = "MonteCarloBrakingIndexResults.txt"

if "data" in sys.argv:
    N = 1000
    T = 5.0e3
    omega0 = 9500
    error = 1e-10
    a0 = 30.0
    #for i in range(N):
    #for chi0 in np.linspace(0, 180, 20):
    for epsA in np.logspace(-8, -11, 20):
        epsI3 = 1.0e-6
        #epsA = 1e-10 #, abs(np.random.normal(1e-9, 1e-10))
        u = np.random.uniform(0, 1)
        chi0 = 10 #np.degrees(np.arccos(1 - 2 * u))
        coeffs = GenerateData(epsI3=epsI3, epsA=epsA, chi0=chi0, T=T, 
                              omega0=omega0, a0=a0, error=error)
        results = {'F2' : coeffs[0],
                   'F1' : coeffs[1],
                   'F0' : coeffs[2],
                   'P0' : coeffs[3],
                   'a0' : a0,
                   'omega0' : omega0,
                   'chi0' : chi0,
                   'epsI3' : epsI3,
                   'epsA' : epsA}
        UpdateResults(file_name, results)

if "plotN" in sys.argv:
    df = ReadResults(file_name)
    df['n'] = df.F2 * df.F0 / df.F1**2
    df['tau'] = -.5 * df.F0 / df.F1

    ax = plt.subplot(111)
    ax.scatter(Secs2Year(df.tau), df.n, marker="o", color="k")
    ax.scatter(Secs2Year(df[df['omega0'] == 9000].tau), df[df['omega0'] == 9000].n, marker="o", color="r")
    ax.scatter(Secs2Year(df[df['omega0'] == 9500].tau), df[df['omega0'] == 9500].n, marker="o", color="b")
    ax.set_xscale('log')
    ax.set_yscale('symlog')

    ax.set_xlabel(r"Characteristic age $\tau$ [years]")
    ax.set_ylabel("Observed braking index $n$")

    print df.sort('n').head(5)
    plt.savefig("nobs_tau.pdf")
    plt.show()

    # Plot the distributions
    fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(10, 5))
    ax1.hist(MagneticField(df.epsA), bins=50)
    ax1.set_xlabel("Magnetic field strength ($B_{0}$)", size=14)
    ax1.set_ylabel("Counts")

    ax2.hist(df.chi0, bins=50)
    ax2.set_xlabel("Magnetic inclination angle ($\chi$)", size=14)
    ax2.set_ylabel("Counts")
    plt.savefig("histograms.pdf")
    plt.show()
