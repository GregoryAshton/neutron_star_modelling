#!/usr/bin/python 

from nsmod.switching_torque_with_Euler import main
from nsmod import Plot
import matplotlib.pyplot as plt
import numpy as np
import argparse



def CompareSpherical(epsI1=0.0, epsI3=3.0e-3, epsA=5.0e-4 , omega0=10,
                    error=1e-12, T=1.0e3 , chi0 = 70.0, AnomTorque=True,
                    a0=20.0, upsilon=0.8, n=10000, cleanup=False):
    """ Compare the body frame components omega, theta and phi for a switched
    and unswitched model
    """

    switching, noswitch = get_data(**locals()) 
    axes = Plot.Spherical_Plot(switching,
                               label="$\upsilon={:1.1f}$".format(upsilon))
    axes = Plot.Spherical_Plot(noswitch, axes=axes, 
                             label=r"$\upsilon=0$")
    plt.show()

def CompareEuler(epsI1=0.0, epsI3=3.0e-3, epsA=5.0e-4 , omega0=10,
                    error=1e-12, T=1.0e3 , chi0 = 70.0, AnomTorque=True,
                    a0=20.0, upsilon=0.8, n=10000, cleanup=False):
    """ Compare the interial frame components Euler angles 
    """

    switching, noswitch = get_data(**locals()) 
    fig, axes = plt.subplots(nrows=3)
    axes = Plot.Euler_Angles(switching, axes=axes, 
                             label="$\upsilon={:1.1f}$".format(upsilon))
    axes = Plot.Euler_Angles(noswitch, axes=axes, 
                             label=r"$\upsilon=0$")
    axes[1].legend(loc=4)
    plt.show()

def CompareResiduals(epsI1=0.0, epsI3=3.0e-3, epsA=5.0e-4 , omega0=10,
                    error=1e-12, T=1.0e3 , chi0 = 70.0, AnomTorque=True,
                    a0=20.0, upsilon=0.8, n=10000, cleanup=False, 
                    order=3):
    """ Compare the timing residual plots """
    main_args = locals()
    del main_args['order'] 
    switching, noswitch = get_data(**main_args) 
    ax = Plot.timing_residual(switching, order=order, 
                              label="$\upsilon={:1.1f}$".format(upsilon))
    ax = Plot.timing_residual(noswitch, order=order, ax=ax, 
                              label="$\upsilon=0.0$")

    plt.legend()
    plt.show()

def CompareNudot(epsI1=0.0, epsI3=2.0e-4, epsA=1.0e-5 , omega0=1e3,
                    error=1e-15, T=100 , chi0=20.0, AnomTorque=True,
                    a0=0.1, upsilon=0.5, n=500000, cleanup=False,
                    divisor=10):
    """ Compare the timing residual plots """
    main_args = locals()
    del main_args['divisor'] 
    switching, noswitch = get_data(**main_args) 
    ax = Plot.nu_dot(switching, divisor=divisor,
                              label="$\upsilon={:1.1f}$".format(upsilon))
    Plot.nu_dot(noswitch, divisor=divisor,  ax=ax,
                              label="$\upsilon=0.0$")

    plt.legend()
    plt.show()
    
def CompareAmplitude(epsI1=0.0, epsI3=3.0e-4, epsA=1e-2 , omega0=1,
                    error=1e-12, T=2e4 , chi0 = 70.0, AnomTorque=True,
                    a0=10.0, upsilon=0.8, n=500000, cleanup=False,
                    Phi0=180, Theta0=50, sigmaPhi=0.3, sigmaTheta=0.3,               
                    ):
    main_args = locals()
    for arg in ['Phi0', 'Theta0', 'sigmaPhi', 'sigmaTheta']:
        del main_args[arg]

    Phi0 = np.radians(Phi0)
    Theta0 = np.radians(Theta0)
    switching, noswitch = get_data(**main_args) 
    axes = Plot.Euler_Angles(switching)
    ax = Plot.Amplitude(switching, Phi0, Theta0, sigmaPhi, sigmaTheta,
                              label="$\upsilon={:1.1f}$".format(upsilon))
    ax = Plot.Amplitude(noswitch, Phi0, Theta0, sigmaPhi, sigmaTheta, ax=ax, 
                              label="$\upsilon=0.0$")
    plt.legend()
    plt.show()

def ComparePulseWidth(epsI1=0.0, epsI3=3.0e-4, epsA=1e-5 , omega0=1,
                    error=1e-12, T=1e5 , chi0 = 70.0, AnomTorque=True,
                    a0=10.0, upsilon=0.8, n=100000, cleanup=False,
                    Phi0=180, Theta0=50, sigmaPhi=0.3, sigmaTheta=0.3, eta=0.01
                    ):
    main_args = locals()
    for arg in ['Phi0', 'Theta0', 'sigmaPhi', 'sigmaTheta', 'eta']:
        del main_args[arg]

    Phi0 = np.radians(Phi0)
    Theta0 = np.radians(Theta0)
    switching, noswitch = get_data(**main_args) 

    ax = Plot.PulseWidth(switching, Phi0, Theta0, sigmaPhi, sigmaTheta,  
                          eta=eta, label="$\upsilon={:1.1f}$".format(upsilon),
                          ls="-",)
    ax = Plot.PulseWidth(noswitch, Phi0, Theta0, sigmaPhi, sigmaTheta, ax=ax, 
                          eta=eta, label="$\upsilon=0.0$", ls="-")
    plt.legend()
    plt.show()

def get_data(epsI1=0.0, epsI3=3.0e-3, epsA=5.0e-4 , omega0=10,
                error=1e-12, T=1.0e3 , chi0 = 70.0, AnomTorque=True,
                a0=20.0, upsilon=0.8, n=10000, cleanup=False):
    args = locals()
    switching = main(**args)
    args['upsilon'] = 0.0
    noswitch = main(**args)
    return switching, noswitch

def _setupArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--spherical",
                        action="store_true", help=CompareSpherical.__doc__)
    parser.add_argument("-e", "--euler", 
                        action="store_true", help=CompareEuler.__doc__)
    parser.add_argument("-r", "--timingresidual", 
                        action="store_true", help=CompareResiduals.__doc__)
    parser.add_argument("-n", "--nudot", 
                        action="store_true", help=CompareNudot.__doc__)
    parser.add_argument("-a", "--amplitude", 
                        action="store_true", help=CompareAmplitude.__doc__)
    parser.add_argument("-p", "--PulseWidth", 
                        action="store_true", help=ComparePulseWidth.__doc__)
    return parser.parse_args()

if __name__ == "__main__":
    args = _setupArgs()
    if args.spherical:
        CompareSpherical()
    if args.euler:
        CompareEuler()
    if args.timingresidual:
        CompareResiduals()
    if args.nudot:
        CompareNudot()
    if args.amplitude:
        CompareAmplitude()
    if args.PulseWidth:
        ComparePulseWidth()
