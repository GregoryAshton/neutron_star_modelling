import matplotlib.pyplot as plt
import numpy as np
from nsmod import Plot
from nsmod import File_Functions
from nsmod.switching_torque_with_Euler import main

def ExamplePlot(epsI1=0.0, epsI3=5.0e-6, epsA=1.0e-8 , omega0=1e4,
                    error=1e-10, T=1.5e3 , chi0=20.0, AnomTorque=True,
                    a0=0.1, upsilon=0.5, n=500000, cleanup=False,
                    divisor=10):
    """ Compare the timing residual plots """
    main_args = locals()
    del main_args['divisor'] 
    switching, noswitch = get_data(**main_args) 

    # Print values
    File_Functions.PrintParameterDictionary(switching)

    fig, (ax1, ax2, ax3) = plt.subplots(nrows=3)
    ax1 = Plot.nu_dot(switching, divisor=divisor, ax=ax1,
                              label="$\upsilon={:1.1f}$".format(upsilon))
    Plot.nu_dot(noswitch, divisor=divisor,  ax=ax1,
                              label="$\upsilon=0.0$")

    ax2 = Plot.timing_residual(switching, ax=ax2,
                              label="$\upsilon={:1.1f}$".format(upsilon))
    ax2 = Plot.timing_residual(noswitch, ax=ax2, 
                              label="$\upsilon=0.0$")    

    Phi0=180; Theta0=50; sigmaPhi=0.3; sigmaTheta=0.3; eta=0.01
    ax3 = Plot.PulseWidth(switching, Phi0, Theta0, sigmaPhi, sigmaTheta, ax=ax3, 
                          eta=eta, label="$\upsilon={:1.1f}$".format(upsilon),
                          ls="-",)
    ax3 = Plot.PulseWidth(noswitch, Phi0, Theta0, sigmaPhi, sigmaTheta, ax=ax3, 
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

ExamplePlot()
