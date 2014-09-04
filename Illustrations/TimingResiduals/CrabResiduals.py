import matplotlib.pylab as plt
import numpy as np
import matplotlib as mpl
mpl.rc_file("/home/greg/Neutron_star_modelling/matplotlibrc")
mpl.rcParams['axes.formatter.limits'] = (-3, 5)

from scipy.integrate import cumtrapz
from scipy.optimize import curve_fit

CrabEphemerisFile = "/home/greg/timing-noise/CrabEphemeris/crab2.txt"

def get_residuals(file_name, start_date, end_date, order=3):
    """ Plot the timing residuals for the Crab timing ephemeris from start_date to end_date in MJD"""
    
    data = np.genfromtxt(file_name, skip_header=6, skip_footer=2, 
                         usecols=(3, 4, 5, 7, 9), autostrip=True)
    
    MJD = data[:, 0]
    
    # Locate index in MJD closest to start_date and end_date
    low_idx = min(range(len(MJD)), key=lambda i: abs(MJD[i]-start_date))
    high_idx = min(range(len(MJD)), key=lambda i: abs(MJD[i]-end_date))

    if low_idx >=high_idx:
        print """The start and end dates provided are either outside the range of the filename
                 or, you put them in the wrong order. 
                 
                 Available dates are {} to {} with {} points inbetween
              """.format(int(MJD[0]), int(MJD[1]), len(MJD))
    
    MJD = data[low_idx:high_idx, 0] # Modifed Julian Date of the Ephemeris 
    TOA_after_midnight = data[low_idx:high_idx, 2] # [s]
    nu = data[low_idx:high_idx, 3] # [s^-1]
    nu_dot = data[low_idx:high_idx, 4]*1e-15 # [s^-2]    
    TOA = MJD*24*3600 + TOA_after_midnight # [s]
   
    try:
        phi = cumtrapz(y=nu, x=TOA, intial=0.0)

    except TypeError:
        phi = cumtrapz(y=nu, x=TOA)
        phi = np.concatenate(([0.0], phi))

    coefs = np.polyfit(TOA, phi, order)
    phi_fit = np.poly1d(coefs)(TOA)
    phase_residual = phi - phi_fit 
    
    return MJD, TOA, phase_residual
    
    
def plot_residuals(file_name, start_date, end_date, order=3, label=None):
    """ Plotting function """
    
    plt.axhline(0, ls="-", color="k", lw=0.5)

    MJD, TOA, phase_residual = get_residuals(file_name, start_date, end_date, order)
    
    plt.plot(MJD, phase_residual, label=label)
        
    plt.xlabel("Modified Julian Date")
    plt.ylabel("Phase residual (cycles)")
    plt.savefig("Crab_residual_phase.pdf")
    plt.show()


plot_residuals(CrabEphemerisFile, start_date=54540, end_date=55519)
