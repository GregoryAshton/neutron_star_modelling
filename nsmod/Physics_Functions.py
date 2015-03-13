#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from math import pi
from numpy import cos, sin, tan
import scipy as sc

from scipy.integrate import cumtrapz



def Cartesian_2_Spherical(x, y, z, Angle_Type="Degrees", fix_varphi=False):
    """ Transform x,y,z to radial,polar and azimuthal vectors"""

    N = len(x)

    radial = [(x[i] * x[i] + y[i] * y[i] + z[i] * z[i]) ** 0.5
                for i in range(N)]
    polar = [np.arccos(z[i] / radial[i]) for i in range(N)]
    azimuth = [np.arctan(y[i] / x[i]) for i in range(N)]

    if "Degrees" in Angle_Type:

        polar_degrees = [p * 180 / pi for p in polar]
        azimuth_degrees = [a * 180 / pi for a in azimuth]

        if fix_varphi:
            azimuth_degrees = Fix_Varphi(azimuth_degrees, Angle_Type="Degrees")

        return (radial, polar_degrees, azimuth_degrees)

    elif Angle_Type in ["Radians", "Radian", "Rads"]:

        if fix_varphi:
            azimuth = Fix_Varphi(azimuth, Angle_Type="Radians")

        return (radial, polar, azimuth)


def Cartesian_2_EBF(x, y, z, beta):
    """Functions to rotate x,z system by angle beta about the y axis"""

    N = len(x)
    Cb = np.cos(beta)
    Sb = np.sin(beta)
    x_prime = [x[i] * Cb - z[i] * Sb for i in xrange(N)]
    z_prime = [z[i] * Cb + x[i] * Sb for i in xrange(N)]
    y_prime = y
    return (x_prime, y_prime, z_prime)

def epsA(Bs, R=1e6, I0=1e45, c=1e10):
    ''' Return epsA from standard equations '''
    return np.power(Bs, 2) * np.power(R, 5) / (4 * I0 * np.power(c, 2))

def Rotational_Kinetic_Energy(Ix, Iy, Iz, omega_x, omega_y, omega_z):
    en_2 = (Ix * pow(omega_x, 2) + Iy * pow(omega_y, 2) + Iz * pow(omega_z, 2))
    return 0.5 * en_2


def Fix_Varphi(varphi, epsilon=170.0, Angle_Type="Degrees"):
    """

    Takes a list of varrphi values and looks for jumps greater than epsilon,
    assuming these reflect a full rotation (2pi-0) it adds a correction
    factor to the subsequent data

    """

    if Angle_Type == "Degrees":
        varphi_fix = []
        fix = 0.0
        varphi_fix.append(varphi[0])
        for i in range(1, len(varphi)):
            if abs(varphi[i] - varphi[i - 1]) > epsilon:
                fix += -1 * np.sign(varphi[i] - varphi[i - 1]) * 180.0
            varphi_fix.append(varphi[i] + fix)
        varphi = varphi_fix
        return varphi

    elif Angle_Type in ["Radians", "Radian", "Rads"]:
        if epsilon == 170.0:
            epsilon = 1.0

        varphi_fix = []
        fix = 0.0
        varphi_fix.append(varphi[0])
        for i in range(1, len(varphi)):
            if abs(varphi[i] - varphi[i - 1]) > epsilon:
                fix += -1 * np.sign(varphi[i] - varphi[i - 1]) * pi
            varphi_fix.append(varphi[i] + fix)
        varphi = varphi_fix
        return varphi


def Beta_Function(epsI, epsA, chi):
    """ Returns beta the rotation of the effective MOI tensor """
    print "BETA FUNCTION IS USED"
    if chi > 2 * pi:
        print "Assuming chi has been given in degrees rather than radians"
        chi = np.radians(chi)

    a = epsA * epsA + epsI * epsI - 2 * epsA * epsI * np.cos(2 * chi)
    beta = (np.arctan2((epsI - epsA * np.cos(2 * chi) - np.sqrt(a)),
                        (2 * epsA * np.sin(chi) * np.cos(chi))))
    return beta


def Inertial_Frame(time, omega, epsI3):
    """ Return the Euler angles at the timesteps of omega """

    def equations(omega, theta, phi, psi):
        theta_dot = omega[0] * np.cos(psi) - omega[1] * np.sin(psi)
        phi_dot = ((omega[0] * np.sin(psi) + omega[1] * np.cos(psi))
                                                      / np.sin(theta))
        psi_dot = omega[2] - phi_dot * np.cos(theta)
        return theta_dot, phi_dot, psi_dot

    # Initial conditions
    #theta = [np.arccos((pow(epsI1, 0.5) * omega[0, 0] * (1 + epsI1)
                      #+ pow(epsI3, 0.5) * omega[2, 0] * (1 + epsI3))
                        #/ (pow(epsI1 + epsI3, 0.5)
                        #* np.norm([omega[0, 0] * (1 + epsI1),
                                   #omega[1, 0],
                                   #omega[2, 0] * (1 + epsI3)])))]
    theta = [np.arccos(omega[2, 0] * (1 + epsI3) /
            (np.norm([omega[0, 0], omega[1, 0], omega[2, 0] * (1. + epsI3)])))]
    phi = [0.0]
    psi = [np.pi / 2]  # Unfinished

    for i in xrange(len(time)-1):
        dt = time[i + 1] - time[i]
        theta_dot, phi_dot, psi_dot = equations(omega[:, i], theta[i],
                                                phi[i], psi[i])
        theta.append(theta[i] + dt * theta_dot)
        phi.append(phi[i] + dt * phi_dot)
        psi.append(psi[i] + dt * psi_dot)

    return theta, phi, psi

# Equations for the Euler angles 

def equations(omega, theta, phi, psi):
    theta_dot = omega[0] * cos(psi) - omega[1] * sin(psi)
    phi_dot = ((omega[0] * sin(psi) + omega[1] * cos(psi))
                                                  / sin(theta))
    psi_dot = omega[2] - phi_dot * cos(theta)
    return theta_dot, phi_dot, psi_dot

def Phi_dot(omega, theta, phi, psi, chi, theta_dot=None, phi_dot=None, 
            psi_dot=None):
    """ 

    Explicitly pass omega=None to use theta_dot, phi_dot
    """
    if omega!=None:
        theta_dot, phi_dot, psi_dot = equations(omega, theta, phi, psi)
    return phi_dot +  sin(chi) * ( 
            (psi_dot * (cos(theta) * sin(chi) - sin(psi) * sin(theta) * cos(chi)) +
             theta_dot * cos(psi) * (sin(psi) * sin(chi) * sin(theta) - cos(chi) * cos(theta))) /
	    (pow(sin(theta) * cos(chi) - cos(theta) * sin(psi) * sin(chi), 2) + pow(cos(psi) * sin(chi), 2)))
    
def Phi(theta, phi, psi, chi, fix=True):
    """ See equation (42) of Jones 2001 """
    A = phi - .5 * np.pi
    B = np.arctan2(
                (cos(psi) * tan(chi)),
                (sin(theta) - cos(theta) * sin(psi) * tan(chi)))
    if fix:
        B = np.unwrap(B)
    return A + B 

def Theta(theta, psi, chi):  
    """ See equation (52) of Jones 2001 """
    return np.arccos(sin(theta) * sin(psi) * sin(chi) + cos(theta) * cos(chi))  

def PhaseResidual(time, w1, w2, w3, theta, phi, psi, chi, order=3,
                    full=False):
    """ 

    Calculate the phase residuals in the inertial frame using Phi_dot the 
    instantaneous electromagnetic frequency. To understand the process it is 
    best to follow the source code.

    Parameters
    ----------
    file_name : str
               hdf5 file to use as source data.
    order : int
            Order of polynomial to fit to the exact phase, options are 2 or 3
    full : bool
           If true return the coefficients of the fit as w;;

    Returns
    -------
    r : tuple
        if full = True : time, T_res, coeffs
        else time, T_res

    """

    if order not in [2, 3]:
        print "ERROR: order must be either 2 or 3 other values are not used"
    #    return

    # Calculate Phi_dot the instantaneous electromagnetic frequency
    #Phi_dot_list = Phi_dot(theta, phi, psi, chi, omega=np.array([w1, w2, w3]))

    # Numerically intergrate Phi_dot to get a phase (initial conditon is Phi=0
    #Phi_list = cumtrapz(y=Phi_dot_list, x=time, initial=0)

    Phi_list = Phi(theta, phi, psi, chi, fix=True)

    # Fit polynomial to Phi or order order
    coefs, V = np.polyfit(time, Phi_list, order, cov=True)

    # poly1d returns the polynomial we then evaluate this at time giving the fitted phi
    Phi_fit = np.poly1d(coefs)(time)

    # Subtract the two to get a residual
    Phi_res = Phi_list - Phi_fit

    if full:
        return Phi_res, coefs
    else:
        return Phi_res

def nu_dot_Lyne(time, w1, w2, w3, theta, phi, psi, chi, tauP, divisor=7):
    """
    
    Calculate the spin down rate using the method precesribed by Lyne 2010 
    
    Parameters:
    -----------
    data: The usual (expand)
    divisor: The fraction of the precession period to use as a segmentation 
             time for calculating the nu_dot values. 

    Returns:
    out: Array of the values of the time and nu_dot

    """



    # Calculate Phi_dot the instantaneous electromagnetic frequency
    #Phi_dot_list = Phi_dot(theta, phi, psi, chi, omega=np.array([w1, w2, w3]))

    ## Numerically intergrate Phi_dot to get a phase (initial conditon is Phi=0
    #Phi_list = cumtrapz(y=Phi_dot_list, x=time, initial=0)

    Phi_list = Phi(theta, phi, psi, chi, fix=True)

    # Convert T into an index range of time, note all time intervals should be uniform
    if tauP < (time[-1] - time[0]):
        T = tauP / divisor
    else:
        T = (time[-1] - time[0]) / divisor
    dt = time[1] - time[0]
    T_index_range = int(T / dt)

    dT = int(0.25 * T_index_range)

    def get_nudot(time, Phi, tref):

        timeprime = time - tref
        poly, pcov = np.polyfit(timeprime, Phi, 2, cov=True)
        return poly[0] / np.pi


    nu_dot_list = []
    time_list = []
    i=0
    tref = time[0]

    while i < len(time)-T_index_range:
        time_mid = 0.5*(time[i] + time[i+T_index_range])
        nu_dot = get_nudot(time[i:i + T_index_range], 
                           Phi_list[i:i + T_index_range],
                           time_mid)
        nu_dot_list.append(nu_dot)
        time_list.append(time_mid) 
    
        i += dT

    return np.array(time_list), np.array(nu_dot_list)

def Ndiff(x, y):
    """ Numeric derivative assuming x is linearly spaced """
    dx = x[1] - x[0]
    return np.gradient(y) / dx

def nu_dot_numeric(time, theta, psi, phi, chi):
    phiddot = Ndiff(time, Ndiff(time, phi))
    psidot = Ndiff(time, Ndiff(time, psi))
    psiddot = Ndiff(time, Ndiff(time, psidot))
    f = sin(chi) * (cos(theta)*sin(chi) - sin(psi)*sin(theta)*cos(chi))/(
                   (sin(theta)*cos(chi) - cos(theta)*sin(psi)*sin(chi))**2 +
                   (cos(psi)*sin(chi))**2
                   )
    dfdpsi = ((2*sin(chi)**3*sin(psi)*sin(theta)*cos(theta) - 
                              sin(chi)**2*sin(psi)**2*sin(theta)**2*cos(chi) - 
                              2*sin(chi)**2*sin(theta)**2*cos(chi) + 
                              sin(chi)**2*cos(chi) - sin(theta)**2*cos(chi)**3
                             )*sin(chi)*sin(theta)*cos(psi)/(
                            (sin(chi)*sin(psi)*cos(theta) - sin(theta)*cos(chi))**2 + 
                             sin(chi)**2*cos(psi)**2)**2
                            )

    return time, phiddot + psiddot*f + psidot**2 * dfdpsi

    
def Amplitude(Phi, Theta, PhiO, ThetaO, sigmaB=0.01, A0=1):
    DeltaSigma = np.arccos(np.sin(Theta)*np.sin(ThetaO) + 
                           np.cos(Theta)*np.cos(ThetaO)*np.cos(np.abs(Phi - PhiO)))
    # Is this real? 
    return A0 *  np.exp(-DeltaSigma**2 / (2*sigmaB**2))

def Wp(Phi_dot, Theta, ThetaO, sigmaB, p=50):
    """ Analytic calculation of the pulse width 
    
    Parameters
    ----------
    Phi_dot, Theta : array_like
        The physical parameters as defined in Jones 2001
    ThetaO : float
        The value, in radians of the observers Theta position
    sigmaB : float
        Measure of the angular beam width
    p : float
        The percentage amount of beam width

    Returns
    -------
    Beam width at p
    """
    A = np.cos(np.sqrt(2 * sigmaB**2 * np.log(100./p))) 
    B = (A-np.sin(Theta) * np.sin(ThetaO)) / (np.cos(Theta) * np.cos(ThetaO))
    C = np.pi * Phi_dot
    return (1 - np.arccos(B) ) / C 
