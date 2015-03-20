""" A collection of signal models for the spin-down """
from numpy import sin, cos, pi

def SignalModelWithGeometric(params, t, geometric=1):
    omega0, epsI, a0, chi, epsA = params

    theta = a0

    c = 3e10
    R = 1e6
    k = 2/3.0 * R/c * epsA

    psi = -epsI*omega0*t + pi/2 + 0.5*k*epsI*sin(chi)**2*omega0**3*t**2

    C = 1 - (cos(theta)*cos(chi))**2 - 0.5*(sin(theta)*sin(chi))**2

    T1 = -k * omega0**3 * C
    T2 =  3*k**2*omega0**5 * C**2 * t

    T3 = 0.5 * k * omega0**3 *(sin(2*theta)*sin(2*chi)*sin(psi) -
                               (sin(theta)*sin(chi))**2*cos(2*psi))

    theta = a0
    psidot = -epsI*omega0
    #GEOMETRIC = psidot**2 * ((2*sin(chi)**3*sin(psi)*sin(theta)*cos(theta) -
    #                          sin(chi)**2*sin(psi)**2*sin(theta)**2*cos(chi) -
    #                          2*sin(chi)**2*sin(theta)**2*cos(chi) +
    #                          sin(chi)**2*cos(chi) - sin(theta)**2*cos(chi)**3
    #                         )*sin(chi)*sin(theta)*cos(psi)/(
    #                        (sin(chi)*sin(psi)*cos(theta) - sin(theta)*cos(chi))**2 +
    #                         sin(chi)**2*cos(psi)**2)**2
    #                        )
    GEOMETRIC = psidot**2 * (theta * cos(chi)/sin(chi) * cos(psi) + 
                             .5*theta**2*(3+cos(2*chi))/(sin(chi))* sin(2*psi)
                             )
    return  (T1 + T2 + T3 + geometric * GEOMETRIC)/(2*pi)
