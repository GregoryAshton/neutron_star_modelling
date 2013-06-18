#!/usr/bin/python

import numpy as np
import pylab as py
from math import pi


def Cartesian_2_Spherical(x, y, z, Angle_Type="Degrees", fix_phi=False):
    """ Transform x,y,z to radial,polar and azimuthal vectors"""

    N = len(x)

    radial = [(x[i] * x[i] + y[i] * y[i] + z[i] * z[i]) ** 0.5
                for i in range(N)]
    polar = [py.arccos(z[i] / radial[i]) for i in range(N)]
    azimuth = [py.arctan(y[i] / x[i]) for i in range(N)]

    if "Degrees" in Angle_Type:

        polar_degrees = [p * 180 / pi for p in polar]
        azimuth_degrees = [a * 180 / pi for a in azimuth]

        if fix_phi:
            azimuth_degrees = Fix_Phi(azimuth_degrees, Angle_Type="Degrees")

        return (radial, polar_degrees, azimuth_degrees)

    elif Angle_Type in ["Radians", "Radian", "Rads"]:

        if fix_phi:
            azimuth = Fix_Phi(azimuth, Angle_Type="Radians")

        return (radial, polar, azimuth)


def Cartesian_2_EBF(x, y, z, beta):
    """Functions to rotate x,z system by angle beta about the y axis"""

    N = len(x)
    Cb = py.cos(beta)
    Sb = py.sin(beta)
    x_prime = [x[i] * Cb - z[i] * Sb for i in xrange(N)]
    z_prime = [z[i] * Cb + x[i] * Sb for i in xrange(N)]
    y_prime = y
    return (x_prime, y_prime, z_prime)


def Rotational_Kinetic_Energy(Ix, Iy, Iz, omega_x, omega_y, omega_z):
    en_2 = (Ix * pow(omega_x, 2) + Iy * pow(omega_y, 2) + Iz * pow(omega_z, 2))
    return 0.5 * en_2


def Fix_Phi(phi, epsilon=170.0, Angle_Type="Degrees"):
    """

    Takes a list of phi values and looks for jumps greater than epsilon,
    assuming these reflect a full rotation (2pi-0) it adds a correction
    factor to the subsequent data

    """

    if Angle_Type == "Degrees":
        phi_fix = []
        fix = 0.0
        phi_fix.append(phi[0])
        for i in range(1, len(phi)):
            if abs(phi[i] - phi[i - 1]) > epsilon:
                fix += -1 * py.sign(phi[i] - phi[i - 1]) * 180.0
            phi_fix.append(phi[i] + fix)
        phi = phi_fix
        return phi

    elif Angle_Type in ["Radians", "Radian", "Rads"]:
        if epsilon == 170.0:
            epsilon = 1.0

        phi_fix = []
        fix = 0.0
        phi_fix.append(phi[0])
        for i in range(1, len(phi)):
            if abs(phi[i] - phi[i - 1]) > epsilon:
                fix += -1 * py.sign(phi[i] - phi[i - 1]) * pi
            phi_fix.append(phi[i] + fix)
        phi = phi_fix
        return phi


def Beta_Function(epsI, epsA, chi):
    """ Returns beta the rotation of the effective MOI tensor """
    if chi > 2 * pi:
        print "Assuming chi has been given in degrees rather than radians"
        chi = chi * pi / 180

    a = epsA * epsA + epsI * epsI - 2 * epsA * epsI * py.cos(2 * chi)
    beta = (py.arctan((epsI - epsA * py.cos(2 * chi) - py.sqrt(a)) /
                        (2 * epsA * py.sin(chi) * py.cos(chi))))
    return beta


def Inertial_Frame(w_1, w_2, w_3):
    """

    Transformation from rotating coordinate system x,y,z to inertial frame

    """

