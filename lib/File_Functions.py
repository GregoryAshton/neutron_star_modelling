#!/usr/bin/python

import numpy as np
from numpy import pi
import pylab as py
import h5py
import os


def Save_Figure(file_name, type_of_plot, format_type=".pdf"):
    """ Saves the currently open figure with an appropriate name"""
    plot_file_name = "{0}_{1}{2}".format(
                        type_of_plot, file_name.rstrip(".hdf5"), format_type
                            )
    py.savefig(plot_file_name)
    print "Saving figure as %s" % plot_file_name


def Parameter_Dictionary(user_input):
    # Could we improve this to give units?
    """

    Function to produce a dictionary with all the values, this reduces the
    amount of code required on the other end. Input can either be a string
    such as the default filenames used in nsmod, or a partial dictionary
    having at least epsI,epsA and omega0 to produce the other values

    """

    if type(user_input) == str:
        # Initiate dictionary
        p_d = {}

        # Remove the file descriptor and the path directory
        f = user_input.rstrip(".hdf5")
        f = f.split("/")[-1]

        # Check if the anomalous torque was used or not
        if "no_anom" in f:
            f = f.lstrip("no_anom_")
            p_d["anom_torque"] = False
        else:
            p_d["anom_torque"] = True

        # Import the rest of the parameters
        f = f.split("_")
        for i in range(0, len(f), 2):
            p_d[f[i]] = f[i + 1]

    elif type(user_input) is dict:
        p_d = user_input

    # Standard values for c , R in cgs
    c = 3e10
    R = 1e6
    I0 = 1e45
    # Compute a couple of often used variabes

    # Test Biaxiality, note that the tauP does not make sense in the triaxial
    # case, need to update this at some point
    try:
        epsI3 = float(p_d["epsI3"])
        epsI1 = float(p_d["epsI1"])
        epsI = max(abs(epsI3), abs(epsI1))
    except KeyError:
        try:
            epsI = float(p_d["epsI"])
            p_d["epsI1"] = "0.0"
            p_d["epsI3"] = str(epsI)
        except KeyError:
            print " ERROR: No epsI specified"
            return

    omega0 = float(p_d["omega0"])
    epsA = float(p_d["epsA"])
    p_d["tauP"] = str(2 * pi * pow(omega0 * epsI, -1))
    if epsA == 0.0:
        pass
        #p_d["tauA"] = float("inf")
    else:
        p_d["tauA"] = str(2 * pi * pow(omega0 * epsA, -1))
        p_d["tauS"] = str(pow(2 * pi, 2) * pow(omega0 ** 2.0 * epsA, -1)
                                                        * 3 * c / (2 * R))
        p_d["Bs"] = str(2 * np.sqrt(epsA * I0 * R * pow(c, 2)) / pow(R, 3))

        # Need to import the beta function
        from Physics_Functions import Beta_Function
        p_d["beta30"] = str(Beta_Function(epsI, epsA, 30 * pi / 180) * 180 / pi)
        p_d["beta75"] = str(Beta_Function(epsI, epsA, 75 * pi / 180) * 180 / pi)

    return p_d


def One_Component_Import(file_name, nmax=None):
    """ Imports time and w1,w2,w3 from file_name """
    # max_int and d_int are obsolete for now
    f = Read_File(file_name)
    time = f['time'].value
    w1 = f['w1'].value
    w2 = f['w2'].value
    w3 = f['w3'].value
    f.close()
    if nmax:
        return (time[0:nmax], w1[0:nmax], w2[0:nmax], w3[0:nmax])
    else:
        return (time, w1, w2, w3)


def Two_Component_Import(file_name):
    """ Imports time and w1,w2,w3,o1,o2,o3 from file_name """
    # max_int and d_int are obsolete for now
    f = Read_File(file_name)
    time = f['time'].value
    w1 = f['w1'].value
    w2 = f['w2'].value
    w3 = f['w3'].value
    o1 = f['o1'].value
    o2 = f['o2'].value
    o3 = f['o3'].value
    f.close()

    return (time, w1, w2, w3, o1, o2, o3)


def Read_File(file_name):
    """

    This opens a .hdf5 file and returns an effective dictionary.
    To access values >>>f['key'].value

    """

    # Check the file name is well formed
    file_type = file_name.split(".")[-1]

    if file_type == "hdf5":
        pass
    else:
        print "File type {} is ill-formed.".format(file_type)
        return

    f = h5py.File(file_name, "r")

    return f


def vprint(verbose, *args):
    """ Function to verbose print """
    if verbose:
        # Print each argument separately so caller doesn't need to
        # stuff everything to be printed into a single string
        for arg in args:
            print arg,
        print
    else:
        vprint = lambda *a: None      # do-nothing function
