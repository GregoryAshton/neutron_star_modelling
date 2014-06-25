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
    time = np.array(f['time'].value)
    w1 = np.array(f['w1'].value)
    w2 = np.array(f['w2'].value)
    w3 = np.array(f['w3'].value)
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


def Euler_Angles_Import(file_name):
    """ Returns tuple of "(time, w1, w2, w3, theta, phi, psi) """    
    f = Read_File(file_name)
    time = np.array(f['time'].value)
    w1 = np.array(f['w1'].value)
    w2 = np.array(f['w2'].value)
    w3 = np.array(f['w3'].value)
    theta = np.array(f['theta'].value)
    phi = np.array(f['phi'].value)
    psi = np.array(f['psi'].value)
    f.close()
    return (time, w1, w2, w3, theta, phi, psi)


def Clean_Data(directory):
    """ Remove all .hdf5 files in directory """

    print ("WARNING: Permanent removal of all '.hdf5' files in the directory {}, do you wish to proceed? \n y/n".format(directory))
    answer = raw_input()
    if answer in['yes', 'y', 'alright']:
        file_type = 'hdf5'
        data_files = [dfile for dfile in os.listdir(directory) if file_type in dfile]
        for dfile in data_files:
            os.remove(dfile)

    else:
        print "Okay no files will be deleted"
        return

def FormatValue(key, val):
    """ Formats value to give sensible file name 
    
    Note: None is returned if the key is not used in the filename, e.g 
          the error.
    """
    if key in ["epsI1", "epsI3", "epsA", "omega0", "T"]:
        formatted_val = "{:.2e}".format(val)
    elif key in ["chi0", "a0"]:
        formatted_val = "{:2.2f}".format(val)
    elif key in ["n"]:
        formatted_val = "{:.0f}".format(val)
    elif key in ["AnomTorque"]:
        formatted_val = "{:.0f}".format(val)
    elif key in ["upsilon"]:
        formatted_val = "{:.3f}".format(val)
    elif key in ["error"]:
        formatted_val = None
    else:
        raise ValueError(
          "The pair {}:{} do not match known values".format(key, val))

    return formatted_val

def FileNamer(**kwargs):
    """ Returns a file name describing the key word arguments """
   
    file_name_list = []
    for key, val in kwargs.iteritems():
        formatted_val = FormatValue(key, val)
        if formatted_val:
            file_name_list.append("{}_{}_".format(key, formatted_val))

    file_name_list[-1] = file_name_list[-1].rstrip("_")
    file_name_list.append(".hdf5")
    file_name = "".join(file_name_list)

    return file_name
