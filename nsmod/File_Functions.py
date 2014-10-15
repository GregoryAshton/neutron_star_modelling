#!/usr/bin/pltthon

import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
import h5py
import os

def RemoveFileSuffix(file_name):
    suffix = ".hdf5"
    if file_name.endswith(suffix):
        return file_name[:-5]

def Save_Figure(file_path, type_of_plot, dest_dir="img", format_type=".pdf"):
    """ Saves the currently open figure with an appropriate name"""
    if not os.path.isdir(dest_dir):
        os.mkdir(dest_dir)

    file_name = file_path.split("/")[-1]
    plot_file_name = "./{0}/{1}_{2}{3}".format(
                                        dest_dir, type_of_plot, 
                                        RemoveFileSuffix(file_name),
                                        format_type)
 
    plt.tight_layout()                  
    plt.savefig(plot_file_name)
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
        f = RemoveFileSuffix(user_input)
        f = f.split("/")[-1]
        
        # Import the rest of the parameters
        f = f.split("_")
        for i in range(0, len(f), 2):
            p_d[f[i]] = float(f[i + 1])

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
        print epsI, epsI3
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
    p_d["tauP"] = 2 * pi * pow(omega0 * epsI, -1)
    if epsA != 0.0:
        p_d["tauA"] = str(2 * pi * pow(omega0 * epsA, -1))
        p_d["tauS"] = str(pow(2 * pi, 2) * pow(omega0 ** 2.0 * epsA, -1)
                                                    * 3 * c / (2 * R))
    Bs = (2 * np.sqrt(epsA * I0 * R * pow(c, 2)) / pow(R, 3))
    p_d["Bs"] = str(Bs)
    chi0 = np.radians(p_d['chi0'])
    a0 = np.radians(p_d['a0'])
    Sx = np.sin(chi0)
    Cx = np.cos(chi0)
    varphi = 0.0
    alpha = np.arccos(Sx * np.sin(a0) * np.cos(varphi) + Cx * np.cos(a0))
    omega_dot0 = -2/3. * omega0**3 * R/c  * np.sin(alpha)**2 * epsA
    p_d['omega_dot0'] = omega_dot0

    p_d['delta_omega_dot0_FP'] = epsI**2 * a0 * np.cos(chi0) * omega0**2 / (
                                   np.sin(chi0))
    # Need to import the beta function
    #from Physics_Functions import Beta_Function
    #p_d["beta30"] = str(Beta_Function(epsI, epsA, 30 * pi / 180) * 180 / pi)
    #p_d["beta75"] = str(Beta_Function(epsI, epsA, 75 * pi / 180) * 180 / pi)

    return p_d

def PrintParameterDictionary(file_name):
    pd = Parameter_Dictionary(file_name)
    print "File: {}".format(file_name)
    for key, val in pd.iteritems():
        try:
            formatted_val = "{:1.4e}".format(float(val))
            print key, ":", formatted_val 
        except ValueError:
            print key, ":", val

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


def Euler_Angles_Import(file_name, time_offset=None):
    """ Returns tuple of "(time, w1, w2, w3, theta, phi, psi) 
    
    Parameters
    ----------
    file_name : str
        String containing file to import from
    time_offset : float [0, 1]
        Removes the first fraction of data
    
    """    

    with Read_File(file_name) as f:
        time = np.array(f['time'].value)
        w1 = np.array(f['w1'].value)
        w2 = np.array(f['w2'].value)
        w3 = np.array(f['w3'].value)
        theta = np.array(f['theta'].value)
        phi = np.array(f['phi'].value)
        psi = np.array(f['psi'].value)
    
    if time_offset:
        idx = np.argmin(np.abs(time-time_offset))
        time = time[idx:]
        w1 = w1[idx:]
        w2 = w2[idx:]
        w3 = w3[idx:]
        theta =theta[idx:]
        phi = phi[idx:]
        psi = psi[idx:]
    
    return (time, w1, w2, w3, theta, phi, psi)


def Clean_Data(directory, ask_user=True):
    """ Remove all .hdf5 files in directory """

    if ask_user:
        print ("WARNING: Permanent removal of all '.hdf5' files in the "
               "directory {}, do you wish to proceed? \n y/n".format(
                                                            directory))
        answer = raw_input()
    else:
        answer="yes"

    if answer in['yes', 'y', 'alright']:
        file_type = 'hdf5'
        data_files = [dfile for dfile in os.listdir(directory) 
                                            if file_type in dfile]
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
    if key in ["epsI1", "epsI3", "epsA", "omega0", "T", "SwitchTime"]:
        formatted_val = "{:.2e}".format(val)
    elif key in ["chi0", "a0"]:
        formatted_val = "{:2.5f}".format(val)
    elif key in ["AnomTorque"]:
        formatted_val = "{:.0f}".format(val)
    elif key in ["upsilon", "eta"]:
        formatted_val = "{:.3f}".format(val)
    elif key in ["n"]:
        if val:
            formatted_val = "{:.0f}".format(val)
        else:
            formatted_val = None  
    elif key in ["error", "cleanup"]:
        formatted_val = None
    else:
        raise ValueError(
          "The pair {}:{} do not match known values".format(key, val))

    return formatted_val

def FileNamer(data_dir="data", **kwargs):
    """ Returns a file name describing the key word arguments """

    # Check if data folder exists
    data_dir = "./"+data_dir+"/"
    if not os.path.isdir(data_dir):
        os.mkdir(data_dir)

    # Create the file name
    file_name_list = []
    for key, val in kwargs.iteritems():
        formatted_val = FormatValue(key, val)
        if formatted_val:
            file_name_list.append("{}_{}_".format(key, formatted_val))

    file_name_list[-1] = file_name_list[-1].rstrip("_")
    file_name_list.append(".hdf5")
    file_name = "".join(file_name_list)

    file_path = data_dir + file_name
    run_sim = True
    if os.path.isfile(file_path):
        if kwargs['cleanup']:
            os.remove(file_path)
        else:
            print("WARNING: File already exists and cleanup=False")
            run_sim = False

    return file_path, run_sim
