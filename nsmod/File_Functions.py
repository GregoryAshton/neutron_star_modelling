#!/usr/bin/pltthon

import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
import h5py
import os
import Physics_Functions, Useful_Tools

def RemoveFileSuffix(file_name):
    suffix = ".hdf5"
    if file_name.endswith(suffix):
        return file_name[:-5]

def Save_Figure(file_path, type_of_plot, dest_dir="img", format_type=".png",
                tight=True):
    """ Saves the currently open figure with an appropriate name"""
    if not os.path.isdir(dest_dir):
        os.mkdir(dest_dir)

    file_name = file_path.split("/")[-1]
    plot_file_name = "./{0}/{1}_{2}{3}".format(
                                        dest_dir, type_of_plot,
                                        RemoveFileSuffix(file_name),
                                        format_type)

    if tight:
        plt.tight_layout()
    plt.savefig(plot_file_name, dpi=500)
    print "Saving figure as %s" % plot_file_name


def Parameter_Dictionary(user_input):
    # Could we improve this to give units?
    """

    Function to produce a dictionary with all the values, this reduces the
    amount of code required on the other end. Input can either be a string
    such as the default filenames used in nsmod, or a partial dictionary
    having at least epsI1, epsI3,epsA and omega0 to produce the other values

    """

    if type(user_input) == str:
        # Initiate dictionary
        p_d = {}

        # Remove the file descriptor and the path directory
        f = RemoveFileSuffix(user_input)
        f = f.split("/")[-1]

        # Import the rest of the parameters
        f = f.split("_")
        p_d['source_script'] = f[0]
        for i in range(1, len(f), 2):
            p_d[f[i]] = float(f[i + 1])

    elif type(user_input) is dict:
        p_d = user_input

    # Standard values for c , R in cgs
    c = 3e10
    R = 1e6
    I0 = 1e45

    # Compute a couple of often used variabes
    epsI3 = float(p_d["epsI3"])
    epsI1 = float(p_d["epsI1"])
    epsA = float(p_d["epsA"])
    omega0 = float(p_d["omega0"])
    chi0 = np.radians(p_d['chi0'])
    a0 = np.radians(p_d['a0'])
    theta0 = a0 / (1+epsI3)
    p_d['theta0'] = theta0

    P = 2*np.pi / omega0
    p_d['P'] = P

    if epsI1 == 0 and epsI3 == 0:
        tauP = np.nan
    else:
        tauP = 2 * pi * pow(omega0 * abs(epsI3), -1) / np.cos(theta0)

    p_d["tauP"] = tauP

    delta_omega_dot0_FP = epsI3**2 * a0 * np.cos(chi0) * omega0**2 / (
                                   np.sin(chi0))
    p_d['delta_omega_dot0_FP'] = delta_omega_dot0_FP


    if epsA != 0.0:
        Sx = np.sin(chi0)
        Cx = np.cos(chi0)
        varphi = 0.0
        alpha = chi0 #np.arccos(Sx * np.sin(a0) * np.cos(varphi) + Cx * np.cos(a0)) # Dodgy

        tauA = 2 * pi * pow(omega0 * epsA, -1)
        p_d["tauA"] = str(tauA)

        tauS =  3 * c / (2 * R * epsA * omega0**2)
        p_d["tauS"] = str(tauS)


        Bs = (2 * np.sqrt(epsA * I0 * pow(c, 2) / pow(R, 5)))
        p_d["Bs"] = Bs

        omega_dot0 = -2 * R /(3. * c) * omega0**3 * np.sin(alpha)**2 * epsA
        p_d['omega_dot0'] = omega_dot0
        nu_dot0 = -1/(3 * np.pi) * (R/c) * omega0**3 * np.sin(alpha)**2 * epsA
        p_d['nu_dot0'] = nu_dot0

        tauE = abs(omega0/omega_dot0)
        p_d["tauE"] = tauE # Differs from tauS by alpha factor


    # Wobble angle calculation

    if epsA != 0 and epsI3 != 0:
        beta = Physics_Functions.Beta_Function(epsI3, epsA, chi0, warning=False)
        wobble_angle_spindown = (P / tauS) * (1 + 1.0/epsI3)
    else:
        beta = 0
        wobble_angle_spindown = 0

    p_d['beta'] = np.degrees(beta)
    p_d['wobble_angle_spindown'] = wobble_angle_spindown

    wobble_angle = a0 - beta + wobble_angle_spindown
    p_d['wobble_angle'] = wobble_angle

    DeltaPhi_49 = wobble_angle / np.tan(chi0)
    p_d['DeltaPhi_49'] = DeltaPhi_49

    if epsA != 0.0:
        p_d['DeltaPhi_75'] = tauP**2 * wobble_angle**2 / (4 * np.pi * tauE * P)

        if p_d.has_key('upsilon'):
            upsilon = p_d['upsilon']
            switching_period = p_d['T']
            DeltaPhi_TS = (1/16.0) * upsilon * omega_dot0 * switching_period**2
            p_d['DeltaPhi_TS'] = DeltaPhi_TS

        p_d['DeltaPhi_49_SpindownTorque'] = wobble_angle_spindown / np.tan(chi0)

        EMtorqueAmplificationfactor = (tauP / P) * (tauP / tauE)
        p_d['EMtorqueAmplificationfactor'] = EMtorqueAmplificationfactor

        #p_d['DeltaPhi_63'] = EMtorqueAmplificationfactor * DeltaPhi_49 / np.pi
        p_d['DeltaPhi_63'] = wobble_angle * EMtorqueAmplificationfactor / np.tan(chi0) / np.pi

        p_d['delta_omega_dot0_FP_EM'] = delta_omega_dot0_FP * EMtorqueAmplificationfactor / np.pi

        p_d['delta_omega_dot0_EM'] = pow(tauS * P, -1)

    # Need to import the beta function
    from Physics_Functions import Beta_Function
    p_d["beta30"] = Beta_Function(epsI3, epsA, 30 * pi / 180) * 180 / pi
    p_d["beta75"] = Beta_Function(epsI3, epsA, 75 * pi / 180) * 180 / pi

    return p_d

def PrintParameterDictionary(file_name):
    pd = Parameter_Dictionary(file_name)
    #print "File: {}".format(file_name)
    for key, val in sorted(pd.iteritems()):
        try:
            formatted_val = "{:1.10e}".format(float(val))
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

def PropertiesTable(file_name, table_name,
                    include=['omega0', 'B', 'chi', 'a0', 'Aem', 'wobble_angle']):
    """ Print a latex table of the NS properties used in a simulation

    Parameters
    ----------
    file_name : str
        Name of the simulation file
    table_name : str
        Name of the file to save data in, file will always be .tex
    include : list
        A list of attributes to save in the properties table. Possible values
        are: omega0, B, chi, a0, Aem, wobble_angle, tauP, tauS, tauA
    """

    table_name = table_name.split(".")[0]

    PD = Parameter_Dictionary(file_name)

    lines = [r"\begin{tabular}{ccl}",
             r"\multicolumn{3}{c}{Simulation parameters} \\",
             "\hline"]
    if 'omega0' in include:
        lines.append(r"$\omega_0$  &=& {} rad/s\\".format(PD['omega0']))
    if 'B' in include:
        try:
            Bs = Useful_Tools.Texify_Float(PD['Bs'], 3)
        except KeyError:
            Bs = 0
        lines.append(r"$B_0$  &=& ${}$ G \\".format(Bs))
    if 'chi' in include:
        lines.append(r"$\chi$  &=& {:2.2f}$^{{\circ}}$ \\".format(PD['chi0']))
    if 'a0' in include:
        lines.append(r"$a_0$ &=& {:2.2f}$^{{\circ}}$ \\".format(PD['a0']))
    if 'wobble_angle' in include:
        lines.append(r"$\tilde{{\theta}}$ &= & {:2.2f}$^{{\circ}}$ \\".format(
                       np.degrees(PD['wobble_angle'])))
    if 'Aem' in include:
        try:
            Aem = Useful_Tools.Texify_Float(PD['EMtorqueAmplificationfactor'])
        except KeyError:
            Aem = 0
        lines.append("$\mathcal{{A}}_{{\mathrm{{EM}}}}$ &= & ${}$".format(Aem))
    if 'tauP' in include:
        tauP = Useful_Tools.Texify_Float(PD['tauP'], 0)
        lines.append(r"$\tau_{{P}}$ &$\approx$& ${}$ \\".format(tauP))
    if 'tauA' in include:
        try:
            tauA = Useful_Tools.Texify_Float(PD['tauA'], 0)
            lines.append(r"$\tau_{{A}}$ &$\approx$& ${}$ \\".format(tauA))
        except KeyError:
            print "WARNING: tauA requested but not in PD. Ignoring."
    if 'tauS' in include:
        try:
            tauS = Useful_Tools.Texify_Float(PD['tauS'], 0)
            lines.append(r"$\tau_{{S}}$ &$\approx$& ${}$ \\".format(tauS))
        except KeyError:
            print "WSRNING: tauS requested but not in PD. Ignoring."


    lines.append("\end{tabular}")

    table = "\n".join(lines)
    with open(table_name+".tex", "w+") as f:
        f.write(table)


def Read_File(file_name):
    """

    This opens a .hdf5 file and returns an effective dictionary.
    To access values >>>f['key'].value

    """
    print file_name

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


def Euler_Angles_Import(file_name, time_offset=None, nmax=None):
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

    return (time[:nmax], w1[:nmax], w2[:nmax], w3[:nmax],
            theta[:nmax], phi[:nmax], psi[:nmax])


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
        formatted_val = "{:1.10e}".format(val)
    elif key in ["AnomTorque", "SpindownTorqueSwitching", "AnomTorqueSwitching"]:
        formatted_val = "{:.0f}".format(val)
    elif key in ["upsilon", "eta"]:
        formatted_val = "{:1.2e}".format(val)
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

def FileNamer(data_dir="data", source_script="Unknown", **kwargs):
    """ Returns a file name describing the key word arguments """

    # Check source-scipt has sensible name
    if "_" in source_script:
        raise ValueError("source_script cannot contain '_' as this will "
                         "interfer with the file naming system. Suggest using - "
                         "instead")
    # Check if data folder exists
    data_dir = "./"+data_dir+"/"
    if not os.path.isdir(data_dir):
        os.mkdir(data_dir)

    # Create the file name
    file_name_list = [source_script+"_"]
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
