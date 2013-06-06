#!/usr/bin/python

import os
import nsmod_cython
import nsmod_one_component_model
import nsmod_two_component_model
import nsmod_two_component_model2
import pynotify


def Run(Input_Dictionary):
    """ Create a generic write file, compile and run in C.

    Keyword arguments must be passed as a dictionary
    chi:float [degrees] -- angle between magnetic dipole and z axis
    epsI:float []-- elastic deformation
    epsA :float []-- magnetic deformation
    omega0 :float [Hz]-- initial spin period
    either
        eta :float []-- Threshold for which to stop the simulation.
             Simulation is initiated using omega**2.0 < eta*omega_0**2.0
    or
        t1:float [s] -- Time for which to run simulation for

    Optional arguments
    err : float[] -- Error argument to pass to the GCC compiler

    For help with the GCC compiler see documentation at
    http://www.gnu.org/software/gsl/manual/
    html_node/Ordinary-Differential-Equations.html
    The generic script is written by Write_File in File_Functions
    """

    #  Required paramaters

    file_name_list = []

    if Input_Dictionary.get('no_anom'):
        file_name_list .append("no_anom_")
        anom_torque = False
    else:
        anom_torque = True

    try:
        chi_degrees = Input_Dictionary['chi']
        file_name_list.append("chi_" + str(Input_Dictionary['chi']))
    except KeyError:
        print " ERROR: You need to specify chi in the input dictionary"
        return

    # Triaxial test
    try:
        epsI1 = str(Input_Dictionary['epsI1'])
        epsI3 = str(Input_Dictionary['epsI3'])
        file_name_list.append("_epsI1_{}".format(Input_Dictionary['epsI1']))
        file_name_list.append("_epsI3_{}".format(Input_Dictionary['epsI3']))
    except KeyError:
        # If not trixial test for biaxial
        try:
            epsI1 = str(Input_Dictionary['epsI'])
            epsI3 = "0"
            file_name_list.append("_epsI_{}".format(Input_Dictionary['epsI']))
        except KeyError:
            print (" ERROR: You need to specify epsI1 and epsI2,"
                  " or epsI (Biaxial case) in the input dictionary")
            return

    try:
        epsA = str(Input_Dictionary['epsA'])
        file_name_list.append("_epsA_" + str(Input_Dictionary['epsA']))
    except KeyError:
        print " ERROR: You need to specify epsA in the input dictionary"
        return

    try:
        omega0 = str(Input_Dictionary['omega0'])
        file_name_list.append("_omega0_" + str(Input_Dictionary['omega0']))
    except KeyError:
        print " ERROR: You need to specify omega0 in the input dictionary"
        return

    try:
        eta = str(Input_Dictionary['eta'])
        file_name_list.append("_eta_" + str(Input_Dictionary['eta']))
        eta_relative = str(float(eta) * pow(float(omega0), 2))
    except KeyError:
        eta = 0.0

    try:
        t1 = str(Input_Dictionary['t1'])
        file_name_list.append("_t1_" + str(Input_Dictionary['t1']))
    except KeyError:
        print "ERROR: t1 not specified using a default "
        t1 = 1e15

#   Additional Arguments
    try:
        error = str(Input_Dictionary['error'])
    except KeyError:
        # Use a default value
        error = 1e-12

    try:
        n = int(Input_Dictionary['n'])
    except KeyError:
        # Don't save at n discrete intervals...
        n = None

    # Create file name
    file_name_list.append(".hdf5")
    file_name = "".join(file_name_list)

    # Tempory hack to remove old files seems to cause issues
    if file_name in os.listdir("."):
        os.remove(file_name)

    nsmod_one_component_model.main(
                      chi_degrees=float(chi_degrees),
                      file_name=file_name,
                      n=n,
                      epsA=float(epsA),
                      epsI1=float(epsI1),
                      epsI3=float(epsI3),
                      omega0=float(omega0),
                      t1=float(t1),
                      anom_torque=anom_torque,
                      error=float(error)
                      )

    pynotify.init("Basic")

    n = pynotify.Notification("Run 1CM complete",
      "output saved as {}".format(file_name)
    )

    n.show()
    return file_name


def Run_Cython_Two_Component(Input_Dictionary):
    """ Function to run the two component model.

    Keyword arguments must be passed as a dictionary
    chi:float [degrees] -- angle between magnetic dipole and z axis
    epsI:float []-- elastic deformation
    epsA :float []-- magnetic deformation
    omega0 :float [Hz]-- initial spin period
    t1:float [s] -- Time for which to run simulation for
    eta: float
    K:float -- Coupling constant between the two components
    Ishell:MOI of the biaxial crust
    Icore:MOI of the spherical core

    Optional arguments
    err : float[] -- Error argument to pass to the GCC compiler

    For help with the GCC compiler see documentation at
    http://www.gnu.org/software/gsl/manual/html_node/
    Ordinary-Differential-Equations.html
    """

    #  Required paramaters

    # Initiate the file_name list
    file_name_list = []

    # Check the anom_torque condition
    if Input_Dictionary.get('no_anom'):
        file_name_list .append("no_anom_")
        anom_torque = False
    else:
        anom_torque = True

    # At the moment these errors do not raise correctly.

    try:
        chi_degrees = Input_Dictionary['chi']
        file_name_list.append("chi_" + str(Input_Dictionary['chi']))
    except KeyError:
        print " ERROR: You need to specify chi in the input dictionary"
        return

    try:
        epsI = str(Input_Dictionary['epsI'])
        file_name_list.append("_epsI_" + str(Input_Dictionary['epsI']))
    except KeyError:
        print " ERROR: You need to specify epsI in the input dictionary"
        return

    try:
        epsA = str(Input_Dictionary['epsA'])
        file_name_list.append("_epsA_" + str(Input_Dictionary['epsA']))
    except KeyError:
        print " ERROR: You need to specify epsA in the input dictionary"
        return

    try:
        omega0 = str(Input_Dictionary['omega0'])
        file_name_list.append("_omega0_" + str(Input_Dictionary['omega0']))
    except KeyError:
        print " ERROR: You need to specify omega0 in the input dictionary"
        return

    try:
        t1 = str(Input_Dictionary['t1'])
        file_name_list.append("_t1_" + str(Input_Dictionary['t1']))
    except KeyError:
        t1 = "1e15"
        print "ERROR: You have not yet specified t1, \
                using default value t1 = {}".format(t1)

    try:
        eta = str(Input_Dictionary['eta'])
        file_name_list.append("_eta_" + str(Input_Dictionary['eta']))
        eta_relative = str(float(eta) * pow(float(omega0), 2))
    except KeyError:
        eta = "0.0"
        print "ERROR: You have not yet specified eta, \
         using default values  eta ={}".format(eta)

    try:
        K = str(Input_Dictionary['K'])
        file_name_list.append("_K_" + str(Input_Dictionary['K']))
    except KeyError:
        print "ERROR: You must specify K"
        return

    try:
        Ishell = str(Input_Dictionary['Ishell'])
        file_name_list.append("_Ishell_" + str(Input_Dictionary['Ishell']))
    except KeyError:
        Ishell = 1e45
        print "ERROR: Ishell not specified \
            using default Ishell={}".format(Ishell)

    try:
        Icore = str(Input_Dictionary['Icore'])
        file_name_list.append("_Icore_" + str(Input_Dictionary['Icore']))
    except KeyError:
        Icore = 1e45
        print "ERROR: Icore not specified using default Icore={}".format(Icore)

    try:
        aw_int = str(Input_Dictionary['aw_int'])
        file_name_list.append("_aw_int_" + str(Input_Dictionary['aw_int']))
    except KeyError:
        aw_int = "50.0"
        print ("ERROR: aw_int not specified"
              " using default aw_int={}".format(aw_int))

    try:
        aW_int = str(Input_Dictionary['aW_int'])
        file_name_list.append("_aW_int_" + str(Input_Dictionary['aW_int']))
    except KeyError:
        aW_int = "50.0"
        print ("ERROR: aW_int not specified"
              " using default aW_int={}".format(aw_int))

#   Additional Arguments
    try:
        error = str(Input_Dictionary['error'])
    except KeyError:
        # Use a default value
        error = 1e-5

    try:
        n = int(Input_Dictionary['n'])
        #if "_t1_" not in file_name_list:
            #print "Warning: Specifying n
            #without t1 means you may get a very course save "
    except KeyError:
        # Don't save at n discrete intervals...
        n = None

    # Create file name
    file_name_list.append(".hdf5")
    file_name = "".join(file_name_list)

    # Tempory hack to remove old files seems to cause issues
    if file_name in os.listdir("."):
        os.remove(file_name)

    nsmod_two_component_model.main(
        chi_degrees=float(chi_degrees),
        file_name=file_name,
        epsA=float(epsA),
        epsI=float(epsI),
        n=n,
        omega0=float(omega0),
        t1=float(t1),
        anom_torque=anom_torque,
        error=float(error),
        Ishell=float(Ishell),
        Icore=float(Icore),
        K=float(K),
        aw_int=float(aw_int),
        aW_int=float(aW_int)
        )

    pynotify.init("Basic")

    n = pynotify.Notification("Run 2CM complete",
      "output saved as {}".format(file_name)
    )

    n.show()
    return file_name


def Run_Cython_Two_Component2(Input_Dictionary):
    """ Function to run the two component model.

    Keyword arguments must be passed as a dictionary
    chi:float [degrees] -- angle between magnetic dipole and z axis
    epsI:float []-- elastic deformation
    epsA :float []-- magnetic deformation
    omega0 :float [Hz]-- initial spin period
    t1:float [s] -- Time for which to run simulation for
    eta: float
    K:float -- Coupling constant between the two components
    Ishell:MOI of the biaxial crust
    Icore:MOI of the spherical core

    Optional arguments
    err : float[] -- Error argument to pass to the GCC compiler

    For help with the GCC compiler see documentation at
    http://www.gnu.org/software/gsl/manual/html_node/
    Ordinary-Differential-Equations.html
    """

    #  Required paramaters

    # Initiate the file_name list
    file_name_list = []

    # Check the anom_torque condition
    if Input_Dictionary.get('no_anom'):
        file_name_list .append("no_anom_")
        anom_torque = False
    else:
        anom_torque = True

    # At the moment these errors do not raise correctly.

    try:
        chi_degrees = Input_Dictionary['chi']
        file_name_list.append("chi_" + str(Input_Dictionary['chi']))
    except KeyError:
        print " ERROR: You need to specify chi in the input dictionary"
        return

    try:
        epsI = str(Input_Dictionary['epsI'])
        file_name_list.append("_epsI_" + str(Input_Dictionary['epsI']))
    except KeyError:
        print " ERROR: You need to specify epsI in the input dictionary"
        return

    try:
        epsA = str(Input_Dictionary['epsA'])
        file_name_list.append("_epsA_" + str(Input_Dictionary['epsA']))
    except KeyError:
        print " ERROR: You need to specify epsA in the input dictionary"
        return

    try:
        t1 = str(Input_Dictionary['t1'])
        file_name_list.append("_t1_" + str(Input_Dictionary['t1']))
    except KeyError:
        t1 = "1e15"
        print "ERROR: You have not yet specified t1, \
                using default value t1 = {}".format(t1)

    try:
        K = str(Input_Dictionary['K'])
        file_name_list.append("_K_" + str(Input_Dictionary['K']))
    except KeyError:
        print "ERROR: You must specify K"
        return

    try:
        Ishell = str(Input_Dictionary['Ishell'])
        file_name_list.append("_Ishell_" + str(Input_Dictionary['Ishell']))
    except KeyError:
        Ishell = 1e45
        print "ERROR: Ishell not specified \
            using default Ishell={}".format(Ishell)

    try:
        Icore = str(Input_Dictionary['Icore'])
        file_name_list.append("_Icore_" + str(Input_Dictionary['Icore']))
    except KeyError:
        Icore = 1e45
        print "ERROR: Icore not specified using default Icore={}".format(Icore)

    # Ugly hack
    core_int_str = Input_Dictionary['core_int']
    core_int = [float(p) for p in
        core_int_str.lstrip("[").rstrip("]").split(",")]

    shell_int_str = Input_Dictionary['shell_int']
    shell_int = [float(p) for p in
        shell_int_str.lstrip("[").rstrip("]").split(",")]

#   Additional Arguments
    try:
        error = str(Input_Dictionary['error'])
    except KeyError:
        # Use a default value
        error = 1e-5

    try:
        n = int(Input_Dictionary['n'])
        #if "_t1_" not in file_name_list:
            #print "Warning: Specifying n
            #without t1 means you may get a very course save "
    except KeyError:
        # Don't save at n discrete intervals...
        n = None

    # Create file name
    file_name_list.append(".hdf5")
    file_name = "".join(file_name_list)

    # Tempory hack to remove old files seems to cause issues
    if file_name in os.listdir("."):
        os.remove(file_name)

    nsmod_two_component_model2.main(
        chi_degrees=float(chi_degrees),
        file_name=file_name,
        epsA=float(epsA),
        epsI=float(epsI),
        n=n,
        t1=float(t1),
        error=float(error),
        Ishell=float(Ishell),
        Icore=float(Icore),
        K=float(K),
        core_int=core_int,
        shell_int=shell_int
        )

    pynotify.init("Basic")

    n = pynotify.Notification("Run 2CM complete",
      "output saved as {}".format(file_name)
    )

    n.show()
    return file_name




