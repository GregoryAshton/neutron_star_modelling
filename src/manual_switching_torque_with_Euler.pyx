"""
Solver for three coupled ODEs which define the single component system along 
with the Euler angle equations

The solver uses GSL odeint to actually solve the odes and Cython_GSL
as a python interface, for help see documentation:

    gnu.org/software/gsl/manual/html_node/Ordinary-Differential-Equations.html

and

    github.com/twiecki/CythonGSL

any changes to these files will require recompiling using the setup.np in
the home directory.

"""

import numpy as np
from cython_gsl cimport *
import h5py
from nsmod.File_Functions import FileNamer


cdef int funcs (double t, double w[], double f[], void *params) nogil:
    """ Function defining the ODEs with the anomalous torque """
    # Define the variables used in the calculation
    cdef double S, Psi

    # Import the constant variables from params
    epsA = (<double *> params)[0]
    chi = (<double *> params)[1]
    epsI1 = (<double *> params)[2]
    epsI3 = (<double *> params)[3]
    AnomTorque = (<double *> params)[4]
    upsilon = (<double *> params)[5]
    SwitchTime = (<double *> params)[6]
    dt = (<double *> params)[7]

    # Calculate the torque
    mx = sin(chi)
    mz = cos(chi)

    pre =  ((pow(w[0], 2) + pow(w[1], 2) + pow(w[2], 2)) * 
             2.0 * pow(10, 6) * pow(9.0 * pow(10,10),-1))

    Tx_sd = pre * epsA * mz * (w[2] * mx - w[0] * mz)
    Ty_sd = -pre * epsA * w[1]
    Tz_sd = pre * epsA * mx * (w[0] * mz - w[2] * mx)

    S = 1.0 

    if t > SwitchTime:
        S = 1.0 - upsilon
    
    
    if AnomTorque == 1:
        Tx = S * (Tx_sd + epsA * (w[0] * mx + w[2] * mz) * w[1] * mz)
        Ty = S * (Ty_sd + epsA * (w[0] * mx + w[2] * mz) * (w[2] * mx - w[0] * mz))
        Tz = S * (Tz_sd - epsA * (w[0] * mx + w[2] * mz) * w[1] * mx)

    else:
        Tx = S * Tx_sd
        Ty = S * Ty_sd
        Tz = S * Tz_sd


    #  Define the three ODEs in f[] as functions of the above variables

    f[0] = Tx * pow(1 + epsI1, -1) - w[1] * w[2] * epsI3 * pow(1 + epsI1, -1)

    f[1] = Ty + w[0] * w[2] * (epsI3 - epsI1)

    f[2] = Tz * pow(1 + epsI3, -1) + w[0] * w[1] * epsI1 * pow(1 + epsI3, -1)
    
    f[3] = w[0] * cos(w[5]) - w[1] * sin(w[5])
    
    f[4] = pow(sin(w[3]), -1) * (w[0] * sin(w[5]) + w[1] * cos(w[5]))

    f[5] = w[2] - pow(sin(w[3]), -1) * (w[0] * sin(w[5]) + w[1] * cos(w[5])) * cos(w[3])
    

    return GSL_SUCCESS


# Currently jac is unused by the ODE solver so is left empty
cdef int jac (double t, double y[], double *dfdy, 
              double dfdt[], void *params) nogil:

    return GSL_SUCCESS


def main (epsI1=0.0, epsI3=1.0e-6, epsA=1.0e-8 , omega0=1.0e1, chi0=30.0,
    a0=50., T=1.0e3, AnomTorque=True , upsilon=0.0, n=10000, error=1e-10,
    cleanup=False, SwitchTime=100, DryRun=False):
    """ One component NS with Euler angles and switching
    
    This solves the Euler equations for a single component NS and the 
    Euler angles to take it into the inertial frame. The body is acted
    on by the Deutsch torque, with the addition of a switching component

    Paramaters
    ----------
    epsI1 : float
        Ellipticity along the x axis
    epsI3 : float
        Ellipticity along the z axis
    epsA : float
        Magnetic deformation [Glampedakis & Jones, 2010]
    omega0 : float
        Initial magnitude of the spin vector
    chi0 : float
        Initial polar angle of the magnetic dipole in degrees
    a0 : float
        Initial polar angle of the spin vector in degrees
    T : float
        Duration of the simulation in seconds
    AnomTorque : bool
        If true, include the anomalous torque
    upsilon : float in [0, 1]
        Fractional reduction in spindown torque during the off-phase   
    n : int
        Number of data points to save
    error : float
        Error passed to the ODE solver
    cleanup = bool
        If true old data files with the same file_name will be removed. If 
        False then the simulation will not run and simply return the file_name
    
    
    """
 
    (file_name, run_sim) = FileNamer(epsI1=epsI1, epsI3=epsI3, epsA=epsA,
                          omega0=omega0, chi0=chi0, a0=a0, T=T, SwitchTime=SwitchTime,
                          AnomTorque=AnomTorque, n=n, upsilon=upsilon,
                          error=error, cleanup=cleanup)
    if not run_sim or DryRun:
        return file_name

   # Convert python Bool to int
    AnomTorque = AnomTorque.real

    # We allow the user to give angles in degrees and convert here
    chi0 = np.deg2rad(chi0)
    a0 = np.deg2rad(a0)

    cdef double ti, dt
    dt = float(T) / n

    # Pass them to params list
    cdef double params[8]
    params[0] = epsA
    params[1] = chi0
    params[2] = epsI1
    params[3] = epsI3
    params[4] = AnomTorque
    params[5] = upsilon
    params[6] = SwitchTime
    params[7] = dt

    # Initial values and calculate eta_relative
    cdef int i
    cdef double t , w[6]
    t = 0.0
    w[0] = omega0*sin(a0)
    w[1] = 0.0
    w[2] = omega0*cos(a0)
    w[3] = (np.arccos(cos(a0) * (1.0 + epsI3) / 
             np.sqrt(pow(sin(a0) * (1.0 + epsI1), 2) + 
                 pow(cos(a0) * (1.0 + epsI3), 2))))
    w[4] = 0.0
    w[5] = 0.5 * np.pi 

    # Inititate the system and define the set of functions
    cdef gsl_odeiv2_system sys

    sys.function = funcs
    sys.jacobian = jac
    sys.dimension = 6
    sys.params = params

    # Setup the solver 
    cdef gsl_odeiv2_driver * d
    d = gsl_odeiv2_driver_alloc_y_new(
                    &sys,  # const gsl_odeiv2_system * sys
                    gsl_odeiv2_step_rkf45, # const gsl_odeiv2_step_type * T
                    1e-15, # const double hstart
                    error,  #const double epsabs,
                    error # const double epsrel
                    )

    cdef int status


    w1 = [w[0]]
    w2 = [w[1]]
    w3 = [w[2]]
    w4 = [w[3]]
    w5 = [w[4]]
    w6 = [w[5]]
    
    
    # Run saving at discrete time values
    for i from 1 <= i <= n:
        ti = i * dt
        status = gsl_odeiv2_driver_apply (d, &t, ti, w)

        if (status != GSL_SUCCESS):
            print("error, return value=%d\n" % status)
            break
        
        w1.append(w[0])
        w2.append(w[1])
        w3.append(w[2])
        w4.append(w[3])
        w5.append(w[4])
        w6.append(w[5])

    gsl_odeiv2_driver_free(d)

    time = np.linspace(0, T, n+1)

    f = h5py.File(file_name, 'w')
    f.create_dataset("time", data=time)
    f.create_dataset("w1", data=w1)
    f.create_dataset("w2", data=w2)
    f.create_dataset("w3", data=w3)
    f.create_dataset("theta", data=w4)
    f.create_dataset("phi", data=w5)
    f.create_dataset("psi", data=w6)
    f.close()
    return file_name
