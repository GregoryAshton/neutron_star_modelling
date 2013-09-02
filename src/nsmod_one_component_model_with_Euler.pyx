"""

Solver for three coupled ODEs which define the single component system

The solver used GSL odeint to actually solve the odes and Cython_GSL
as a npthon interface documentation can be found at

    http://www.gnu.org/software/gsl/manual/
    html_node/Ordinary-Differential-Equations.html

and

    https://github.com/twiecki/CythonGSL

any changes to these files will require recompiling using the setup.np in
the home directory.

"""

from numpy import arccos
import numpy as np
import math
from cython_gsl cimport *
import h5py


cdef int funcs (double t, double w[], double f[], void *params) nogil:
    """ Function defining the ODEs with the anomalous torque """

    # Import the constant variables from params
    epsA = (<double *> params)[0]
    chi = (<double *> params)[1]
    epsI1 = (<double *> params)[2]
    epsI3 = (<double *> params)[3]
    anom_torque_b = (<double *> params)[4]

    # Calculate the torque
    mx = sin(chi)
    mz = cos(chi)

    pre =  ( pow(w[0],2) + pow(w[1],2) + pow(w[2],2) ) * 2.0 * pow(10,6) * pow(9.0 * pow(10,10),-1)

    Tx_sd = pre * epsA * mz * (w[2] * mx - w[0] * mz)
    Ty_sd = -pre * epsA * w[1]
    Tz_sd = pre * epsA * mx * (w[0] * mz - w[2] * mx)

    if anom_torque_b == 1:
        Tx = Tx_sd + epsA * (w[0] * mx + w[2] * mz) * w[1] * mz
        Ty = Ty_sd + epsA * (w[0] * mx + w[2] * mz) * (w[2] * mx - w[0] * mz)
        Tz = Tz_sd - epsA * (w[0] * mx + w[2] * mz) * w[1] * mx

    else:
        Tx = Tx_sd
        Ty = Ty_sd
        Tz = Tz_sd


    #  Define the three ODEs in f[] as functions of the above variables

    #f[0] = Tx * pow(1 + epsI1,-1) - w[1] * w[2] * epsI3 * pow(1 + epsI1, -1)

    #f[1] = Ty - w[0] * w[2] * (epsI1 - epsI3)

    #f[2] = Tz * pow(1 + epsI3,-1) + w[0] * w[1] * epsI1 * pow(1 + epsI3, -1)
    
    #f[3] = w[0] * cos(w[5]) - w[1] * sin(w[5])
    
    #f[4] = pow(sin(w[3]), -1) * (w[0] * sin(w[5]) + w[1] * cos(w[5]))
 
    #f[5] = w[2] - pow(sin(w[3]), -1) * (w[0] * sin(w[5]) + w[1] * cos(w[5])) * cos(w[3])
    
    f[0] = Tx - w[1] * w[2] * epsI3

    f[1] = Ty + w[0] * w[2] * epsI3

    f[2] = Tz
    
    f[3] = w[0] * cos(w[5]) - w[1] * sin(w[5])
    
    f[4] = pow(sin(w[3]), -1) * (w[0] * sin(w[5]) + w[1] * cos(w[5]))
 
    f[5] = w[2] - pow(sin(w[3]), -1) * (w[0] * sin(w[5]) + w[1] * cos(w[5])) * cos(w[3])

    return GSL_SUCCESS


# Currently jac is unused by the ODE solver so is left empty
cdef int jac (double t, double y[], double *dfdy, double dfdt[], void *params) nogil:

    return GSL_SUCCESS


def main (epsI1=0.0, epsI3=1.0e-6, epsA=1.0e-8 , omega0=1.0e1,
    error=1e-10, t1=1.0e3 ,  chi = 30.0, anom_torque=True ,
    a_int=50.0, file_name="generic.hdf5", n=None):
    
    """ Solve the one component model  using gsl_odeiv2_step_rk8pd """

    if epsI1 != 0.0:
        print "Triaxiality is currently not implemented, you need to define the initial conditions correctly"
        return

   # Test if the anomalous torque is required or not
    if anom_torque:
        anom_torque_b = 1
    else:
        anom_torque_b = 0

    # We allow the user to give chi in degrees and convert here
    chi = np.deg2rad(chi)
    a_int = np.deg2rad(a_int)

    # Pass them to params list
    cdef double params[5]
    params[0] = epsA
    params[1] = chi
    params[2] = epsI1
    params[3] = epsI3
    params[4] = anom_torque_b

    # Initial values and calculate eta_relative
    cdef int i
    cdef double t , w[6] ,h 
    h = 1e-15   # Initial step size
    t = 0.0
    w[0] = omega0*sin(a_int)
    w[1] = 0.0
    w[2] = omega0*cos(a_int)
    w[3] = np.arccos(cos(a_int) * (1.0 + epsI3) * pow(pow(sin(a_int) * (1.0 + epsI1), 2) + pow(cos(a_int) * (1.0 + epsI3), 2), -0.5))
    w[4] = 0.0
    w[5] = np.arccos(sin(a_int) * (1.0 + epsI1) * pow(1 - pow(cos(a_int) * (1.0 + epsI3), 2), -0.5))


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
                    1e-10, # const double hstart
                    error,  #const double epsabs,
                    error # const double epsrel
                    )

    cdef int status
    cdef double ti

    time = [0.0]
    w1 = [w[0]]
    w2 = [w[1]]
    w3 = [w[2]]
    w4 = [w[3]]
    w5 = [w[4]]
    w6 = [w[5]]


    # Run saving at discrete time values
    for i from 1 <= i <= n:
        ti = i * t1 / n
        status = gsl_odeiv2_driver_apply (d, &t, ti, w)

        if (status != GSL_SUCCESS):
            print("error, return value=%d\n" % status)
            break

        time.append(t)
        w1.append(w[0])
        w2.append(w[1])
        w3.append(w[2])
        w4.append(w[3])
        w5.append(w[4])
        w6.append(w[5])

    gsl_odeiv2_driver_free(d)

    f = h5py.File(file_name, 'w')
    f.create_dataset("time", data=time)
    f.create_dataset("w1", data=w1)
    f.create_dataset("w2", data=w2)
    f.create_dataset("w3", data=w3)
    f.create_dataset("theta", data=w4)
    f.create_dataset("phi", data=w5)
    f.create_dataset("psi", data=w6)
    return file_name
