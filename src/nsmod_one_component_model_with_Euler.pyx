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

#!python
#cython: boundscheck=False
#cython: wraparound=False 
#cython: nonecheck=False

from cython_gsl cimport *
import h5py
import numpy as np


cdef int funcs (double t, double w[], double f[], void *params) nogil:
    """ Function defining the ODEs with the anomalous torque """

    cdef double mx, mz, Tx_sd, Ty_sd, Tz_sd, sin_theta, cos_theta, sin_psi, cos_psi, epsA, chi_rad, epsI1, epsI3
    
    # Import the constant variables from params
    epsA = (<double *> params)[0]
    chi_rad = (<double *> params)[1]
    epsI1 = (<double *> params)[2]
    epsI3 = (<double *> params)[3]
    anom_torque_b = (<double *> params)[4]

    # Precalculations
    mx = sin(chi_rad)
    mz = cos(chi_rad)

    pre =  epsA * (pow(w[0],2) + pow(w[1],2) + pow(w[2],2)) * pow(45000, -1)

    Tx_sd = pre * mz * (w[2] * mx - w[0] * mz)
    Ty_sd = -pre * w[1]
    Tz_sd = pre * mx * (w[0] * mz - w[2] * mx)
    
    if anom_torque_b == 1:
        Tx = Tx_sd + epsA * (w[0] * mx + w[2] * mz) * w[1] * mz
        Ty = Ty_sd + epsA * (w[0] * mx + w[2] * mz) * (w[2] * mx - w[0] * mz)
        Tz = Tz_sd - epsA * (w[0] * mx + w[2] * mz) * w[1] * mx

    else:
        Tx = Tx_sd
        Ty = Ty_sd
        Tz = Tz_sd

    #  Define the three ODEs in f[] as functions of the above variables

    sin_theta = sin(w[4])
    cos_theta = cos(w[4])
    sin_psi = sin(w[5])
    cos_psi = cos(w[5])
    
    f[0] = pow(1 + epsI1, -1) * (Tx  - w[1] * w[2] * epsI3) 

    f[1] = Ty + w[0] * w[2] * (epsI3 - epsI1)

    f[2] = pow(1 + epsI3, -1) * (Tz + w[0] * w[1] * epsI1 )
    
    f[3] = pow(sin_theta, -1) * (w[0] * sin_psi + w[1] * cos_psi)
    
    f[4] = (w[0] - f[3] * sin_theta * sin_psi) * pow(cos_psi, -1)  #cos_psi *(w[0] - w[1] * tan(w[5])) 
    
    f[5] = w[2] - f[3] * cos_theta
    
    return GSL_SUCCESS


# Currently jac is unused by the ODE solver so is left empty
cdef int jac (double t, double y[], double *dfdy, double dfdt[], void *params) nogil:

    return GSL_SUCCESS


def main(epsI1=0.0, epsI3=1.0e-6, epsA=1.0e-8 , omega0=1.0e1,
    error=1e-10, t1=1.0e3 ,  chi = 30.0, anom_torque=True ,
    a_int=50.0, file_name="generic.hdf5", n=None):
    
    """ Solve the one component model  using gsl_odeiv2_step_rk8pd """

    cdef long double chi_rad, a_int_rad
   # Test if the anomalous torque is required or not
    if anom_torque:
        anom_torque_b = 1
    else:
        anom_torque_b = 0

    # We allow the user to give chi in degrees and convert here
    chi_rad = chi * M_PI / 180.0
    a_int_rad = a_int * M_PI / 180.0

    # Pass them to params list
    cdef double params[5]
    params[0] = epsA
    params[1] = chi_rad
    params[2] = epsI1
    params[3] = epsI3
    params[4] = anom_torque_b

    # Initial values and calculate eta_relative
    cdef int i
    cdef double t , w[6] 
    t = 0.0
    w[0] = omega0*sin(a_int_rad)
    w[1] = 0.0
    w[2] = omega0*cos(a_int_rad)
    w[3] = 0.0
    w[4] = acos(cos(a_int_rad) * (1.0 + epsI3) / sqrt(pow(sin(a_int_rad) * (1.0 + epsI1), 2) + pow(cos(a_int_rad) * (1.0 + epsI3), 2)))
    w[5] = 0.5 * M_PI


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
    cdef double ti, dt

    w1 = [w[0]]
    w2 = [w[1]]
    w3 = [w[2]]
    w4 = [w[3]]
    w5 = [w[4]]
    w6 = [w[5]]

    dt = t1 / n
    
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

    time = np.linspace(0, t1, n+1)
    f = h5py.File(file_name, 'w')
    f.create_dataset("time", data=time)
    f.create_dataset("w1", data=w1)
    f.create_dataset("w2", data=w2)
    f.create_dataset("w3", data=w3)
    f.create_dataset("phi", data=w4)
    f.create_dataset("theta", data=w5)
    f.create_dataset("psi", data=w6)
    return file_name
