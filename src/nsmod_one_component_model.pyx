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

import numpy as np
import math
from cython_gsl cimport *
import h5py

ctypedef struct vector:
    double one
    double two
    double three

cdef double Torque_over_Io (double w[], void *params) nogil:
    """ Returns the Goldreich torque

    Note: it is important that chi is given in radians and not degrees,
    this is not checked.

    """

    cdef double chi, epsA, anom_torque, pre, mx, mz, T_o_I[3]

    epsA = (<double *> params)[0]
    chi = (<double *> params)[1]
    anom_torque = (<double *> params)[4]

    mx = sin(chi)
    mz = cos(chi)

    pre =  ( pow(w[0],2) + pow(w[1],2) + pow(w[2],2) ) * 2.0 * pow(10,6) * pow(9.0 * pow(10,10),-1)

    Tx_sd = pre * epsA * mz * (w[2] * mx - w[0] * mz)
    Ty_sd = -pre * epsA * w[2]
    Tz_sd = pre * epsA * mx * (w[0] * mz - w[2] * mx)

    #if anom_torque == 1:
    Tx_an = epsA * (w[0] * mx + w[2] * mz) * w[1] * mz
    Ty_an = epsA * (w[0] * mx + w[2] * mz) * (w[2] * mx - w[0] * mz)
    Tz_an = -epsA * (w[0] * mx + w[2] * mz) * w[1] * mx

    #T_o_I = {Tx_sd + Tx_an, Tx_sd + Tx_an, Tx_sd + Tx_an}

    #T_o_I[1] = Ty_sd + Ty_an
    #T_o_I[2] = Tz_sd + Tz_an

    #else

    #T_o_I[0] = Tx_sd
    # T_o_I[1] = Ty_sd
    #T_o_I[2] = Tz_sd

    return Tx_sd


cdef int funcs (double t, double w[], double f[], void *params) nogil:
    """ Function defining the ODEs with the anomalous torque """
    # Define the variables used in the calculation
    cdef double epsI1,epsI3

    # Import the constant variables from params
    epsA = (<double *> params)[0]
    chi = (<double *> params)[1]
    epsI1 = (<double *> params)[2]
    epsI3 = (<double *> params)[3]
    anom_torque = (<double *> params)[4]

    # Calculate the torque
    mx = sin(chi)
    mz = cos(chi)

    pre =  ( pow(w[0],2) + pow(w[1],2) + pow(w[2],2) ) * 2.0 * pow(10,6) * pow(9.0 * pow(10,10),-1)

    Tx_sd = pre * epsA * mz * (w[2] * mx - w[0] * mz)
    Ty_sd = -pre * epsA * w[2]
    Tz_sd = pre * epsA * mx * (w[0] * mz - w[2] * mx)

    if anom_torque == 1:
        Tx = Tx_sd + epsA * (w[0] * mx + w[2] * mz) * w[1] * mz
        Ty = Ty_sd + epsA * (w[0] * mx + w[2] * mz) * (w[2] * mx - w[0] * mz)
        Tz = Tz_sd -epsA * (w[0] * mx + w[2] * mz) * w[1] * mx

    else:
        Tx = Tx_sd
        Ty = Ty_sd
        Tz = Tz_sd


    #  Define the three ODEs in f[] as functions of the above variables

    f[0] = Tx * pow(1 + epsI1,-1) - w[1] * w[2] * epsI3 * pow(1 + epsI1, -1)

    f[1] = Ty - w[0] * w[2] * (epsI1 - epsI3)

    f[2] = Tz * pow(1 + epsI3,-1) + w[0] * w[1] * epsI1 * pow(1 + epsI3, -1)

    return GSL_SUCCESS


cdef int no_anom_torque (double t, double y[], double f[], void *params) nogil:
    """ Function defining the ODEs with the anomalous torque """
    # Define the variables used in the calculation
    cdef double wx,wy,wz,w_2,epsI,epsA,chi,Lambda

    # Import the three time dependant variables from y[]
    wx=y[0]
    wy=y[1]
    wz=y[2]
    w_2 = pow(wx,2)+pow(wy,2)+pow(wz,2)

    # Import the constant variables from params
    Lambda = (<double *> params)[0]    #= 2R/3C */
    epsI= (<double *> params)[1]
    epsA = (<double *> params)[2]
    chi = (<double *> params)[3]

    # Calculate the angular parts of the three equations to avoid repeated calculation
    cdef double Cx,Sx
    Cx = cos(chi)
    Sx = sin(chi)

    #  Define the three ODEs in f[] as functions of the above variables
    f[0] = epsA*(Lambda*w_2*Cx*(wz*Sx-wx*Cx)) - wy*wz*epsI

    f[1] = epsA*(-Lambda*w_2*wy) + wx*wz*epsI

    f[2] = epsA*pow(1+epsI,-1) * (Lambda*w_2*Sx*(wx*Cx - wz*Sx))
    return GSL_SUCCESS

# Currently jac is unused by the ODE solver so is left empty
cdef int jac (double t, double y[], double *dfdy, double dfdt[], void *params) nogil:

    return GSL_SUCCESS


def main (epsI1=-1.0e-6, epsI3=1.0e-6, epsA=1.0e-8 , omega0=1.0e1,
    error=1e-10, t1=1.0e3 , eta=0.0 ,chi_degrees = 30.0 ,anom_torque=True ,
    a_int=0.8736, file_name="generic.hdf5",n=None):
    """ Solve the one component model  using gsl_odeiv2_step_rk8pd """

    # Define default variables
    cdef double chi

   # Test if the anomalous torque is required or not
    if anom_torque:
        anom_torque = 1
    else:
        anom_torque = 0

    # We allow the user to give chi in degrees and convert here
    chi=math.pi*chi_degrees/180.0

    # Pass them to params list
    cdef double params[5]
    params[0] = epsA
    params[1] = chi
    params[2] = epsI1
    params[3] = epsI3
    params[4] = anom_torque

    # Initial values and calculate eta_relative
    cdef int i
    cdef double t , w[3] ,h , eta_relative
    eta_relative = eta*pow(omega0,2)
    h = 1e-15   # Initial step size
    t = 0.0
    w[0] = omega0*sin(a_int)
    w[1] = 0.0
    w[2] = omega0*cos(a_int)

    # Inititate the system and define the set of functions
    cdef gsl_odeiv2_system sys

    sys.function = funcs
    sys.jacobian = jac
    sys.dimension = 3
    sys.params = params

    # Setup the solver ~ Note not all of these cdefs are always used.
    # It seems cython won't accept a cdef inside an if statement.

    # For n
    cdef gsl_odeiv2_driver * d
    d = gsl_odeiv2_driver_alloc_y_new(
    &sys, gsl_odeiv2_step_rk8pd,
    error, error, 0.0)

    # Setup the solver alternative to n
    cdef gsl_odeiv2_step_type * T
    T = gsl_odeiv2_step_rk8pd

    cdef gsl_odeiv2_step * s
    s = gsl_odeiv2_step_alloc (T, 3)
    cdef gsl_odeiv2_control * c
    c = gsl_odeiv2_control_y_new (error, error)
    cdef gsl_odeiv2_evolve * e
    e = gsl_odeiv2_evolve_alloc (3)

    cdef int status
    cdef double ti

    time = []
    w1=[]
    w2=[]
    w3=[]

    if n :
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

        gsl_odeiv2_driver_free(d)

    else:

        while (pow(w[0],2)+pow(w[1],2)+pow(w[2],2) > eta_relative and t < t1 ):
            status = gsl_odeiv2_evolve_apply (e, c, s, &sys, &t, t1, &h, w)

            if (status != GSL_SUCCESS):
                break

            time.append(t)
            w1.append(w[0])
            w2.append(w[1])
            w3.append(w[2])

        gsl_odeiv2_evolve_free (e)
        gsl_odeiv2_control_free (c)
        gsl_odeiv2_step_free (s)

    f = h5py.File(file_name,'w')
    f.create_dataset("time",data=time)
    f.create_dataset("w1",data=w1)
    f.create_dataset("w2",data=w2)
    f.create_dataset("w3",data=w3)
    return file_name
