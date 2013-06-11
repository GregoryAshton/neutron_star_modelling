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

from cython_gsl cimport *
import math
import numnp as np
import h5py

# Need to add in anom_torque

cdef Torque_over_Io(double y[], void *params):
    """ Returns the Goldreich torque for omega=[w1,w2,w3]"""

    cdef double chi, epsA, anom_torque, R, c, mx, my, mz, wx, wy, wz

    R = 1e6
    c = 3e10

    epsA = (<double *> params)[0]
    chi = (<double *> params)[1]
    anom_torque = (<double *> params)[4]

    # Check that chi is coming in radians
    if chi > 2 * np.pi:
        print " Check that chi is imported as radians and not degrees"

    wx = y[0]
    wy = y[1]
    wz = y[2]

    mx = sin(chi)
    my = 0.0
    mz = cos(chi)

    omega_squared = sum([pow(w, 2) for w in omega])

    T1 = (2 * R * pow(3 * c, -1) * epsA * omega_squared *
                        np.cross(np.cross(omega, m_vec), m_vec))

    T2 = epsA * np.dot(omega, m_vec) * np.cross(omega, m_vec)

    if anom_torque == 1:
        return T1 + T2
    elif anom_torque == 0:
        return T1

#def Torque_over_Io(wx, wy, wz, epsA, chi, anom_torque):
    #""" Returns the Goldreich torque for omega=[w1,w2,w3]"""

    #R = 1e6
    #c = 3e10

    ## Check that chi is coming in radians
    #if chi > 2 * np.pi:
        #print " Check that chi is imported as radians and not degrees"

    ##cdef double  T1[3], T2[3], omega_squared

    #m_vec = [sin(chi), 0.0, cos(chi)]
    #omega_squared = pow(wx, 2) + pow(wy, 2) + pow(wz, 2)
    #omega = [wx, wy, wz]

    #T1 = (2 * R * pow(3 * c, -1) * epsA * omega_squared *
                        #np.cross(np.cross(omega, m_vec), m_vec))

    #T2 = epsA * np.dot(omega, m_vec) * np.cross(omega, m_vec)

    #if anom_torque == 1:
        #return T1 + T2
    #elif anom_torque == 0:
        #return T1

#cdef int funcs (double t, double y[], double f[], void *params) nogil:
    #""" Function defining the ODEs with the anomalous torque """
    ## Define the variables used in the calculation
    #cdef double wx,wy,wz,w_2,epsI1,epsI3,epsA,chi
    ## Import the three time dependant variables from y[]
    #wx = y[0]
    #wy = y[1]
    #wz = y[2]

    #omega = [wx, wy, wz]
    ## Import the constant variables from params
    #epsA = (<double *> params)[0]
    #chi = (<double *> params)[1]

    #epsI1 = (<double *> params)[2]
    #epsI3 = (<double *> params)[3]

    #anom_torque = (<double *> params)[4]

    #T_o_I = Torque_over_Io(y, params)


    ##  Define the three ODEs in f[] as functions of the above variables

    #f[0] = (T_o_I[0] * pow(1+epsI1,-1)
            #- wy * wz * epsI3 * pow(1 + epsI1, -1))

    #f[1] = (T_o_I[1]
            #- wx * wz * (epsI1 - epsI3))

    #f[0] = (T_o_I[2] * pow(1+epsI3,-1)
            #+ wx * wy * epsI1 * pow(1 + epsI3, -1))

    #return GSL_SUCCESS

cdef int funcs (double t, double y[], double f[], void *params) nogil:
    """ Function defining the ODEs with the anomalous torque """
    # Define the variables used in the calculation
    cdef double wx,wy,wz,w_2,epsI1,epsI3,epsA,chi
    # Import the three time dependant variables from y[]
    wx = y[0]
    wy = y[1]
    wz = y[2]

    # Import the constant variables from params
    epsA = (<double *> params)[0]
    chi = (<double *> params)[1]

    epsI1 = (<double *> params)[2]
    epsI3 = (<double *> params)[3]

    anom_torque = (<double *> params)[4]


    T_o_I = Torque_over_Io(wx, wy, wz, epsA, chi, anom_torque)


    #  Define the three ODEs in f[] as functions of the above variables

    f[0] = (T_o_I[0] * pow(1+epsI1,-1)
            - wy * wz * epsI3 * pow(1 + epsI1, -1))

    f[1] = (T_o_I[1]
            - wx * wz * (epsI1 - epsI3))

    f[0] = (T_o_I[2] * pow(1+epsI3,-1)
            + wx * wy * epsI1 * pow(1 + epsI3, -1))

    return GSL_SUCCESS

cdef int no_anom_torque (double t, double y[], double f[], void *params) nogil:
    """ Function defining the ODEs with the anomalous torque """
    # Define the variables used in the calculation
    cdef double wx,wy,wz,w_2,epsI1,epsI3,epsA,chi,Lambda

    # Import the three time dependant variables from y[]
    wx=y[0]
    wy=y[1]
    wz=y[2]
    w_2 = pow(wx,2)+pow(wy,2)+pow(wz,2)

    # Import the constant variables from params
    Lambda = (<double *> params)[0]    #= 2R/3C */
    epsA = (<double *> params)[1]
    chi = (<double *> params)[2]
    epsI1 = (<double *> params)[3]
    epsI3 = (<double *> params)[4]

    # Calculate the angular parts of the three equations to avoid repeated calculation
    cdef double Cx,Sx
    Cx = cos(chi)
    Sx = sin(chi)

    #  Define the three ODEs in f[] as functions of the above variables
    f[0] = epsA * (Lambda * w_2 * Cx * (wz * Sx - wx * Cx)) * pow(1 + epsI1, -1) - wy * wz * epsI3 * pow(1 +epsI1, -1)

    f[1] = epsA * (- Lambda * w_2 * wy) - wx * wz * (epsI1 - epsI3)

    f[2] = epsA * pow(1 + epsI3, -1) * (Lambda * w_2 * Sx * (wx * Cx - wz * Sx)) + wx * wy * epsI1 * pow(1 + epsI3, -1)
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
    cdef double t , y[3] ,h , eta_relative
    eta_relative = eta*pow(omega0,2)
    h = 1e-15   # Initial step size
    t = 0.0
    y[0] = omega0*sin(a_int)
    y[1] = 0.0
    y[2] = omega0*cos(a_int)
    print y[0],y[1],y[2]

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
            status = gsl_odeiv2_driver_apply (d, &t, ti, y)

            if (status != GSL_SUCCESS):
                print("error, return value=%d\n" % status)
                break

            time.append(t)
            w1.append(y[0])
            w2.append(y[1])
            w3.append(y[2])

        gsl_odeiv2_driver_free(d)

    else:

        while (pow(y[0],2)+pow(y[1],2)+pow(y[2],2) > eta_relative and t < t1 ):
            status = gsl_odeiv2_evolve_apply (e, c, s, &sys, &t, t1, &h, y)

            if (status != GSL_SUCCESS):
                break

            time.append(t)
            w1.append(y[0])
            w2.append(y[1])
            w3.append(y[2])

        gsl_odeiv2_evolve_free (e)
        gsl_odeiv2_control_free (c)
        gsl_odeiv2_step_free (s)



    f = h5py.File(file_name,'w')
    f.create_dataset("time",data=time)
    f.create_dataset("w1",data=w1)
    f.create_dataset("w2",data=w2)
    f.create_dataset("w3",data=w3)
    return file_name

