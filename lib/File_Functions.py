#!/usr/bin/python 

from datetime import datetime
import numpy as np
import pylab as py

def Save_Figure(file_name,type_of_plot,format_type=".png"):
	plot_file_name = type_of_plot+"_"+file_name.rstrip(".txt")+format_type
	py.savefig(plot_file_name)
	print "Saving figure as %s" % plot_file_name

def Parameter_Dictionary(user_input):
	# Could we improve this to give units?
	"""Function to produce a dictionary with all the values, this reduces the amount of code required on the other end. Input can either be a string such as the default filenames used in nsmod, or a partial dictionary having at least epsI,epsA and omega0 to produce the other values"""
	if type(user_input)==str:
		# Initiate dictionary 
		p_d={}

		# Remove the file descriptor and the path directory
		f = user_input.rstrip(".txt")
		f = f.split("/")[-1]	

		# Check if the anomalous torque was used or not
		if "no_anom" in f: 
			f = f.lstrip("no_anom_")
			p_d["no_anom"] = True
		else : p_d["no_anom"] = False
	
		# Import the rest of the parameters
		f = f.split("_")
		for i in range(0,len(f),2):
			p_d[f[i]]=f[i+1]


	elif type(user_input)==dict:
		p_d = user_input

	# Standard values for c , R in cgs
	c = 3e10 
	R = 1e6
	I0=1e45
	# Compute a couple of often used variabes	
	omega0 = float(p_d["omega0"])
	epsI = float(p_d["epsI"])
	epsA = float(p_d["epsA"])
	p_d["tauP"] = str(pow(omega0*epsI,-1))
	p_d["tauA"] = str(pow(omega0*epsA,-1))
	p_d["tauS"] = str(pow(omega0**2.0 *epsA,-1)*3*c/(2*R))
	p_d["Bs"] = str(2*np.sqrt(epsA*I0*R*pow(c,2))/pow(R,3))
	
	# Need to import the beta function
	from Physics_Functions import Beta_Function
	p_d["beta30"] = str(Beta_Function(epsI,epsA,30*np.pi/180)*180/np.pi)
	p_d["beta75"] = str(Beta_Function(epsI,epsA,75*np.pi/180)*180/np.pi)
	return p_d

def Import_Data(file_name,max_int=-1,d_int=1):
	if max_int != -1 : max_int = 4 *max_int # Note it is 4 since each line contains 4 columns so the total data set is 4*nmax
	data = np.fromfile(file_name,dtype=float,count=max_int,sep=" ") 
	n = len(data)
	time = data[0:n:4]
	x = data[1:n+1:4]
	y = data[2:n+2:4]
	z = data[3:n+3:4]
	print len(time),len(x)

	return (time,x,y,z)


def Import_Data_Array_Method(file_name,max_int=None,*args):
	t1 = datetime.now() # Used for testing time 
	""" Return the contents of 'file_name' as 4 lists + a list of params"""
	params = open(file_name,"r").readline().split()
	data = np.genfromtxt(file_name, skip_header=1,delimiter=" ")

	time = data[:,0]
	x = data[:,1]
	y = data[:,2]
	z = data[:,3]

	t2 = datetime.now()
	if "timing" in args:
		print "Time for operations = ",t2-t1	

	return (time,x,y,z,params)

def Write_File_t1(chi,eps_I,eps_A,omega_0,t1,err=1.0e-12,no_anom=False):
	"""Write a generic script to be compiles and run in C with paramaters chi,eps_I,eps_A,omega_0,t1.err=1.0e-12, if "no_anom" in *args the generated code will not contain the anomalous torque"""

	if "no_anom":
		text= r'''
/* THIS IS AN AUTOMATED FILE AND SHOULD NOT BE EDITED */

/********************************************************************/
/*                                                                  */
/* Evolution of omega , a and phi for ODEs from Goldreich           */
/*                                                                  */
/********************************************************************/

#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <math.h>

/* Function containing the ODEs*/

int
func (double t, const double y[], double f[] , void *params)
{
	/* Import the three time dependant variables from y[]*/
	double wx=y[0]; 
	double wy=y[1];
	double wz=y[2];

	/* Calculate w*w labelled as w_2 */
	double w_2 = pow(wx,2.0)+pow(wy,2.0)+pow(wz,2.0);

	/* Import the constant variables from params*/
	double *P = (double *) params;
	double lambda = P[0];    /* = 2R/3C */
	double eps_I= P[1];
	double eps_A = P[2];
	double chi = P[3];      

	/* Calculate the angular parts of the three equations to avoid repeating the same calculation */
	double Cx = cos(chi) ;
	double Sx = sin(chi) ;

	/* Define the three ODEs in f[] as functions of the above variables*/
	f[0] = eps_A*(lambda*w_2*Cx*(wz*Sx-wx*Cx)) - wy*wz*eps_I;

	f[1] = eps_A*(-lambda*w_2*wy) + wx*wz*eps_I;

	f[2] = (eps_A/(1+eps_I)) * (lambda*w_2*Sx*(wx*Cx - wz*Sx)) ;

	return GSL_SUCCESS;
}

/* Currently jac is unused by the ODE solver so is left empty*/
int
jac (double t, const double y[], double *dfdy, 
	 double dfdt[], void *params)
{

	return GSL_SUCCESS;
}

int
main (void)
{
	/* Input parameters*/
	double chi=M_PI*KEY_CHI/180; 		// Radians	
	double R = 1.0e6;               	// cm
	double c_speed = 3e10 ;         	// cm/s
	double lambda = 2*R / (3*c_speed) ;     // s

	double eps_I=KEY_EPS_I;           
	double eps_A=KEY_EPS_A;

	/* Initial value of omega, a_int is the initial angle against the z-axis*/    
	double omega_0 = KEY_OMEGA_0 ; 
	double a_int=M_PI*50.0/180; 
	double y[3] = { omega_0*sin(a_int), 0.0 , omega_0*cos(a_int) };

	// Start and stop times for integration:
	double t = 0.0, t1 = KEY_t1;

	double params[4];	
	params[0] = lambda;
	params[1] = eps_I ; 
	params[2] = eps_A;
	params[3] = chi;


	/*************************************************************/
	/* GSL ODE desclarations                                     */
	/*************************************************************/

	// Choose the stepping alogrithm:
	const gsl_odeiv2_step_type * T  = gsl_odeiv2_step_rk8pd;

	 	// Specify number of dependent variables:
	gsl_odeiv2_step * s = gsl_odeiv2_step_alloc (T, 3);

	// Specify absolute and relative errors (actual allowed error is a linear combination of these):
	gsl_odeiv2_control * c = gsl_odeiv2_control_y_new (KEY_err, KEY_err);

	// Created evolution function:
	gsl_odeiv2_evolve * e  = gsl_odeiv2_evolve_alloc (3);

	// Tell GSL the dimension of the system and the names of the functions containing the ODEs:
	gsl_odeiv2_system sys = {func, jac, 3, &params};

	// Trial step size for first step
	double h = 1e-10;

	/*************************************************************/


	while (t < t1)
	{
		int status = gsl_odeiv2_evolve_apply (e, c, s, &sys, &t, t1, &h, y);
		if (status != GSL_SUCCESS)
	 	{
			printf ("error, return value=%d\n", status);
			break;
	 	}

		printf ("%.16e %.16e %.16e %.16e\n", t, y[0], y[1] , y[2]);
	}

	   gsl_odeiv2_evolve_free (e);
	   gsl_odeiv2_control_free (c);
	   gsl_odeiv2_step_free (s);
	return 0;
}

	'''

	else : 

		text= r'''
/* THIS IS AN AUTOMATED FILE AND SHOULD NOT BE EDITED */

/********************************************************************/
/*                                                                  */
/* Evolution of omega , a and phi for ODEs from Goldreich           */
/*                                                                  */
/********************************************************************/

#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <math.h>

/* Function containing the ODEs*/

int
func (double t, const double y[], double f[] , void *params)
{
	/* Import the three time dependant variables from y[]*/
	double wx=y[0]; 
	double wy=y[1];
	double wz=y[2];

	/* Calculate w*w labelled as w_2 */
	double w_2 = pow(wx,2.0)+pow(wy,2.0)+pow(wz,2.0);

	/* Import the constant variables from params*/
	double *P = (double *) params;
	double lambda = P[0];    /* = 2R/3C */
	double eps_I= P[1];
	double eps_A = P[2];
	double chi = P[3];      

	/* Calculate the angular parts of the three equations to avoid repeating the same calculation */
	double Cx = cos(chi) ;
	double Sx = sin(chi) ;

	/* Define the three ODEs in f[] as functions of the above variables*/
	f[0] = eps_A*(lambda*w_2*Cx*(wz*Sx-wx*Cx)+(Sx*wx+Cx*wz)*wy*Cx) - wy*wz*eps_I;

	f[1] = eps_A*(-lambda*w_2*wy+(Sx*wx+Cx*wz)*(wz*Sx-wx*Cx)) + wx*wz*eps_I;

	f[2] = (eps_A/(1+eps_I)) * (lambda*w_2*Sx*(wx*Cx - wz*Sx) -(Sx*wx+Cx*wz)*wy*Sx) ;

	return GSL_SUCCESS;
}

/* Currently jac is unused by the ODE solver so is left empty*/
int
jac (double t, const double y[], double *dfdy, 
	 double dfdt[], void *params)
{

	return GSL_SUCCESS;
}

int
main (void)
{
	/* Input parameters*/
	double chi=M_PI*KEY_CHI/180; 		// Radians	
	double R = 1.0e6;               	// cm
	double c_speed = 3e10 ;         	// cm/s
	double lambda = 2*R / (3*c_speed) ;     // s

	double eps_I=KEY_EPS_I;           
	double eps_A=KEY_EPS_A;

	/* Initial value of omega, a_int is the initial angle against the z-axis*/    
	double omega_0 = KEY_OMEGA_0 ; 
	double a_int=M_PI*50.0/180; 
	double y[3] = { omega_0*sin(a_int), 0.0 , omega_0*cos(a_int) };

	// Start and stop times for integration:
	double t = 0.0, t1 = KEY_t1;

	double params[4];	
	params[0] = lambda;
	params[1] = eps_I ; 
	params[2] = eps_A;
	params[3] = chi;


	/*************************************************************/
	/* GSL ODE desclarations                                     */
	/*************************************************************/

	// Choose the stepping alogrithm:
	const gsl_odeiv2_step_type * T  = gsl_odeiv2_step_rk8pd;

	 	// Specify number of dependent variables:
	gsl_odeiv2_step * s = gsl_odeiv2_step_alloc (T, 3);

	// Specify absolute and relative errors (actual allowed error is a linear combination of these):
	gsl_odeiv2_control * c = gsl_odeiv2_control_y_new (KEY_err, KEY_err);

	// Created evolution function:
	gsl_odeiv2_evolve * e  = gsl_odeiv2_evolve_alloc (3);

	// Tell GSL the dimension of the system and the names of the functions containing the ODEs:
	gsl_odeiv2_system sys = {func, jac, 3, &params};

	// Trial step size for first step
	double h = 1e-10;

	/*************************************************************/

	while (t < t1)
	{
		int status = gsl_odeiv2_evolve_apply (e, c, s, &sys, &t, t1, &h, y);
		if (status != GSL_SUCCESS)
	 	{
			printf ("error, return value=%d\n", status);
			break;
	 	}

		printf ("%.16e %.16e %.16e %.16e\n", t, y[0], y[1] , y[2]);
	}

	   gsl_odeiv2_evolve_free (e);
	   gsl_odeiv2_control_free (c);
	   gsl_odeiv2_step_free (s);
	return 0;
}

'''


	text=text.replace("KEY_CHI",str(chi))
	text=text.replace("KEY_EPS_I",str(eps_I))
	text=text.replace("KEY_EPS_A",str(eps_A))
	text=text.replace("KEY_OMEGA_0",str(omega_0))
	text=text.replace("KEY_t1",str(t1))
	text=text.replace("KEY_err",str(err))

	file_name="generic_script.c"
	write_file=open(file_name,"w+")
	write_file.write(text)
	write_file.close()


def Write_File_eta(chi,eps_I,eps_A,omega_0,eta_relative,err=1.0e-12,no_anom=False):
	"""Write a generic script to be compiles and run in C with paramaters chi,eps_I,eps_A,omega_0,eta,err=1.0e-12, if "no_anom" in the generated code will not contain the anomalous torque"""

	if "no_anom" :
		text= r'''
/* THIS IS AN AUTOMATED FILE AND SHOULD NOT BE EDITED */

/********************************************************************/
/*                                                                  */
/* Evolution of omega , a and phi for ODEs from Goldreich           */
/*                                                                  */
/********************************************************************/

#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <math.h>

/* Function containing the ODEs*/

int
func (double t, const double y[], double f[] , void *params)
{
	/* Import the three time dependant variables from y[]*/
	double wx=y[0]; 
	double wy=y[1];
	double wz=y[2];

	/* Calculate w*w labelled as w_2 */
	double w_2 = pow(wx,2.0)+pow(wy,2.0)+pow(wz,2.0);

	/* Import the constant variables from params*/
	double *P = (double *) params;
	double lambda = P[0];    /* = 2R/3C */
	double eps_I= P[1];
	double eps_A = P[2];
	double chi = P[3];      

	/* Calculate the angular parts of the three equations to avoid repeating the same calculation */
	double Cx = cos(chi) ;
	double Sx = sin(chi) ;

	/* Define the three ODEs in f[] as functions of the above variables*/
	f[0] = eps_A*(lambda*w_2*Cx*(wz*Sx-wx*Cx)) - wy*wz*eps_I;

	f[1] = eps_A*(-lambda*w_2*wy) + wx*wz*eps_I;

	f[2] = (eps_A/(1+eps_I)) * (lambda*w_2*Sx*(wx*Cx - wz*Sx)) ;

	return GSL_SUCCESS;
}

/* Currently jac is unused by the ODE solver so is left empty*/
int
jac (double t, const double y[], double *dfdy, 
	 double dfdt[], void *params)
{

	return GSL_SUCCESS;
}

int
main (void)
{
	/* Input parameters*/
	double chi=M_PI*KEY_CHI/180; 		// Radians	
	double R = 1.0e6;               	// cm
	double c_speed = 3e10 ;         	// cm/s
	double lambda = 2*R / (3*c_speed) ;     // s

	double eps_I=KEY_EPS_I;           
	double eps_A=KEY_EPS_A;

	/* Initial value of omega, a_int is the initial angle against the z-axis*/    
	double omega_0 = KEY_OMEGA_0 ; 
	double a_int=M_PI*50.0/180; 
	double y[3] = { omega_0*sin(a_int), 0.0 , omega_0*cos(a_int) };

	// Start and stop times for integration in the automated version we choose a maximum value of t1:
	double t = 0.0, t1 = 1e15; 

	double params[4];	
	params[0] = lambda;
	params[1] = eps_I ; 
	params[2] = eps_A;
	params[3] = chi;


	/*************************************************************/
	/* GSL ODE desclarations                                     */
	/*************************************************************/

	// Choose the stepping alogrithm:
	const gsl_odeiv2_step_type * T  = gsl_odeiv2_step_rk8pd;

	 	// Specify number of dependent variables:
	gsl_odeiv2_step * s = gsl_odeiv2_step_alloc (T, 3);

	// Specify absolute and relative errors (actual allowed error is a linear combination of these):
	gsl_odeiv2_control * c = gsl_odeiv2_control_y_new (KEY_err, KEY_err);

	// Created evolution function:
	gsl_odeiv2_evolve * e  = gsl_odeiv2_evolve_alloc (3);

	// Tell GSL the dimension of the system and the names of the functions containing the ODEs:
	gsl_odeiv2_system sys = {func, jac, 3, &params};

	// Trial step size for first step
	double h = 1e-10;

	/*************************************************************/
	
	/* square value of spin magnitude to cut of simulation */
	double eta = KEY_ETA_REL ;

	while (pow(y[0],2)+pow(y[1],2)+pow(y[2],2) > eta && t < t1)
	{
		int status = gsl_odeiv2_evolve_apply (e, c, s, &sys, &t, t1, &h, y);
		if (status != GSL_SUCCESS)
	 	{
			printf ("error, return value=%d\n", status);
			break;
	 	}

		printf ("%.16e %.16e %.16e %.16e\n", t, y[0], y[1] , y[2]);
	}

	   gsl_odeiv2_evolve_free (e);
	   gsl_odeiv2_control_free (c);
	   gsl_odeiv2_step_free (s);
	return 0;
}

	'''

	else : 

		text= r'''
/* THIS IS AN AUTOMATED FILE AND SHOULD NOT BE EDITED */

/********************************************************************/
/*                                                                  */
/* Evolution of omega , a and phi for ODEs from Goldreich           */
/*                                                                  */
/********************************************************************/

#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <math.h>

/* Function containing the ODEs*/

int
func (double t, const double y[], double f[] , void *params)
{
	/* Import the three time dependant variables from y[]*/
	double wx=y[0]; 
	double wy=y[1];
	double wz=y[2];

	/* Calculate w*w labelled as w_2 */
	double w_2 = pow(wx,2.0)+pow(wy,2.0)+pow(wz,2.0);

	/* Import the constant variables from params*/
	double *P = (double *) params;
	double lambda = P[0];    /* = 2R/3C */
	double eps_I= P[1];
	double eps_A = P[2];
	double chi = P[3];      

	/* Calculate the angular parts of the three equations to avoid repeating the same calculation */
	double Cx = cos(chi) ;
	double Sx = sin(chi) ;

	/* Define the three ODEs in f[] as functions of the above variables*/
	f[0] = eps_A*(lambda*w_2*Cx*(wz*Sx-wx*Cx)+(Sx*wx+Cx*wz)*wy*Cx) - wy*wz*eps_I;

	f[1] = eps_A*(-lambda*w_2*wy+(Sx*wx+Cx*wz)*(wz*Sx-wx*Cx)) + wx*wz*eps_I;

	f[2] = (eps_A/(1+eps_I)) * (lambda*w_2*Sx*(wx*Cx - wz*Sx) -(Sx*wx+Cx*wz)*wy*Sx) ;

	return GSL_SUCCESS;
}

/* Currently jac is unused by the ODE solver so is left empty*/
int
jac (double t, const double y[], double *dfdy, 
	 double dfdt[], void *params)
{

	return GSL_SUCCESS;
}

int
main (void)
{
	/* Input parameters*/
	double chi=M_PI*KEY_CHI/180; 		// Radians	
	double R = 1.0e6;               	// cm
	double c_speed = 3e10 ;         	// cm/s
	double lambda = 2*R / (3*c_speed) ;     // s

	double eps_I=KEY_EPS_I;           
	double eps_A=KEY_EPS_A;

	/* Initial value of omega, a_int is the initial angle against the z-axis*/    
	double omega_0 = KEY_OMEGA_0 ; 
	double a_int=M_PI*50.0/180; 
	double y[3] = { omega_0*sin(a_int), 0.0 , omega_0*cos(a_int) };

	// Start and stop times for integration in the automated version we choose a maximum value of t1:
	double t = 0.0, t1 = 1e15;

	double params[4];	
	params[0] = lambda;
	params[1] = eps_I ; 
	params[2] = eps_A;
	params[3] = chi;


	/*************************************************************/
	/* GSL ODE desclarations                                     */
	/*************************************************************/

	// Choose the stepping alogrithm:
	const gsl_odeiv2_step_type * T  = gsl_odeiv2_step_rk8pd;

	 	// Specify number of dependent variables:
	gsl_odeiv2_step * s = gsl_odeiv2_step_alloc (T, 3);

	// Specify absolute and relative errors (actual allowed error is a linear combination of these):
	gsl_odeiv2_control * c = gsl_odeiv2_control_y_new (KEY_err, KEY_err);

	// Created evolution function:
	gsl_odeiv2_evolve * e  = gsl_odeiv2_evolve_alloc (3);

	// Tell GSL the dimension of the system and the names of the functions containing the ODEs:
	gsl_odeiv2_system sys = {func, jac, 3, &params};

	// Trial step size for first step
	double h = 1e-10;

	/*************************************************************/

	/* square value of spin magnitude to cut of simulation */
	double eta = KEY_ETA_REL ;

	while (pow(y[0],2)+pow(y[1],2)+pow(y[2],2) > eta && t < t1) 
	{
		int status = gsl_odeiv2_evolve_apply (e, c, s, &sys, &t, t1, &h, y);
		if (status != GSL_SUCCESS)
	 	{
			printf ("error, return value=%d\n", status);
			break;
	 	}

		printf ("%.16e %.16e %.16e %.16e\n", t, y[0], y[1] , y[2]);
	}

	   gsl_odeiv2_evolve_free (e);
	   gsl_odeiv2_control_free (c);
	   gsl_odeiv2_step_free (s);
	return 0;
}

'''


	text=text.replace("KEY_CHI",str(chi))
	text=text.replace("KEY_EPS_I",str(eps_I))
	text=text.replace("KEY_EPS_A",str(eps_A))
	text=text.replace("KEY_OMEGA_0",str(omega_0))
	text=text.replace("KEY_ETA_REL",str(eta_relative))
	text=text.replace("KEY_err",str(err))

	file_name="generic_script.c"
	write_file=open(file_name,"w+")
	write_file.write(text)
	write_file.close()




