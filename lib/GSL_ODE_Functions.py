#!/usr/bin/python 

import sys, os

text_no_anomalous_torque= r'''
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

	/* Some computed time scales */
	double ts = pow(lambda*eps_A*omega_0*omega_0,-1);
	double ta = pow(eps_A*omega_0,-1);
	double tp = pow(eps_I*omega_0,-1);

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
	gsl_odeiv2_control * c = gsl_odeiv2_control_y_new (1e-12, 1.0e-12);

	// Created evolution function:
	gsl_odeiv2_evolve * e  = gsl_odeiv2_evolve_alloc (3);

	// Tell GSL the dimension of the system and the names of the functions containing the ODEs:
	gsl_odeiv2_system sys = {func, jac, 3, &params};

	// Trial step size for first step
	double h = 1e-10;

	/*************************************************************/

	/* Print some computed initial data to file */
	printf("chi=%.5e eps_A=%.5e eps_I=%.5e  t_end=%.2e ts=%.2e ta=%.2e tp=%.2e\n", chi , eps_A, eps_I, t1 ,ts,ta,tp);

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
 

text_with_anomalous_torque= r'''
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

	/* Some computed time scales */
	double ts = pow(lambda*eps_A*omega_0*omega_0,-1);
	double ta = pow(eps_A*omega_0,-1);
	double tp = pow(eps_I*omega_0,-1);

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

	/* Print some computed initial data to file */
	printf("chi=%.5e eps_A=%.5e eps_I=%.5e  t_end=%.2e ts=%.2e ta=%.2e tp=%.2e\n", chi , eps_A, eps_I, t1 ,ts,ta,tp);

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

def Write_Generic_Script(chi,eps_I,eps_A,omega_0,t1,err=1.0e-12,*args):

	# We allow two different types of code, one without anomalous torque 
	if "no_anom" in args:
		text = text_no_anomalous_torque
	else : 
		text = text_with_anomalous_torque

	# Replace the keys in the required text 
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



