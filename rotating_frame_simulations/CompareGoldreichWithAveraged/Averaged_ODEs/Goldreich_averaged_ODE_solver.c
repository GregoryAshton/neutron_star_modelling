
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
	double w=y[0]; 
	double a=y[1];

	/* Import the constant variables from params*/
	double *P = (double *) params;
	double lambda = P[0];    /* = 2R/3C */
	double eps_I= P[1];
	double eps_A = P[2];
	double chi = P[3];      

	/* Calculate the angular parts of the three equations to avoid repeating the same calculation */
	/*double Cx = cos(chi) ; */
	double Sx = sin(chi) ;

	double Sa = sin(a) ;
	double Ca = cos(a) ;

	/* Define the three ODEs in f[] as functions of the above variables*/
	f[0] = -lambda*eps_A*pow(w,3)*(pow(Sx,2)+pow(Sa,2)*(1-1.5*pow(Sx,2.0)));

	f[1] = -lambda*eps_A*pow(w,2)*Sa*Ca*(1-1.5*pow(Sx,2.0));



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
	double chi=M_PI*75.0/180; 		// Radians	
	double R = 1.0e6;               	// cm
	double c_speed = 3e10 ;         	// cm/s
	double lambda = 2*R / (3*c_speed) ;     // s

	double eps_I=1e-9;           
	double eps_A=5e-11;

	/* Initial value of omega, a_int is the initial angle against the z-axis*/    
	double omega_0 = 1e4 ; 
	double a_int=M_PI*50.0/180; 
	double y[2] = { omega_0, a_int} ;

	/* Some computed time scales */
	double ts = pow(lambda*eps_A*omega_0*omega_0,-1);
	double ta = pow(eps_A*omega_0,-1);
	double tp = pow(eps_I*omega_0,-1);

	// Start and stop times for integration:
	double t = 0.0 ;
	double t1 = 8.0e10;
	/* double t1 = 1.8e11; time for chi = 30 */

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
	gsl_odeiv2_step * s = gsl_odeiv2_step_alloc (T, 2);

	// Specify absolute and relative errors (actual allowed error is a linear combination of these):
	gsl_odeiv2_control * c = gsl_odeiv2_control_y_new (1e-12, 1.0e-12);

	// Created evolution function:
	gsl_odeiv2_evolve * e  = gsl_odeiv2_evolve_alloc (2);

	// Tell GSL the dimension of the system and the names of the functions containing the ODEs:
	gsl_odeiv2_system sys = {func, jac, 2, &params};

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

		printf ("%.16e %.16e  %.16e\n", t, y[0], y[1] );
	}

       gsl_odeiv2_evolve_free (e);
       gsl_odeiv2_control_free (c);
       gsl_odeiv2_step_free (s);
	return 0;
}

