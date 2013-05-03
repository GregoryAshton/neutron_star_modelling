#!/usr/bin/python 

import File_Functions
import os 

def Run(Input_Dictionary):
	""" Create a generic write file, compile and run in C. 

	Keyword arguments must be passed as a dictionary
	chi:float [degrees] -- angle between magnetic dipole and z axis
	epsI:float []-- elastic deformation 
	epsA :float []-- magnetic deformation 
	omega0 :float [Hz]-- initial spin period 
	either 
		eta :float []-- Threshold for which to stop the simulation. Simulation is initiated using omega**2.0 < eta*omega_0**2.0 
	or
		t1:float [s] -- Time for which to run simulation for		

	Optional arguments
	err : float[] -- Error argument to pass to the GCC compiler 

	For help with the GCC compiler see documentation at http://www.gnu.org/software/gsl/manual/html_node/Ordinary-Differential-Equations.html 
	The generic script is written by Write_File in File_Functions
	"""

	#  Required paramaters

	# Check if the anomalous torque should be used or not and initiate the file_name
	file_name_list = []
	if Input_Dictionary.get('no_anom'):
		file_name_list .append("no_anom_")
		no_anom=True
	else :
		no_anom=False
		
	# At the moment these erros do not raise correctly. 

	try :
		chi = Input_Dictionary['chi']
		file_name_list.append("chi_"+str(Input_Dictionary['chi']))
	except KeyError:
		print " ERROR: You need to specify chi in the input dictionary"
		return

	try :
		epsI= str(Input_Dictionary['epsI'])
		file_name_list.append("_epsI_"+str(Input_Dictionary['epsI']))
	except KeyError:
		print " ERROR: You need to specify epsI in the input dictionary"
		return

	try :
		epsA = str(Input_Dictionary['epsA'])
		file_name_list.append("_epsA_"+str(Input_Dictionary['epsA']))
	except KeyError:
		print " ERROR: You need to specify epsA in the input dictionary"
		return

	try :
		omega0 = str(Input_Dictionary['omega0'])
		file_name_list.append("_omega0_"+str(Input_Dictionary['omega0']))
	except KeyError:
		print " ERROR: You need to specify omega0 in the input dictionary"
		return

	# Take either t1 or eta but not both
	try :
		eta = str(Input_Dictionary['eta'])
		file_name_list.append("_eta_"+str(Input_Dictionary['eta']))
		eta_relative = str(float(eta)*pow(float(omega0),2)) 
	except KeyError:
		eta = None
	try :
		t1 = str(Input_Dictionary['t1'])
		file_name_list.append("_t1_"+str(Input_Dictionary['t1']))
	except KeyError:
		t1 = None

	#  Additional Arguments
	try :
		err = str(Input_Dictionary['err'])
	except KeyError:
		# Use a default value 
		err = 1e-12

	#  Write the generic write_file
	if t1==eta : 
		print " ERROR: You have not specified either eta or t1"
		return 
	if t1 != None:
		# Ugly fix for the moment 
		if 'tri' in Input_Dictionary:
			epsI3 = epsI
			epsI1 = Input_Dictionary['epsI1']
			File_Functions.Write_File_t1_TRIAXIAL(chi,epsI1,epsI3,epsA,omega0,t1)
			file_name_list.reverse()
			file_name_list.append("TRIAXIAL_")
			file_name_list.reverse()
			print "something"
		else :
			File_Functions.Write_File_t1(chi,epsI,epsA,omega0,t1,err,no_anom)
	if eta != None:
		File_Functions.Write_File_eta(chi,epsI,epsA,omega0,eta_relative,err,no_anom)
	
	# Create file name 
	file_name_list.append(".txt")
	file_name = "".join(file_name_list)

	os.system("gcc -Wall -I/usr/local/include -c generic_script.c")
	os.system("gcc -static generic_script.o -lgsl -lgslcblas -lm")
	os.system("./a.out >  %s" % (file_name) )
	os.system("rm generic_script.*  a.out")
	if "verbose" in Input_Dictionary:
		print " Run is complete for this data, the output is saved in the file "+file_name 
	return file_name

def Run_Cython(Input_Dictionary):
	"""

	"""

	from scipy.integrate import ode
	from pylab import cos,sin,pi


	#  Required paramaters

	# Check if the anomalous torque should be used or not and initiate the file_name
	file_name_list = []
	if Input_Dictionary.get('no_anom'):
		file_name_list .append("no_anom")
		no_anom=True
	else :
		no_anom=False
		
	# At the moment these erros do not raise correctly. 

	try :
		chi_degrees = Input_Dictionary['chi']
		file_name_list.append("chi_"+str(Input_Dictionary['chi']))
	except KeyError:
		print " ERROR: You need to specify chi in the input dictionary"
		return

	try :
		epsI= str(Input_Dictionary['epsI'])
		file_name_list.append("_epsI_"+str(Input_Dictionary['epsI']))
	except KeyError:
		print " ERROR: You need to specify epsI in the input dictionary"
		return

	try :
		epsA = str(Input_Dictionary['epsA'])
		file_name_list.append("_epsA_"+str(Input_Dictionary['epsA']))
	except KeyError:
		print " ERROR: You need to specify epsA in the input dictionary"
		return

	try :
		omega0 = str(Input_Dictionary['omega0'])
		file_name_list.append("_omega0_"+str(Input_Dictionary['omega0']))
	except KeyError:
		print " ERROR: You need to specify omega0 in the input dictionary"
		return

	# Take either t1 or eta but not both 
	# Note: For now we will use only t1 
	try :
		eta = str(Input_Dictionary['eta'])
		file_name_list.append("_eta_"+str(Input_Dictionary['eta']))
		eta_relative = str(float(eta)*pow(float(omega0),2)) 
	except KeyError:
		eta = None
	try :
		t1 = str(Input_Dictionary['t1'])
		file_name_list.append("_t1_"+str(Input_Dictionary['t1']))
	except KeyError:
		t1 = None

#	#  Additional Arguments
	try :
		err = str(Input_Dictionary['err'])
	except KeyError:
		# Use a default value 
		err = 1e-12

	#------------------------------------

	import scipy.weave as weave

	def run(chi_degrees,epsI,epsA,omega0,t1,err):

		assert(type(chi_degrees) == float)
		assert(type(epsI) == float)
		assert(type(epsA) == float)
		assert(type(omega0) == float)
		assert(type(t1) == float)
		assert(type(err) == float)
		
		code = '''
/********************************************************************/
/*                                                                  */
/* Evolution of omega , a and phi for ODEs from Goldreich           */
/*                                                                  */
/********************************************************************/


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
	double chi=M_PI*chi_degrees/180; 		// Radians	
	double R = 1.0e6;               	// cm
	double c_speed = 3e10 ;         	// cm/s
	double lambda = 2*R / (3*c_speed) ;     // s

	double eps_I ;           
	double eps_A ;

	/* Initial value of omega, a_int is the initial angle against the z-axis*/    
	double omega_0; 
	double a_int=M_PI*50.0/180; 
	double y[3] = { omega_0*sin(a_int), 0.0 , omega_0*cos(a_int) };

	// Start and stop times for integration:
	double t = 0.0, t1 ;

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
	gsl_odeiv2_control * c = gsl_odeiv2_control_y_new (err, err);

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
		return weave.inline(code,['chi_degrees','epsI','epsA','omega0','t1','err'], compiler = 'gcc',headers = ['<stdio.h>','<gsl/gsl_errno.h>','<gsl/gsl_matrix.h>','<gsl/gsl_odeiv2.h>','<math.h>'])

	run(float(chi_degrees),float(epsI),float(epsA),float(omega0),float(t1),float(err))
	





