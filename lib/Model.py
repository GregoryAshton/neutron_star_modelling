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
		file_name_list .append("no_anom")
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
