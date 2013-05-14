#!/usr/bin/python 

import os 
import nsmod_cython
import nsmod_two_component_model 
import pynotify

# Obsolete code to be deleted by 9/June/2012

#def Run(Input_Dictionary):
#	""" Create a generic write file, compile and run in C. 

#	Keyword arguments must be passed as a dictionary
#	chi:float [degrees] -- angle between magnetic dipole and z axis
#	epsI:float []-- elastic deformation 
#	epsA :float []-- magnetic deformation 
#	omega0 :float [Hz]-- initial spin period 
#	either 
#		eta :float []-- Threshold for which to stop the simulation. Simulation is initiated using omega**2.0 < eta*omega_0**2.0 
#	or
#		t1:float [s] -- Time for which to run simulation for		

#	Optional arguments
#	err : float[] -- Error argument to pass to the GCC compiler 

#	For help with the GCC compiler see documentation at http://www.gnu.org/software/gsl/manual/html_node/Ordinary-Differential-Equations.html 
#	The generic script is written by Write_File in File_Functions
#	"""

#	#  Required paramaters

#	# Check if the anomalous torque should be used or not and initiate the file_name
#	file_name_list = []
#	if Input_Dictionary.get('no_anom'):
#		file_name_list .append("no_anom_")
#		no_anom=True
#	else :
#		no_anom=False
#		
#	# At the moment these erros do not raise correctly. 

#	try :
#		chi = Input_Dictionary['chi']
#		file_name_list.append("chi_"+str(Input_Dictionary['chi']))
#	except KeyError:
#		print " ERROR: You need to specify chi in the input dictionary"
#		return

#	try :
#		epsI= str(Input_Dictionary['epsI'])
#		file_name_list.append("_epsI_"+str(Input_Dictionary['epsI']))
#	except KeyError:
#		print " ERROR: You need to specify epsI in the input dictionary"
#		return

#	try :
#		epsA = str(Input_Dictionary['epsA'])
#		file_name_list.append("_epsA_"+str(Input_Dictionary['epsA']))
#	except KeyError:
#		print " ERROR: You need to specify epsA in the input dictionary"
#		return

#	try :
#		omega0 = str(Input_Dictionary['omega0'])
#		file_name_list.append("_omega0_"+str(Input_Dictionary['omega0']))
#	except KeyError:
#		print " ERROR: You need to specify omega0 in the input dictionary"
#		return

#	# Take either t1 or eta but not both
#	try :
#		eta = str(Input_Dictionary['eta'])
#		file_name_list.append("_eta_"+str(Input_Dictionary['eta']))
#		eta_relative = str(float(eta)*pow(float(omega0),2)) 
#	except KeyError:
#		eta = None
#	try :
#		t1 = str(Input_Dictionary['t1'])
#		file_name_list.append("_t1_"+str(Input_Dictionary['t1']))
#	except KeyError:
#		t1 = None

#	#  Additional Arguments
#	try :
#		err = str(Input_Dictionary['err'])
#	except KeyError:
#		# Use a default value 
#		err = 1e-12

#	#  Write the generic write_file
#	if t1==eta : 
#		print " ERROR: You have not specified either eta or t1"
#		return 
#	if t1 != None:
#		# Ugly fix for the moment 
#		if 'tri' in Input_Dictionary:
#			epsI3 = epsI
#			epsI1 = Input_Dictionary['epsI1']
#			File_Functions.Write_File_t1_TRIAXIAL(chi,epsI1,epsI3,epsA,omega0,t1)
#			file_name_list.reverse()
#			file_name_list.append("TRIAXIAL_")
#			file_name_list.reverse()
#			print "something"
#		else :
#			File_Functions.Write_File_t1(chi,epsI,epsA,omega0,t1,err,no_anom)
#	if eta != None:
#		File_Functions.Write_File_eta(chi,epsI,epsA,omega0,eta_relative,err,no_anom)
#	
#	# Create file name 
#	file_name_list.append(".txt")
#	file_name = "".join(file_name_list)

#	os.system("gcc -Wall -I/usr/local/include -c generic_script.c")
#	os.system("gcc -static generic_script.o -lgsl -lgslcblas -lm")
#	os.system("./a.out >  %s" % (file_name) )
#	os.system("rm generic_script.*  a.out")
#	if "verbose" in Input_Dictionary:
#		print " Run is complete for this data, the output is saved in the file "+file_name 
#	return file_name

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

	# We ask the user to turn of the anom torque, it is by default, on, in the program. 
	if Input_Dictionary.get('no_anom'):
		file_name_list .append("no_anom_")
		anom_torque=False
	else :
		anom_torque=True
		
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

	# Two params can be used to specify the end of the run, eta and t1. The run will continue until the first is satisfied.
	try :
		eta = str(Input_Dictionary['eta'])
		file_name_list.append("_eta_"+str(Input_Dictionary['eta']))
		eta_relative = str(float(eta)*pow(float(omega0),2)) 
	except KeyError:
		eta = 0.0

	try :
		t1 = str(Input_Dictionary['t1'])
		file_name_list.append("_t1_"+str(Input_Dictionary['t1']))
	except KeyError:
		print "ERROR: t1 not specified using a default "
		t1 = 1e15

#	#  Additional Arguments
	try :
		error = str(Input_Dictionary['error'])
	except KeyError:
		# Use a default value 
		error = 1e-12

	try :
		n = int(Input_Dictionary['n'])
		#if "_t1_" not in file_name_list: 
			#print "Warning: Specifying n without t1 means you may get a very course save "
	except KeyError:
		# Don't save at n discrete intervals...
		n=None


	# Create file name 
	file_name_list.append(".hdf5")
	file_name = "".join(file_name_list)

	# Tempory hack to remove old files seems to cause issues
	if file_name in os.listdir("."): os.remove(file_name)
	nsmod_cython.main(chi_degrees=float(chi_degrees),file_name=file_name, n = n,
						 epsA = float(epsA) , epsI=float(epsI) , 
						omega0=float(omega0) , t1 = float(t1) ,
						anom_torque = anom_torque , error=float(error))
	
	pynotify.init("Basic")

	n = pynotify.Notification("Run 1CM complete",
	  "output saved as {}".format(file_name)
	)

	n.show()
	return file_name
	
def Run_Cython_Two_Component(Input_Dictionary):
	""" Create a generic write file, compile and run in C for the two component model. 

	Keyword arguments must be passed as a dictionary
	chi:float [degrees] -- angle between magnetic dipole and z axis
	epsI:float []-- elastic deformation 
	epsA :float []-- magnetic deformation 
	omega0 :float [Hz]-- initial spin period 
	t1:float [s] -- Time for which to run simulation for	
	eta: float 	
	K:float -- Coupling constant between the two components
	Ishell:MOI of the biaxial crust
	Icore:MOI of the spherical core	

	Optional arguments
	err : float[] -- Error argument to pass to the GCC compiler 

	For help with the GCC compiler see documentation at http://www.gnu.org/software/gsl/manual/html_node/Ordinary-Differential-Equations.html 
	"""

	import nsmod_cython

	#  Required paramaters

	# Check if the anomalous torque should be used or not and initiate the file_name
	file_name_list = []

	# We ask the user to turn of the anom torque, it is by default, on, in the program. 
	if Input_Dictionary.get('no_anom'):
		file_name_list .append("no_anom_")
		anom_torque=False
	else :
		anom_torque=True
		
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


	try :
		t1 = str(Input_Dictionary['t1'])
		file_name_list.append("_t1_"+str(Input_Dictionary['t1']))
	except KeyError:
		t1 = "1e15"
		print "ERROR: You have not yet specified t1, using default value t1 = {}".format(t1)
		

	try :
		eta = str(Input_Dictionary['eta'])
		file_name_list.append("_eta_"+str(Input_Dictionary['eta']))
		eta_relative = str(float(eta)*pow(float(omega0),2)) 
	except KeyError:
		eta = "0.0"
		print "ERROR: You have not yet specified eta, using default values  eta ={}".format(eta)
		

	
	try : 
		K = str(Input_Dictionary['K'])
		file_name_list.append("_K_"+str(Input_Dictionary['K']))	
	except KeyError:
		print "ERROR: You must specify K"
		return

	try : 
		Ishell = str(Input_Dictionary['Ishell'])
		file_name_list.append("_Ishell_"+str(Input_Dictionary['Ishell']))	
	except KeyError:
		Ishell = 1e45
		print "ERROR: Ishell not specified using default Ishell={}".format(Ishell)



	try : 
		Icore = str(Input_Dictionary['Icore'])
		file_name_list.append("_Icore_"+str(Input_Dictionary['Icore']))	
	except KeyError:
		Icore = 1e45
		print "ERROR: Icore not specified using default Icore={}".format(Icore)



#	#  Additional Arguments
	try :
		error = str(Input_Dictionary['error'])
	except KeyError:
		# Use a default value 
		error = 1e-5

	try :
		n = int(Input_Dictionary['n'])
		#if "_t1_" not in file_name_list: 
			#print "Warning: Specifying n without t1 means you may get a very course save "
	except KeyError:
		# Don't save at n discrete intervals...
		n=None

	# Create file name 
	file_name_list.append(".hdf5")
	file_name = "".join(file_name_list)

	# Tempory hack to remove old files seems to cause issues
	if file_name in os.listdir("."): os.remove(file_name)
	nsmod_two_component_model.main(chi_degrees=float(chi_degrees) , file_name = file_name,
						 epsA = float(epsA) , epsI=float(epsI) , n=n,
						omega0=float(omega0) , t1 = float(t1) ,
						anom_torque = anom_torque , error=float(error),Ishell=float(Ishell),Icore=float(Icore),K=float(K))
	pynotify.init("Basic")

	n = pynotify.Notification("Run 2CM complete",
	  "output saved as {}".format(file_name)
	)

	n.show()
	return file_name
	





