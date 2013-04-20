#!/usr/bin/python 

import optparse, sys, os
import pylab as py

# Import internal modules
import lib.File_Functions as File_Functions

def Magnetic_Field_to_Epsilon_A(Bs):
	R = 1e6 #cm
	c=3e10 #cm/s
	I0=1e45
	m=0.5*Bs*pow(R,3)
	epsA=pow(m,2)/(I0*R*pow(c,2))
	print "Bs="+str(Bs)+" Epsilon_A="+str(epsA)

def Epsilon_A_to_Magnetic_Field(Epsilon_A,Option_Dictionary={}):
	R = 1e6 #cm
	c=3e10 #cm/s
	I0=1e45
	m = py.sqrt(Epsilon_A*I0*R*pow(c,2))
	Bs = 2*m/pow(R,3)
	if Option_Dictionary['verbose']==True:
		print "Bs="+str(Bs)+" Epsilon_A="+str(Epsilon_A)
	else :
		return  Bs

def Run(Input_Dictionary):
	""" Create a generic write file, compile and run in C. 

	Keyword arguments must be passed as a dictionary
	chi=float [degrees] -- angle between magnetic dipole and z axis
	epsI=float []-- elastic deformation 
	epsA =float []-- magnetic deformation 
	omega_0 =float [Hz]-- initial spin period 
	either 
		eta =float []-- Threshold for which to stop the simulation. Simulation is initiated using omega**2.0 < eta*omega_0**2.0 
	or
		t1=float [s] -- Time for which to run simulation for		

	Optional arguments
	err = float[] -- Error argument to pass to the GCC compiler 

	For help with the GCC compiler see documentation at http://www.gnu.org/software/gsl/manual/html_node/Ordinary-Differential-Equations.html 
	The generic script is written by Write_File in File_Functions
	"""

	#  Required paramaters

	# Check if the anomalous torque should be used or not and initiate the file_name
	file_name_list = []
	if Input_Dictionary.get('no_anom'):
		file_name_list .append("no_anom")
		args="no_anom"
		
	try :
		chi = Input_Dictionary['chi']
		file_name_list.append("_chi_"+str(Input_Dictionary['chi']))
	except KeyError:
		print " ERROR: You need to specify chi in the input dictionary"

	try :
		epsI= str(Input_Dictionary['epsI'])
		file_name_list.append("_epsI_"+str(Input_Dictionary['epsI']))
	except KeyError:
		print " ERROR: You need to specify epsI in the input dictionary"

	try :
		epsA = str(Input_Dictionary['epsA'])
		file_name_list.append("_epsA_"+str(Input_Dictionary['epsA']))
	except KeyError:
		print " ERROR: You need to specify epsA in the input dictionary"

	try :
		omega0 = str(Input_Dictionary['omega0'])
		file_name_list.append("_omega0_"+str(Input_Dictionary['omega0']))
	except KeyError:
		print " ERROR: You need to specify omega0 in the input dictionary"

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

	if t1 and eta : 
		print " ERROR: You have not specified either eta or t1"
	else : 	
		print "ERROR: You have specified both eta and t1, this is incorrect"
		return 
	
	# Optional Arguments
	try :
		err = str(Input_Dictionary['err'])
	except KeyError:
		# Use a default value 
		err = 1e-12


	# Create file name 
	file_name = "".join(file_name_list)
	print file_name
#	if Option_Dictionary.has_key('eta'):
#		if Option_Dictionary['eta']!=False:
#			eta = str(Option_Dictionary['eta'])
#			# Caluclate the relative eta 
#			

#			if Option_Dictionary.has_key('no_anom') and Option_Dictionary['no_anom']==True:
#				print " Running code WITHOUT the anomalous torque"
#				file_name = "no_anom_chi_%s_epsI_%s_epsA_%s_omega0_%s_eta_%s.txt" % (chi,epsI,epsA,omega0,eta) 
#				args="no_anom"
#			else :
#				print " Running code WITH the anomalous torque"
#				file_name = "chi_%s_epsI_%s_epsA_%s_omega0_%s_eta_%s.txt" % (chi,epsI,epsA,omega0,eta) 
#				args = None

#			File_Functions.Write_File_Automatic(chi,epsI,epsA,omega0,eta_relative,err,args)

#	if Option_Dictionary.has_key('t1'):
#		t1 = str(Option_Dictionary['t1'])

#		if Option_Dictionary.has_key('no_anom') and Option_Dictionary['no_anom']==True:
#			print " Running code WITHOUT the anomalous torque"
#			file_name = "no_anom_chi_%s_epsI_%s_epsA_%s_omega0_%s_t1_%s.txt" % (chi,epsI,epsA,omega0,t1) 
#			args="no_anom"
#		else :
#			print " Running code WITH the anomalous torque"
#			file_name = "chi_%s_epsI_%s_epsA_%s_omega0_%s_t1_%s.txt" % (chi,epsI,epsA,omega0,t1) 
#			args = None

#		File_Functions.Write_File(chi,epsI,epsA,omega0,t1,err,args) # Note this is not the automatic writer..consider relabelling

	if 'file_name' not in locals(): print "You have not specified either eta or t1"
	try :
		file_name
	except NameError:
		print "You have not specified either eta or t1" # This needs to be fixed
		
	os.system("gcc -Wall -I/usr/local/include -c generic_script.c")
	os.system("gcc -static generic_script.o -lgsl -lgslcblas -lm")
	os.system("./a.out >  %s" % (file_name) )
	os.system("rm generic_script.*  a.out")
	print " Run is complete for this data, the output is saved in the file "+file_name
	print 
	return file_name

def Print_Parameters(file_name):
	""" Print a list of parameters about the file to the terminal """
	from lib.File_Functions import Parameter_Dictionary
	Parameter_Dictionary = Parameter_Dictionary(file_name)
	for item in Parameter_Dictionary:
		print " %s = %s " % (item,Parameter_Dictionary[item])
	

def Create_Dictionary(opts):
	Option_Dictionary={}
	for item in opts.split(","):
		if "=" in item: Option_Dictionary[item.split("=")[0]] = float(item.split("=")[1])
	else : Option_Dictionary[item] = True
	return Option_Dictionary
	

def main():
	def parse_command_line(argvs):

		parser = optparse.OptionParser()
	
		parser.add_option("--b_2_e", help="Return the epsilon_A value given a surface magnetic field (canonical value of R=1e6 cm is assumed)",metavar="Bs" )

		parser.add_option("--e_2_b", help="Return the Bs value given a value of Epsilon_A (canonical value of R=1e6 cm is assumed)",metavar="Epsilon_A" )

		parser.add_option("-B", "--beta" ,help="Takes as input =sign,epsI,epsA,chi and returns the angle Beta through which the xz coordinates rotate", metavar="FILE")

		parser.add_option("-r", "--run" , help = Run.__doc__)

		parser.add_option("-p","--print_parameters",help = Print_Parameters.__doc__)

		# Additional arguments are passed to opts 
		parser.add_option("-o","--opts")
	
		# Set the verbose/quite options, this will be passed to the option dictionary
		parser.add_option("-v", action="store_true", dest="verbose", default=True)
		parser.add_option("-q", action="store_false", dest="verbose")

		(options,arguments) = parser.parse_args(argvs)
		return options, arguments

	options, arguments = parse_command_line(sys.argv)

	# Create the options dictionary 
	if options.opts : Option_Dictionary = Create_Dictionary(options.opts)
	else : Option_Dictionary = {}

	# Add the verbosity to the Option Dictionary
	Option_Dictionary['verbose'] = options.verbose

	if options.b_2_e : Magnetic_Field_to_Epsilon_A(options.b_2_e,Option_Dictionary)

	if options.e_2_b : Epsilon_A_to_Magnetic_Field(options.e_2_b,Option_Dictionary)

	if options.beta : Print_Beta(options.beta)

	# For the Run command the argument should be a dictionary of values. We create this from the argument, need to document
	Input_Dictionary = Create_Dictionary(options.run)
	if options.run : Run(Input_Dictionary)

	if options.print_parameters : Print_Parameters(options.print_parameters)


if __name__ == "__main__":
    main()


