#!/usr/bin/python 

import lib.File_Functions as File_Functions

import optparse, sys, os
import pylab as py

# Define the functions

def Beta_func(epsI,epsA,chi):
	if epsI>=0:
		sign=-1.0
	else: sign=1.0

	b=py.sqrt(epsA*epsA+epsI*epsI-2*epsA*epsI*py.cos(2*chi))
	return py.arctan((epsI+epsA*(1-2*pow(py.cos(chi),2.0))+sign*b)/(2*epsA*py.sin(chi)*py.cos(chi)))

def Print_Beta(foptions):
	"""This functions takes as input the three parameters epsI,epsA and chi and returns Beta the angle between the principle axis of the new MOI tensor which aligns with the original z axis, and the original z axis."""
	foptions = foptions.split(",")
	sign=float(foptions[0])
	epsI=sign*pow(10,float(foptions[1]))
	epsA=pow(10,float(foptions[2]))
	chi=float(foptions[3])*pi/180.0
	print "epsI="+str(epsI)+"  epsA="+str(round_to_n(epsA,2))+"  chi="+str(foptions[3])
	print "beta=" , Beta_func(epsI,epsA,chi)*180.0/pi, " degrees"

def Magnetic_Field_to_Epsilon_A(Bs):
	R = 1e6 #cm
	c=3e10 #cm/s
	I0=1e45
	m=0.5*Bs*pow(R,3)
	epsA=pow(m,2)/(I0*R*pow(c,2))
	print "Bs="+str(Bs)+" Epsilon_A="+str(epsA)

def Epsilon_A_to_Magnetic_Field(Epsilon_A):
	R = 1e6 #cm
	c=3e10 #cm/s
	I0=1e45
	m = py.sqrt(Epsilon_A*I0*R*pow(c,2))
	Bs = 2*m/pow(R,3)
	print "Bs="+str(Bs)+" Epsilon_A="+str(Epsilon_A)


# Top level function to sort the dictionary of parameters write a .c script compile and output to a suitable .txt file
# may be able to use http://docs.python.org/2/tutorial/controlflow.html instead eg **dictionary
def Run(Option_Dictionary):
	""" Takes as input a dictionary of parameters chi,epsI,epsA,omega_0,eta,err=1.0e-12,*args from the input. chi must be unit degrees and not radians"""

	# Minimum number of input parameters
	if Option_Dictionary.has_key('chi') : chi = str(Option_Dictionary['chi'])
	else : print "ERROR you have not specified chi" 

	if Option_Dictionary.has_key('epsI') : epsI = str(Option_Dictionary['epsI'])
	else : print "ERROR you have not specified epsI"

	if Option_Dictionary.has_key('epsA') : epsA = str(Option_Dictionary['epsA'])
	else : print "ERROR you have not specified epsA"

	if Option_Dictionary.has_key('omega0') : omega0 = str(Option_Dictionary['omega0'])
	else : print "ERROR you have not specified omega0"

	if Option_Dictionary.has_key('eta') : eta = str(Option_Dictionary['eta'])
	else : print "ERROR you have not specified eta"

	if Option_Dictionary.has_key('err') : err = str(Option_Dictionary['err'])
	else : 
		#print " Using default error value of 1e-12" 
		err = 1e-12

	# Caluclate the relative eta 
	eta_relative = str(float(eta)*pow(float(omega0),2))

	# Create file name 
	if Option_Dictionary.has_key('no_anom') and Option_Dictionary['no_anom']==True:
		print " Running code WITHOUT the anomalous torque"
		file_name = "no_anom_chi_%s_epsI_%s_epsA_%s_omega0_%s_eta_%s.txt" % (chi,epsI,epsA,omega0,eta) 
		args="no_anom"
	else :
		print " Running code WITH the anomalous torque"
		file_name = "chi_%s_epsI_%s_epsA_%s_omega0_%s_eta_%s.txt" % (chi,epsI,epsA,omega0,eta) 
		args = None

	File_Functions.Write_File_Automatic(chi,epsI,epsA,omega0,eta_relative,err,args)
	os.system("gcc -Wall -I/usr/local/include -c generic_script.c")
	os.system("gcc -static generic_script.o -lgsl -lgslcblas -lm")
	os.system("./a.out >  %s" % (file_name) )
	os.system("rm generic_script.*  a.out")
	print " Run is complete for this data, the output is saved in the file "+file_name
	print 
	return file_name

def Create_Option_Dictionary(opts):
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

		parser.add_option("-r", "--run" , help = 
									"""Create the generic write file that can be compiled
									 and run in C. The input should be of the form of a dictionary including 

									chi,epsI,epsA,omega_0,eta,err=1.0e-12,*args 

									where
									chi = angle between magnetic dipole and z axis
									epsI = elastic deformation 
									epsA magnetic deformation 
									omega_0 = initial spin period 
									eta = Simultion will stop when omega**2.0 < eta*omega_0**2.0 """)

		parser.add_option("-o","--opts")


		# Options to be called with above

		(options,arguments) = parser.parse_args(argvs)
		return options, arguments

	options, arguments = parse_command_line(sys.argv)

	if options.opts : Option_Dictionary = Create_Option_Dictionary(options.opts)
	else : Option_Dictionary = {}

	if options.b_2_e : Magnetic_Field_to_Epsilon_A(options.b_2_e)

	if options.e_2_b : Epsilon_A_to_Magnetic_Field(options.e_2_b)

	if options.beta : Print_Beta(options.beta)

	# For the fun command the argument should be a dictionary of values. No Option_Dictionary exists
	if options.run : Run(options.run)


if __name__ == "__main__":
    main()


