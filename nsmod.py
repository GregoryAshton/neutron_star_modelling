#!/usr/bin/python 

import argparse, sys, os
import pylab as py

# Import external modules
import lib.Model as Model
import lib.NLD_Functions as NLD_Functions
import lib.Plot as Plot
import lib.Useful_Tools as Useful_Tools



#def Magnetic_Field_to_Epsilon_A(Bs):
#	Bs=float(Bs)
#	R = 1e6 #cm
#	c = 3e10 #cm/s
#	I0=1e45
#	m=0.5*Bs*pow(R,3)
#	epsA=pow(m,2)/(I0*R*pow(c,2))
#	print "Bs="+str(Bs)+" Epsilon_A="+str(epsA)

#def Epsilon_A_to_Magnetic_Field(Epsilon_A,Option_Dictionary={}):
#	R = 1e6 #cm
#	c=3e10 #cm/s
#	I0=1e45
#	m = py.sqrt(Epsilon_A*I0*R*pow(c,2))
#	Bs = 2*m/pow(R,3)
#	if Option_Dictionary['verbose']==True:
#		print "Bs="+str(Bs)+" Epsilon_A="+str(Epsilon_A)
#	else :
#		return  Bs



def Create_Dictionary(opts):
	"""

	Takes as input a string and returns a dictionary

	Each key:value should be deliminated by a comma, all values will be passed as strings e.g
	'x:10.0,y:True,z=:10.0,2.0'
	will return 
	{'x':'10.0','y':'True','z':'10.0,2.0'}
	Terms with no value will be treated as 'True'

	"""

	Option_Dictionary={}
	for item in opts.split(","):
		if ":" in item: Option_Dictionary[item.split(":")[0]] = item.split(":")[1]
	else : Option_Dictionary[item] = 'True'

	return Option_Dictionary
	
def parse_command_line(argvs):

	parser = argparse.ArgumentParser() 

#		parser.add_option("-B", "--beta" ,help="Takes as input =sign,epsI,epsA,chi and returns the angle Beta through which the xz coordinates rotate", metavar="FILE")

	# Model arguments
	parser.add_argument("-r", "--run" , help = Model.Run.__doc__)

	parser.add_argument("-pp", "--print_parameters",help = Useful_Tools.Print_Parameters.__doc__)

	# Plotting arguments

	parser.add_argument("-p","--plot",help=Plot.Simple_Plot.__doc__, metavar="FILE")

	parser.add_argument("-s","--splot",help= Plot.Spherical_Plot.__doc__, metavar="FILE")

	parser.add_argument("-a", "--alpha",help = Plot.Alpha_Plot.__doc__, metavar="FILE")

	parser.add_argument("-3", "--threeDplot", help = Plot.ThreeD_Plot_Cartesian.__doc__ ) 

	parser.add_argument("-g", "--angle_space", help = Plot.Angle_Space_Plot.__doc__ )

	parser.add_argument("-l", "--plot_beta_transform", help = Plot.Simple_Plot_Transform.__doc__ )

	parser.add_argument("-z", "--splot_beta_transform", help = Plot.Spherical_Plot_Transform.__doc__ )

	parser.add_argument("-n", "--param_space_plot",help = NLD_Functions.Parameter_Space_Plot.__doc__)


	# Additional arguments are passed to opts 
	parser.add_argument("-o","--options")

	# Set the verbose/quite options, this will be passed to the option dictionary
	parser.add_argument("-v", action="store_true", dest="verbose", default=True)
	parser.add_argument("-q", action="store_false", dest="verbose")

	arguments = parser.parse_args()
	return arguments

def main():

	arguments = parse_command_line(sys.argv)

	# Create the options dictionary 
	if arguments.options : Option_Dictionary = Create_Dictionary(arguments.options)
	else : Option_Dictionary = {}

	# Add the verbosity to the Option Dictionary
	Option_Dictionary['verbose'] = arguments.verbose

	if arguments.run : 
		Input_Dictionary = Create_Dictionary(arguments.run)
		Model.Run(Input_Dictionary)

	if arguments.print_parameters : Useful_Tools.Print_Parameters(arguments.print_parameters)

	if arguments.plot: Plot.Simple_Plot(arguments.plot,Option_Dictionary)

	if arguments.splot: Plot.Spherical_Plot(arguments.splot,Option_Dictionary)

	if arguments.alpha : Plot.Alpha_Plot(arguments.alpha,Option_Dictionary)

	if arguments.threeDplot :  Plot.ThreeD_Plot_Cartesian(arguments.threeDplot,Option_Dictionary)
	
	if arguments.angle_space : Plot.Angle_Space_Plot(arguments.angle_space,Option_Dictionary)

	if arguments.plot_beta_transform: Plot.Simple_Plot_Transform(arguments.plot_beta_transform,Option_Dictionary)

	if arguments.splot_beta_transform: Plot.Spherical_Plot_Transform(arguments.splot_beta_transform,Option_Dictionary)

	if arguments.param_space_plot: NLD_Functions.Parameter_Space_Plot(arguments.param_space_plot,Option_Dictionary)

if __name__ == "__main__":
    main()


