#!/usr/bin/python 
import optparse, sys, os
import numpy as np
import pylab as py 
import optparse, sys, os
from math import pi 
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d

# Import functions external files
import lib.NLD_Functions as NLD_Functions
import lib.GSL_ODE_Functions as GSL_ODE_Functions
import lib.File_Functions as File_Functions
import lib.Physics_Functions as Physics_Functions 
import lib.Plotting_Functions as Plotting_Functions




def main():
	def parse_command_line(argvs):

		parser = optparse.OptionParser()
	

		parser.add_option("-p","--plot",help="Plot the last three colums in FILE against the first column it is assumed this refer to wx,wy and wz", metavar="FILE")

		parser.add_option("-s","--splot",help="Plot the last three colums in FILE against the first column converting to spherical polars", metavar="FILE")



		parser.add_option("-c", "--Cplot", help ="Same as plot but transforming coordinates to principle axis of the effective moment of inertia tensor ", metavar="FILE")

		parser.add_option("-d", "--Csplot", help ="Same as splot but transforming coordinates to principle axis of the effective moment of inertia tensor", metavar="FILE")

		parser.add_option("-3", "--ThreeD", help ="Plot the input in 3d of input file ", metavar="FILE")

		parser.add_option("-a", "--alpha",help = "Plot the allignment of omega with the magnetic dipole axis ", metavar="FILE")
	
		parser.add_option("-x", "--param_space_plot", help = "Plot omega_dot, a_dot, phi_dot ", metavar="FILE")

		parser.add_option("-z", "--correlation_sum" , help ="Compute the correlated sum of input file", metavar="dotted_variables_XXX")

		parser.add_option("-r", "--run" , help = "Create the generic write file that can be compiled and run in C. The input should be of the form  chi,epsI,epsA,omega_0,eta,err=1.0e-12,*args  where \n chi = angle between magnetic dipole and z axis \n epsI = elastic deformation \n epsA magnetic deformation \n omega_0 = initial spin period \n eta = Simultion will stop when omega**2.0 < eta*omega_0**2.0 ")

		parser.add_option("-i", "--plot_angular_momentum_and_energy" ) 


		# Options to be called with above

		parser.add_option("-e", "--end_val", help = "Print the average end value from the given input data ", action='store_true')

		parser.add_option("-o", "--opts", help = "Options argument used in some particular tasks " ,metavar="o_1,o_2,..")

		parser.add_option("--save_data", help = "Save the data from the task in an appropriate way", action='store_true', metavar="FILE")

		parser.add_option("--save_fig" , help = "Save the figure in an appropriate way", action='store_true' ,metavar="FILE")

		(options,arguments) = parser.parse_args(argvs)
		return options, arguments

	options, arguments = parse_command_line(sys.argv)


	if options.b_2_e: nsmod_functions.Magnetic_Field_to_Epsilon_A(options)

	if options.e_2_b: nsmod_functions.Epsilon_A_to_Magnetic_Field(options)

	if options.plot: nsmod_functions.Simple_Plot(options)

	if options.splot: nsmod_functions.Spherical_Plot(options.splot,options.opts,)

	if options.Beta: nsmod_functions.Print_Beta(options)

	if options.Cplot: nsmod_functions.Simple_Plot_Transform(options.Cplot,options)

	if options.Csplot: nsmod_functions.Spherical_Plot_Transform(options.Csplot,options)

	if options.ThreeD : nsmod_functions.ThreeD_Plot_Cartesian(options)

	if options.alpha : nsmod_functions.Alpha_Plot(options)


	if options.run : nsmod_functions.Run(options.run,options.opts)

	if options.plot_angular_momentum_and_energy: nsmod_functions.Plot_Angular_Momentum_and_Energy(options)

if __name__ == "__main__":
    main()
