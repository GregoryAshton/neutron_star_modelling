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
#import lib.NLD_Functions as NLD_Functions
#import lib.GSL_ODE_Functions as ODE_Functions
import lib.File_Functions as File_Functions
import lib.Physics_Functions as Physics_Functions 
import lib.Plotting_Functions as Plotting_Functions

def Save_Figure(file_name,type_of_plot,format_type=".png"):
	plot_file_name = type_of_plot+"_"+file_name+format_type
	py.savefig(plot_file_name)
	print "Saving figure as %s" % plot_file_name

def Defaults():
	# Set the default font for all plots 
	from matplotlib import rc
	rc('font',**{'family':'serif','serif':['Computer MOption_Dictionaryern Roman']})
	rc('text', usetex=True)

	# Set the defaults for axis
	py.rcParams['axes.color_cycle'] = ['k', 'r', 'cyan']
	py.rcParams['font.size'] = 18
	py.rcParams[ 'axes.labelsize'] = 30	
	py.rcParams['lines.linewidth'] = 2
	py.rcParams['axes.grid']=True
	py.rcParams['figure.figsize']= (10.0, 8.0)
	#py.subplots_adjust(left=0.13, right=0.9, top=0.9, bottom=0.12,hspace=0.0)




def Beta_Function(epsI,epsA,chi):
	if abs(chi) < 2*pi : print "Please check that the chi used in Beta_Function is in degrees and not radians"
 
	if epsI>=0:
		sign=-1.0
	else: sign=1.0

	b=py.sqrt(epsA*epsA+epsI*epsI-2*epsA*epsI*py.cos(2*chi))
	return py.arctan((epsI+epsA*(1-2*pow(py.cos(chi),2.0))+sign*b)/(2*epsA*py.sin(chi)*py.cos(chi)))


# Plotting functions

def Simple_Plot(file_name,Option_Dictionary):
	""" Plots the given file as a function of time, it is assumed the columns of data are given by time , omega_x , omega_y and omega_z"""
		
	# Import the data 
	(time,x,y,z) = File_Functions.Import_Data(file_name)

	# Handle any additional options which are in the dictionary
	if Option_Dictionary.has_key('tmax') : tmax = Option_Dictionary['tmax'] 
	else : tmax = max(time)
	if Option_Dictionary.has_key('tmin') : tmin  = Option_Dictionary['tmin']
	else : tmin = min(time)

	fig1=py.subplot(3,1,1) ; fig1.set_xticklabels([])
	fig1.plot(time,x)
	py.ylabel(r"$\omega_{x}$",rotation="horizontal")
	py.yticks(fig1.get_yticks()[1:-1])
	py.xlim(tmin,tmax)

	fig2=py.subplot(3,1,2) ; fig2.set_xticklabels([])
	fig2.plot(time,y)
	py.yticks() 
	py.ylabel("$\omega_{y}$",rotation="horizontal")
	py.yticks(fig2.get_yticks()[1:-1])
	py.xlim(tmin,tmax)

	fig2=py.subplot(3,1,3)
	fig2.plot(time,z)
	py.ylabel("$\omega_{z} $",rotation="horizontal")
	py.xlim(tmin,tmax)

	py.xlabel(r"$t$")

	py.show()


def Spherical_Plot(file_name,Option_Dictionary):
	""" Plot the input data after transforming to spherical polar coordinates 
	The opts dictionary may contain 
	nmax ~ limit the data from 0:nmax
	tmax , tmin ~ limit the xaxis
	end_val ~ print the average of the last 100 points
	save_fig ~ 
"""
	# Default settings 
	labelx = -0.1  # x position of the yaxis labels

	# Handle any additional options which are in the dictionary
	if Option_Dictionary.has_key("nmax"):
		max_n = int(Option_Dictionary['nmax']) 
	else : max_n = -1

	(time,omega_x,omega_y,omega_z)= File_Functions.Import_Data(file_name,max_n) 

	if Option_Dictionary.has_key('tmax') : tmax = Option_Dictionary['tmax'] 
	else : tmax = max(time)
	if Option_Dictionary.has_key('tmin') : tmin  = Option_Dictionary['tmin']
	else : tmin = 0.0


	print "Number of points in plot = ",len(time)
	# Transform to spherical polar coordinates
	(omega,a,phi) = Physics_Functions.Transform_Cartesian_2_Spherical(omega_x,omega_y,omega_z)
	
	# Function to help scale the x-axis
	(t_scaled,scale_val) = Plotting_Functions.Sort_Out_Some_Axis(time)	

	# Plot omega(t)
	fig = py.figure()
	ax1=fig.add_subplot(3,1,1) ; ax1.set_xticklabels([])
	ax1.plot(t_scaled,omega)
	ax1.set_xlim(tmin*pow(10,-scale_val),tmax*pow(10,-scale_val))

	ax1.set_ylim(0,1.1*max(omega))
	#py.yticks(fig1.get_yticks()[1:-1])
	ax1.set_ylabel(r"$\omega$  [Hz] ",rotation="vertical")
	ax1.yaxis.set_label_coords(labelx, 0.5)

	# Plot a(t)
	ax2=fig.add_subplot(3,1,2) ; ax2.set_xticklabels([])
	ax2.plot(t_scaled,a)
	#py.axhline(90,ls="--",color="k")

	ax2.set_ylim(0,105)
	#py.yticks(fig2.get_yticks()[0:-2])
	ax2.set_yticks(py.arange(0,105,15))
	ax2.set_ylabel("$a $ [deg]",rotation="vertical")
	ax2.yaxis.set_label_coords(labelx, 0.5)
	ax2.set_xlim(tmin*pow(10,-scale_val),tmax*pow(10,-scale_val))

	# Plot phi(t)
	ax3=fig.add_subplot(3,1,3) 
	phi=Physics_Functions.Fix_Phi(phi) # By default we assume the phi is broken so fix it
	ax3.plot(t_scaled,phi)

	#Ploptions
	#ax3.set_ylim(0,110)
	#ax3.set_yticks(py.arange(0,105,15))
	ax3.set_yticks(ax3.get_yticks()[0:-1])
	ax3.set_ylabel("$\phi$ [deg]",rotation="vertical")
	ax3.yaxis.set_label_coords(labelx, 0.5)
	ax3.set_xlabel(r"time  [$1\times 10^{"+str(scale_val)+"}$ s]")
	ax3.set_xlim(tmin*pow(10,-scale_val),tmax*pow(10,-scale_val))
	if Option_Dictionary.has_key('end_val') :
		print " Data on the end value of the spherical components of omega"
		omega_end = omega[-100:-1]
		print " Average of |omega| :  %s s^-1  \n Range of omega : %s"	% (py.average(omega_end),max(omega_end)-min(omega_end))
		a_end = a[-100:-1]
		print " Average of a  : %s s^-1 \n Range of a : %s"	% (py.average(a_end),max(a_end)-min(a_end))
		phi_end = omega[-100:-1]
		print " Average of phi :  %s s^-1 \n Range of phi : %s"	% (py.average(phi_end),max(phi_end)-min(phi_end))

	py.subplots_adjust(left=0.13, right=0.9, top=0.9, bottom=0.12,hspace=0.0)
	if Option_Dictionary.has_key('save_fig') and Option_Dictionary['save_fig'] == True :
		Save_Figure(file_name,"Spherical_Plot")
	else:
		py.show()

def Alpha_Plot(file_name,Option_Dictionary):
	""" Plots the alignment of the input file [t,omega_x , omega_y and omega_z] against the magnetic dipole"""

	# Handle any additional options which are in the dictionary
	if Option_Dictionary.has_key("nmax"):
		max_n = int(Option_Dictionary['nmax']) 
	else : max_n = -1

	# Import the data in components x,y,z
	(time,omega_x,omega_y,omega_z)= File_Functions.Import_Data(file_name,max_n) 

	# Get the paramters of the run
	Parameter_Dictionary = File_Functions.Params_From_File_Name(file_name)


	# Extract some parameters from the first line in the file 
	epsI=Parameter_Dictionary["epsI"]
	epsA=Parameter_Dictionary["epsA"]
	chi=Parameter_Dictionary["chi"] # This will be returned in radians
	
	# Transform to spherical polar coordinates specifying that we wish the angles to be in Radians rather than degrees
	(omega,a,phi) = Physics_Functions.Transform_Cartesian_2_Spherical(omega_x,omega_y,omega_z,"Radians")
	
	# Function to help scale the t-axis
	(t_scaled,scale_val) = Plotting_Functions.Sort_Out_Some_Axis(time)	

	# Calculate the angle made with the magnetic dipole assumed to lie at chi to the z axis in the x-z plane
	chi_radians = chi * pi/180
	Sx=py.sin(chi_radians) ; Cx=py.cos(chi_radians)
	alpha = [py.arccos(Sx*py.sin(a[i])*py.cos(phi[i]) + Cx*py.cos(a[i]))*180/pi for i in range(len(a))]

	fig = py.figure()
	ax1=fig.add_subplot(111)
	ax1.plot(t_scaled,alpha,lw=2)
	ax1.set_ylabel(r"$\alpha$ [deg]",rotation="horizontal")
	ax1.set_xlabel(r"time  [$1\times 10^{"+str(scale_val)+"}$ s]")

	if Option_Dictionary.has_key('save_fig') and Option_Dictionary['save_fig'] == True :
		Save_Figure(file_name,"Alpha")
	else:
		py.show()

	if Option_Dictionary.has_key('end_val') :
		print " Average of the last "+str(5)+" points is "+str(py.average(alpha[-6:-1]))+" degrees"

def ThreeD_Plot_Cartesian(file_name,Option_Dictionary):
	""" 

	Plots the components of input file in 3D Option_Dictionary take the following 
	# Integers to retrict the number of plotted points
	start = int  
	stop = int 

	# Colouring commands
	power = float 
	"""

	# Import the data in components x,y,z
	(time,omega_x,omega_y,omega_z)= File_Functions.Import_Data(file_name) 

	# Get the paramters of the run
	Parameter_Dictionary = File_Functions.Params_From_File_Name(file_name)

	# Set the defaults and then overide if they exist in Option_Dictionary
	start = 0 
	stop = -1

	if Option_Dictionary.has_key("start") : start =  int(Option_Dictionary["start"])
	if Option_Dictionary.has_key("stop") : stop =  int(Option_Dictionary["stop"])

	if Option_Dictionary['verbose'] ==True :
		print 
		print " Reducing plotted data size from %s to %s" % (len(time) , stop-start)
		print " The observation time  is from t=%s to %s seconds" % ( time[start] , time[stop])
		print 
	
	# Reduce the number of points
	time = time[start:stop] ; omega_x = omega_x[start:stop] ; omega_y = omega_y[start:stop] ; omega_z = omega_z[start:stop]

	# Extract some parameters from the dictionary
	epsI=Parameter_Dictionary["epsI"]
	epsA=Parameter_Dictionary["epsA"]
	chi =Parameter_Dictionary["chi"]

	# Compute the rotation angle of primed axis and transform the solution omega_{xyz} into these coordinates
#	beta=Beta_Function(epsI,epsA,chi)

#	Cb=py.cos(beta) ; Sb=py.sin(beta)

#	xprime=[x[i]*Cb - z[i]*Sb for i in range(len(x))]
#	zprime=[z[i]*Cb + x[i]*Sb for i in range(len(x))]
#	yprime = y

	#fig = py.figure()
	x = omega_x ; y = omega_y ; z = omega_z

	ax = py.subplot(111, projection='3d')

	# Create and label the primed axis
#	ax.plot(py.zeros(100),py.zeros(100),py.linspace(-max(z),max(z),100),color="k")
#	ax.text(0,0,max(z)*1.1,"$z'$")
#	ax.plot(py.zeros(100),py.linspace(max(y),min(y),100),py.zeros(100),color="k")
#	ax.text(0,max(y)*1.1,0,"$y'$")
#	ax.plot(py.linspace(max(x),min(x),100),py.zeros(100),py.zeros(100),color="k")
#	ax.text(max(x)*1.1,0,0,"$x'$")

	# Compute same variables used for colouring and plot the x',y' and z' transforming the colour as time changes
	if Option_Dictionary.has_key('power') : 
		power = float(Option_Dictionary['power'])
#		n=len(x)
#		d = int(Plotting_Functions.Round_To_n(n,0))/100
#		s=n/d 	
#		for i in range(1,d-1):
#			ax.plot(x[s*i:s*i+s],y[s*i:s*i+s],z[s*i:s*i+s],color=(0.0,1-float(pow(i,power)*1.0/pow(d,power)),0.8),alpha=0.5)

		n = len(time)
		for i in range(0,n-10,10):
			ax.plot(x[i:i+11],y[i:i+11],z[i:i+11],color=(0.0,1-float(pow(i,power)*1.0/pow(n,power)),float(pow(i,power)*1.0/pow(n,power))),alpha=1.0)		
	else :
		ax.plot(x,y,z)



	#py.ylabel(r"$\omega_{y}'$")
	#py.zlabel(r"$\omega_{z}'$")

	ax.set_xticklabels([])
	ax.set_yticklabels([])
	ax.set_zticklabels([])

	ax.set_zlabel(r"$\hat{z}$",rotation='horizontal')
	ax.grid(False)

	if Option_Dictionary.has_key('save_fig') and Option_Dictionary['save_fig'] == True :
		Save_Figure(file_name,"Spherical_Plot")
	else:
		py.show()

#	if options.save_data:
#		file_name=options.save_data
#		write_file=open(file_name,"w+")
#		for i in range (len(xprime)):
#			write_file.write(str(time[i])+" "+str(xprime[i])+" "+str(yprime[i])+" "+str(zprime[i])+"\n")
#		write_file.close()


def Angle_Space_Plot(file_name,Option_Dictionary):
	"""Experimental! USed to plot the angular components against each other of the spin vector"""	
	
	# Handle any additional options which are in the dictionary
	if Option_Dictionary.has_key("nmax"):
		max_n = int(Option_Dictionary['nmax']) 
	else : max_n = -1
	# Import the data in components x,y,z
	(time,omega_x,omega_y,omega_z)= File_Functions.Import_Data(file_name,max_n) 

	# Transform to spherical coordinates
	(omega,a,phi) = Physics_Functions.Transform_Cartesian_2_Spherical(omega_x,omega_y,omega_z)

	phi=Physics_Functions.Fix_Phi(phi) # By default we assume the phi is broken so fix it

	fig = py.figure()
	ax1 = fig.add_subplot(111)
	n = len(time)
	list_vals = py.linspace(0,1,n)
	#for i in range(0,n-12,10):
	#	ax1.plot(phi[i:i+11],a[i:i+11],color=(list_vals[i],0.0,list_vals[i]))
	ax1.plot(phi,a)
	ax1.set_xlabel("$\phi$ [deg]")
	ax1.set_ylabel("$a$ [deg]")	
	#ax1.set_xlim(0,180)
	#ax1.set_ylim(0,180)

	if Option_Dictionary.has_key('arrow'):	
		if Option_Dictionary['arrow'] == True :
			for i in range(0,n,int(n/20.0)):
				py.arrow(phi[i],a[i],0.1*(phi[i+1]-phi[i]),0.1*(a[i+1]-a[i]),color="k",fill=True, head_width=1.0,head_starts_at_zero=True,alpha=0.8)

		else :
			try :			
				i = int(Option_Dictionary['arrow'])
				py.arrow(phi[i],a[i],0.1*(phi[i+1]-phi[i]),0.1*(a[i+1]-a[i]),color="k",fill=True, head_width=1.5,head_starts_at_zero=True,alpha=0.8)
			except ValueError :
				pass
	
	if Option_Dictionary.has_key('save_fig') and Option_Dictionary['save_fig'] == True :
		Save_Figure(file_name,"Angle_Space_Plot")
	else:
		py.show()
def main():

	# A function which translates the options.opts input into a dictionary of parameters
	def Create_Option_Dictionary(opts):
		Option_Dictionary={}
		for item in opts.split(","):
			if "=" in item: Option_Dictionary[item.split("=")[0]] = float(item.split("=")[1])
		
		else : Option_Dictionary[item] = True
		return Option_Dictionary

	def parse_command_line(argvs):

		parser = optparse.OptionParser()
	

		parser.add_option("-p","--plot",help=Simple_Plot.__doc__, metavar="FILE")

		parser.add_option("-s","--splot",help= Spherical_Plot.__doc__, metavar="FILE")


#		parser.add_option("-c", "--Cplot", help ="Same as plot but transforming coordinates to principle axis of the effective moment of inertia tensor ", metavar="FILE")

#		parser.add_option("-d", "--Csplot", help = , metavar="FILE")

		parser.add_option("-a", "--alpha",help = Alpha_Plot.__doc__, metavar="FILE")

		parser.add_option("-3", "--threeDplot", help = ThreeD_Plot_Cartesian.__doc__ ) 

		parser.add_option("-g", "--angle_space", help = Angle_Space_Plot.__doc__)

		parser.add_option("-o", "--opts")

		# Set the verbose/quite options, this will be passed to the option dictionary
		parser.add_option("-v", action="store_true", dest="verbose", default=True)
		parser.add_option("-q", action="store_false", dest="verbose")


		(options,arguments) = parser.parse_args(argvs)
		return options, arguments

	options, arguments = parse_command_line(sys.argv)

	if options.opts : Option_Dictionary = Create_Option_Dictionary(options.opts)
	else : Option_Dictionary = {}

	# Add the verbosity to the Option Dictionary
	Option_Dictionary['verbose'] = options.verbose

	if options.plot: Simple_Plot(options.plot,Option_Dictionary)

	if options.splot: Spherical_Plot(options.splot,Option_Dictionary)

	if options.alpha : Alpha_Plot(options.alpha,Option_Dictionary)

	if options.threeDplot :  ThreeD_Plot_Cartesian(options.threeDplot,Option_Dictionary)
	
	if options.angle_space : Angle_Space_Plot(options.angle_space,Option_Dictionary)

#	if options.Cplot: Simple_Plot_Transform(options.Cplot,options)

#	if options.Csplot: Spherical_Plot_Transform(options.Csplot,options)




if __name__ == "__main__":
	Defaults()
	main()
