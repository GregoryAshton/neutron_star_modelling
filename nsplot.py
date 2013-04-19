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
from lib.Physics_Functions import Beta_Function #?

def Save_Figure(file_name,type_of_plot,format_type=".png"):
	plot_file_name = type_of_plot+"_"+file_name.rstrip(".txt")+format_type
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



# OLD beta function
#def Beta_Function(epsI,epsA,chi):
#	if abs(chi) < 2*pi : print "Please check that the chi used in Beta_Function is in degrees and not radians"
# 
#	if epsI>=0:
#		sign=-1.0
#	else: sign=1.0

#	b=py.sqrt(epsA*epsA+epsI*epsI-2*epsA*epsI*py.cos(2*chi))
#	return py.arctan((epsI+epsA*(1-2*pow(py.cos(chi),2.0))+sign*b)/(2*epsA*py.sin(chi)*py.cos(chi)))

### Beta function
#def Beta_Function(epsI,epsA,chi):
#	if chi>2*pi :
#		print "Assuming chi has been given in degrees rather than radians, we now transform"
#		chi = chi*pi/180
#	a=epsA*epsA+epsI*epsI-2*epsA*epsI*py.cos(2*chi)
#	return py.arctan((epsI-epsA*py.cos(2*chi)-py.sqrt(a))/(2*epsA*py.sin(chi)*py.cos(chi)))

# Plotting functions

def Simple_Plot(file_name,Option_Dictionary):
	""" Plots the given file as a function of time, it is assumed the columns of data are given by time , omega_x , omega_y and omega_z"""

	# Handle any additional options which are in the dictionary
	if Option_Dictionary.has_key("nmax"):
		max_n = int(Option_Dictionary['nmax']) 
	else : max_n = -1
		
	# Import the data 
	(time,x,y,z) = File_Functions.Import_Data(file_name,max_n)

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
	nmax=int limit the data from 0:nmax
	tmax=float, tmin=float ~ limit the xaxis
	end_val=True ~ print the average of the last 100 points
	save_fig=True ~ saves the figure
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
	phi=Physics_Functions.Fix_Phi(phi,epsilon=70) # By default we assume the phi is broken so fix it
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
	
#	ax = py.subplot(111)
#	chi = float(file_name.split("_")[3])
#	R = [(py.sin(a[i])*py.cos(a[i])*py.sin(chi)*py.cos(chi)-pow(py.cos(a[i])*py.sin(chi),2)) for i in range(len(a))]
#	ax.plot(R)
#	py.show()

def Alpha_Plot(file_name,Option_Dictionary):
	""" Plots the alignment of the input file [t,omega_x , omega_y and omega_z] against the magnetic dipole"""

	# Handle any additional options which are in the dictionary
	if Option_Dictionary.has_key("nmax"):
		max_n = int(Option_Dictionary['nmax']) 
	else : max_n = -1

	# Import the data in components x,y,z
	(time,omega_x,omega_y,omega_z)= File_Functions.Import_Data(file_name,max_n) 

	# Get the paramters of the run
	Parameter_Dictionary = File_Functions.Parameter_Dictionary(file_name)


	# Extract some parameters from the first line in the file 
	epsI=float(Parameter_Dictionary["epsI"])
	epsA=float(Parameter_Dictionary["epsA"])
	chi=float(Parameter_Dictionary["chi"])*pi/180 # This will be returned in radians
	
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
	ax1.set_ylabel(r"$\alpha$ [deg]",rotation="vertical")
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
	Parameter_Dictionary = File_Functions.Parameter_Dictionary(file_name)

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
	epsI=float(Parameter_Dictionary["epsI"])
	epsA=float(Parameter_Dictionary["epsA"])
	chi=float(Parameter_Dictionary["chi"])*pi/180 # This will be returned in radians

	# Compute the rotation angle of primed axis and transform the solution omega_{xyz} into these coordinates
#	beta=Beta_Function(epsI,epsA,chi)

#	Cb=py.cos(beta) ; Sb=py.sin(beta)

#	xprime=[x[i]*Cb - z[i]*Sb for i in range(len(x))]
#	zprime=[z[i]*Cb + x[i]*Sb for i in range(len(x))]
#	yprime = y

	#fig = py.figure()
	x = omega_x ; y = omega_y ; z = omega_z

	ax = py.subplot(111, projection='3d')
	elev = 15.0
	azim = -134.0
	ax.view_init(elev, azim) 

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

	#ax.set_xticklabels([])
	#ax.set_yticklabels([])
	#ax.set_zticklabels([])

	ax.set_xlabel(r"$\hat{x}$",rotation='horizontal')
	ax.set_ylabel(r"$\hat{y}$",rotation='horizontal')
	ax.set_zlabel(r"$\hat{z}$",rotation='horizontal')
	ax.grid(False)

	if Option_Dictionary.has_key('save_fig') and Option_Dictionary['save_fig'] == True :
		Save_Figure(file_name,"ThreeD_Plot_Cartesian")
	else:
		py.show()

#	if options.save_data:
#		file_name=options.save_data
#		write_file=open(file_name,"w+")
#		for i in range (len(xprime)):
#			write_file.write(str(time[i])+" "+str(xprime[i])+" "+str(yprime[i])+" "+str(zprime[i])+"\n")
#		write_file.close()


def Angle_Space_Plot(file_name,Option_Dictionary):
	"""Experimental! Used to plot the angular components against each other of the spin vector. Option_Dictionary takes the following as arguments 
	nmax : int ~ take only the first nmax points from the file
	2D : True ~ This will plot the angular components phi and a in normal plot
	3D : True ~ This will plot the angular components phi and a projected onto the unit sphere
		split=n1/n2/n3/n4/... For use with 3D this will segment the data into chunks starting and stopping at the integers n1,n2...  
	"""	
	
	# Handle any additional options which are in the dictionary
	if Option_Dictionary.has_key("nmax"):
		max_n = int(Option_Dictionary['nmax']) 
	else : max_n = -1
	# Import the data in components x,y,z
	(time,omega_x,omega_y,omega_z)= File_Functions.Import_Data(file_name,max_n) 

	if Option_Dictionary.has_key('beta'):
		Parameter_Dictionary = File_Functions.Parameter_Dictionary(file_name)
		epsI=float(Parameter_Dictionary["epsI"])
		epsA=float(Parameter_Dictionary["epsA"])
		chi=float(Parameter_Dictionary["chi"])*pi/180 # This will be returned in radians
		beta=Beta_Function(epsI,epsA,chi)
		(omega_x,omega_y,omega_z) = Physics_Functions.Transform_Cartesian_Body_Frame_2_Effective_Body_Frame(omega_x,omega_y,omega_z,beta)
		# Note we have chosen to the omega_x...here when this is actually in the primed coordinates. This is to avoid confusion but should we perhaps use some generic names instead?
		
	if Option_Dictionary.has_key('2D'):

		# Transform to spherical coordinates
		(omega,a,phi) = Physics_Functions.Transform_Cartesian_2_Spherical(omega_x,omega_y,omega_z)

		phi=Physics_Functions.Fix_Phi(phi,epsilon=70) # By default we assume the phi is broken so fix it

		fig = py.figure()
		ax1 = fig.add_subplot(111)
		n = len(time)
		list_vals = py.linspace(0,1,n)

		ax1.plot(phi,a)
		if Option_Dictionary.has_key('beta'):
			ax1.set_xlabel("$\phi'$ [deg]")
			ax1.set_ylabel("$a'$ [deg]")	
		else:
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
			Save_Figure(file_name,"Angle_Space_Plot_2D")
		else:
			py.show()

	if Option_Dictionary.has_key('3D'):
		
		ax = py.subplot(111, projection='3d')
		# Set the viewing position
		elev = 15.0
		azim = -134.0
		ax.view_init(elev, azim) 

	
		# Transform to spherical coordinates
		(omega,a,phi) = Physics_Functions.Transform_Cartesian_2_Spherical(omega_x,omega_y,omega_z,Angle_Type="Radians")

		phi=Physics_Functions.Fix_Phi(phi,Angle_Type="Radians") # By default we assume the phi is broken so fix it	

		
		if Option_Dictionary.has_key('split'):
			# Import the values to split by
			values=[int(it) for it in Option_Dictionary['split'].split("/")]
		
			# Define color library for the splits
			colors = ["b","r","k","c"]
			for j in range(0,len(values),2):
				# Trajectory of path
				low = values[j] ; high = values[j+1] 
				x=[1.0*py.sin(a[i])*py.cos(phi[i]) for i in range(low,high)]
				y=[1.0*py.sin(a[i])*py.sin(phi[i]) for i in range(low,high)]
				z=[1.0*py.cos(a[i]) for i in range(low,high)]
				elevation = 25.0
				azimuth = 25.0
				Plotting_Functions.ThreeD_Sphere(ax,elevation,azimuth,x,y,z,ls="-",lw=0.8,color=colors[j/2])
				

		else :	
			# Trajectory of path
			x=[1.0*py.sin(a[i])*py.cos(phi[i]) for i in range(len(time))]
			y=[1.0*py.sin(a[i])*py.sin(phi[i]) for i in range(len(time))]
			z=[1.0*py.cos(a[i]) for i in range(len(time))]

			elevation = 25.0
			azimuth = 25.0
			Plotting_Functions.ThreeD_Sphere(ax,elevation,azimuth,x,y,z,ls="-",lw=0.8,color="k")

		# Sphere of unit radius
		u = np.linspace(0, 2 * np.pi , 100)
		v = np.linspace(0, np.pi, 100)
		x = 1.0 * np.outer(np.cos(u), np.sin(v))
		y = 1.0 * np.outer(np.sin(u), np.sin(v))
		z = 1.0 * np.outer(np.ones(np.size(u)), np.cos(v))
		
		ax.plot_surface(x, y,z,rstride=4, cstride=9, color='white',alpha=0.9,lw=0.1)	

		if Option_Dictionary.has_key('beta'):
			ax.set_xlabel(r"$e_{1}$")
			ax.set_ylabel(r"$e_{2}$")
			ax.set_zlabel(r"$e_{3}$",rotation="horizontal")
			# These axis labels may prove to be wrong in the case of epsI<0
		else :
			ax.set_xlabel(r"$\hat{\omega_{x}}$")
			ax.set_ylabel(r"$\hat{\omega_{y}}$")
			ax.set_zlabel(r"$\hat{\omega_{z}}$")

		ax.set_xticklabels([])
		ax.set_yticklabels([])
		ax.set_zticklabels([])


		ax.grid(False)

		if Option_Dictionary.has_key('save_fig') and Option_Dictionary['save_fig'] == True :
			Save_Figure(file_name,"Angle_Space_Plot_3D")
		else:
			py.show()

def Simple_Plot_Transform(file_name,Option_Dictionary):
	""" Same as Simple_Plot() except transform to the effective MOI tensor principle axis, this has limited functionality as it is generally used for checking data is correct and not final plots """
	
	# Handle any additional options which are in the dictionary
	if Option_Dictionary.has_key("nmax"):
		max_n = int(Option_Dictionary['nmax']) 
	else : max_n = -1
	# Import the data in components x,y,z
	(time,omega_x,omega_y,omega_z)= File_Functions.Import_Data(file_name,max_n) 

	# Get the paramters of the run
	Parameter_Dictionary = File_Functions.Parameter_Dictionary(file_name)

	epsI=float(Parameter_Dictionary["epsI"])
	epsA=float(Parameter_Dictionary["epsA"])
	chi=float(Parameter_Dictionary["chi"])*pi/180 # This will be returned in radians
	beta=Beta_Function(epsI,epsA,chi)

	(omega_x_prime,omega_y_prime,omega_z_prime) = Physics_Functions.Transform_Cartesian_Body_Frame_2_Effective_Body_Frame(omega_x,omega_y,omega_z,beta)


	fig1=py.subplot(3,1,1)
	fig1.plot(time,omega_x_prime)
	py.ylabel(r"$\omega_{x}' $",fontsize=20,rotation="horizontal")

	fig2=py.subplot(3,1,2)
	fig2.plot(time,omega_y_prime)
	py.yticks() 
	py.ylabel("$\omega_{y}' $",fontsize=20,rotation="horizontal")


	fig2=py.subplot(3,1,3)
	fig2.plot(time,omega_z_prime)
	py.xlabel(r"$t$",fontsize=20)
	py.ylabel("$\omega_{z}' $",fontsize=20,rotation="horizontal")
	py.subplots_adjust(left=0.13, right=0.9, top=0.9, bottom=0.12,hspace=0.0)
	py.show()

def Spherical_Plot_Transform(file_name,Option_Dictionary):
	""" Plot the input data after transforming to spherical polar coordinates in the primed coordinates 
	options should be passed as a dictionary with the following arguments:
	opts : max_len ~ A maximum number of points to include
	end_val : True ~ Prints an average of the last 100 points in the plot
	save_fig : True ~ Save the figure in an appropriate way """


	# Default settings 
	labelx = -0.1  # x position of the yaxis labels

	# Handle any additional options which are in the dictionary
	if Option_Dictionary.has_key("nmax"):
		max_n = int(Option_Dictionary['nmax']) 
	else : max_n = -1
	# Import the data in components x,y,z
	(time,omega_x,omega_y,omega_z)= File_Functions.Import_Data(file_name,max_n) 

	# Get the paramters of the run
	Parameter_Dictionary = File_Functions.Parameter_Dictionary(file_name)

	epsI=float(Parameter_Dictionary["epsI"])
	epsA=float(Parameter_Dictionary["epsA"])
	chi=float(Parameter_Dictionary["chi"])*pi/180 # This will be returned in radians
	beta=Beta_Function(epsI,epsA,chi)
	print 	
	print "Beta  = %s degrees for %s" % (beta*180/pi,file_name)

	(omega_x_prime,omega_y_prime,omega_z_prime) = Physics_Functions.Transform_Cartesian_Body_Frame_2_Effective_Body_Frame(omega_x,omega_y,omega_z,beta)

	# Transform to the spherical polar coordinates
	(omega_prime , a_prime ,phi_prime) = Physics_Functions.Transform_Cartesian_2_Spherical(omega_x_prime,omega_y_prime,omega_z_prime)
	

	if Option_Dictionary.has_key('no_omega'):

			# Function to help scale the x-axis
		(t_scaled,scale_val) = Plotting_Functions.Sort_Out_Some_Axis(time)	

		fig=py.figure()

		ax1 = fig.add_subplot(111)

		# Plot a_prime(t)
		ax1.plot(t_scaled,a_prime,lw=1.0)
		#py.axhline(90,ls="--",color="k")

		# Ploptions
		y_max = 105
		py.ylim(0,y_max)
		#py.yticks(fig2.get_yticks()[0:-2])
		py.yticks(py.arange(0,y_max,15))
		py.ylabel("$a' \;[^{\circ}]$",rotation="horizontal",fontsize=18)


		# Plot phi_prime(t)
		ax2 = ax1.twinx()

		phi_prime = Physics_Functions.Fix_Phi(phi_prime)
		if abs(phi_prime[-1])>100:
			def Scale_Axis(axis):
				max_item = max(axis) ; min_item = min(axis)
				if abs(max_item) < abs(min_item): max_item = abs(min_item)
				scale = round_to_n(max_item,0)
				axis_scaled = [ ai/scale for ai in axis]
				return (axis_scaled , scale)

			(phi_prime_scaled , scale) = Scale_Axis(phi_prime)
			ax2.plot(t_scaled,phi_prime_scaled,color="b",lw=1.0)

			py.ylabel(r"$\phi' \; 1\times 10^{"+str(int(py.log10(scale)))+"} [^{\circ}]$",rotation="vertical")

		else : 
			ax2.plot(t_scaled,phi_prime,color="b",lw=1.0)
			py.ylabel("$\phi' \; [^{\circ}]$",rotation="vertical")
		#Ploptions
		ax2.tick_params(axis='y', colors='blue')
		ax2.yaxis.label.set_color('blue')
		py.xlabel(r"time  [$1\times 10^{"+str(scale_val)+"}$ s]",fontsize=16)
	
		


	
	else : 

		# Function to help scale the x-axis
		(t_scaled , scale_val) = Plotting_Functions.Sort_Out_Some_Axis(time)	

		# Plot omega_prime(t)
		fig = py.figure()		
		ax1 = fig.add_subplot(311) 
		ax1.set_xticklabels([])
		ax1.plot(t_scaled,omega_prime)

		ax1.set_ylim(0,1.1*max(omega_prime))
		#py.yticks(fig1.get_yticks()[1:-1])
		ax1.set_ylabel(r"$\omega'$ [Hz]",rotation="vertical")
		ax1.yaxis.set_label_coords(labelx, 0.5)

		# Plot a_prime(t)
		ax2 = fig.add_subplot(312) 
		ax2.set_xticklabels([])
		ax2.plot(t_scaled,a_prime)
		#py.axhline(90,ls="--",color="k")

		y_max = 120
		ax2.set_ylim(0,y_max)
		#py.yticks(fig2.get_yticks()[0:-2])
		ax2.set_yticks(py.arange(0,y_max,15))
		ax2.set_ylabel("$a'$ [deg]",rotation="vertical")
		ax2.yaxis.set_label_coords(labelx, 0.5)

		# Plot phi_prime(t)
		ax3 = fig.add_subplot(313) 

		# Check and fix rotations of 2pi in phi
		phi_prime = Physics_Functions.Fix_Phi(phi_prime)

		# Often phi becomes very large in which case we scale the axis
		if abs(phi_prime[-1])>1000:

			(phi_prime_scaled , scale) = Plotting_Functions.Sort_Out_Some_Axis(phi_prime)
			ax3.plot(t_scaled,phi_prime_scaled)
			ax3.set_ylabel(r"$\phi' [\;1\times 10^{"+str(int(py.log10(scale)))+"} $deg]",rotation="vertical")

		else : 
			ax3.plot(t_scaled,phi_prime)
			ax3.set_ylabel("$\phi'$  [deg]",rotation="vertical")

		#Ploptions
		ax3.set_xlabel(r"time  [$1\times 10^{"+str(scale_val)+"}$ s]",fontsize=16)		
		ax3.yaxis.set_label_coords(labelx, 0.5)
		ax3.set_yticks(ax3.get_yticks()[0:-1])

	if Option_Dictionary.has_key('end_val') :
		print " Data on the end value of the spherical components of omega"
		omega_end = omega_prime[-100:-1]
		print " Average of |omega| :  %s s^-1  \n Range of omega : %s"	% (py.average(omega_end),max(omega_end)-min(omega_end))
		a_end = a_prime[-100:-1]
		print " Average of a  : %s degrees \n Range of a : %s"	% (py.average(a_end),max(a_end)-min(a_end))
		phi_end = phi_prime[-100:-1]
		print " Average of phi :  %s  degrees \n Range of phi : %s"	% (py.average(phi_end),max(phi_end)-min(phi_end))

	py.subplots_adjust(left=0.13, right=0.9, top=0.9, bottom=0.12,hspace=0.0)
	if Option_Dictionary.has_key('save_fig'):
		Save_Figure(file_name,"Spherical_Plot_Transform")
	else : 
		py.show()

def main():

	# A function which translates the options.opts input into a dictionary of parameters
	def Create_Option_Dictionary(opts):
		Option_Dictionary={}
		for item in opts.split(","):
			if "=" in item: Option_Dictionary[item.split("=")[0]] = item.split("=")[1]
			else : Option_Dictionary[item] = True
		return Option_Dictionary
	
	def parse_command_line(argvs):

		parser = optparse.OptionParser()
	

		parser.add_option("-p","--plot",help=Simple_Plot.__doc__, metavar="FILE")

		parser.add_option("-s","--splot",help= Spherical_Plot.__doc__, metavar="FILE")

		parser.add_option("-a", "--alpha",help = Alpha_Plot.__doc__, metavar="FILE")

		parser.add_option("-3", "--threeDplot", help = ThreeD_Plot_Cartesian.__doc__ ) 

		parser.add_option("-g", "--angle_space", help = Angle_Space_Plot.__doc__ )

		parser.add_option("-l", "--plot_beta_transform", help = Simple_Plot_Transform.__doc__ )

		parser.add_option("-z", "--splot_beta_transform", help = Spherical_Plot_Transform.__doc__ )

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

	if options.plot_beta_transform: Simple_Plot_Transform(options.plot_beta_transform,Option_Dictionary)

	if options.splot_beta_transform: Spherical_Plot_Transform(options.splot_beta_transform,Option_Dictionary)




if __name__ == "__main__":
	Defaults()
	main()
