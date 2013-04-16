#!/usr/bin/python 

import numpy as np
import pylab as py 
from math import pi


def Transform_Cartesian_2_Spherical(x,y,z,Angle_Type="Degrees"):

	if "Degrees" in Angle_Type:
		""" Transform x,y,z to spherical coordinates returning omega,a,phi in degrees. This is the default"""
		print "Transform using angle type Degrees"
		N=len(x)
		radial=[(x[i]*x[i]+y[i]*y[i]+z[i]*z[i])**0.5 for i in range(N)]
		polar=[py.arccos(z[i]/radial[i])*180/pi for i in range(N)]
		azimuth=[py.arctan(y[i]/x[i])*180/pi for i in range(N)]
		return (radial,polar,azimuth)

	elif Angle_Type in ["Radians","Radian","Rads"]:
		""" Transform x,y,z to spherical coordinates returning omega,a,phi in radians"""
		print "Transform using angle type Radians"
		N=len(x)
		radial=[(x[i]*x[i]+y[i]*y[i]+z[i]*z[i])**0.5 for i in range(N)]
		polar=[py.arccos(z[i]/radial[i]) for i in range(N)]
		azimuth=[py.arctan(y[i]/x[i]) for i in range(N)]
		return (radial,polar,azimuth)

def Transform_Cartesian_Body_Frame_2_Effective_Body_Frame(x,y,z,Beta):
	"""Functions to rotate x,z system by angle Beta w.r.t z axis"""
	Cb=py.cos(Beta) ; Sb=py.sin(Beta)
	x_prime=[x[i]*Cb - z[i]*Sb for i in range(len(x))]
	z_prime=[z[i]*Cb + x[i]*Sb for i in range(len(x))]
	y_prime = y
	return (x_prime,y_prime,z_prime)

def Rotational_Kinetic_Energy(Ix,Iy,Iz,omega_x,omega_y,omega_z):
	return 0.5*(Ix*pow(omega_x,2)+Iy*pow(omega_y,2)+Iz*pow(omega_z,2))

def Fix_Phi(phi,epsilon=170.0,Angle_Type="Degrees"):
	""" Takes a list of phi values and looks for jumps greater than epsilon, assuming these reflect a full rotation (2pi-0) it adds a correction factor to the subsequent data"""
	
	if Angle_Type == "Degrees" :
		phi_fix=[] ; fix = 0.0 ;  phi_fix.append(phi[0])
		for i in range(1,len(phi)):
			if abs(phi[i]-phi[i-1])> epsilon:
				fix+= -1*py.sign(phi[i]-phi[i-1])*180.0
			phi_fix.append(phi[i]+fix)
		phi=phi_fix
		return phi

	elif Angle_Type in ["Radians","Radian","Rads"] :
		# Change the default value of epsilon
		if epsilon==170.0 : epsilon = 1.0 

		phi_fix=[] ; fix = 0.0 ;  phi_fix.append(phi[0])
		for i in range(1,len(phi)):
			if abs(phi[i]-phi[i-1])> epsilon:
				fix+= -1*py.sign(phi[i]-phi[i-1])*pi
			phi_fix.append(phi[i]+fix)
		phi=phi_fix
		return phi

