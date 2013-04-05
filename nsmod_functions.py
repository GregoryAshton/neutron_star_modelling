#!/usr/bin/python 

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

def Magnetic_Field_to_Epsilon_A(Bs):
	Bs=float(Bs)
	R = 1e6 #cm
	c=3e10 #cm/s
	I0=1e45
	m=0.5*Bs*pow(R,3)
	epsA=pow(m,2)/(I0*R*pow(c,2))
	print "Bs="+str(Bs)+" Epsilon_A="+str(epsA)

def Epsilon_A_to_Magnetic_Field(Epsilon_A):
	Epsilon_A = float(Epsilon_A)
	R = 1e6 #cm
	c=3e10 #cm/s
	I0=1e45
	m = py.sqrt(Epsilon_A*I0*R*pow(c,2))
	Bs = 2*m/pow(R,3)
	print "Bs="+str(Bs)+" Epsilon_A="+str(Epsilon_A)
