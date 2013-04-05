#!/usr/bin/python 
""" Some commonly used functions in plotting called from tools"""

import pylab as py 
from math import floor

# Generic functions 

def Sort_Out_Some_Axis(x):
	largest_value = max(x)
	largest_value_0sf = Round_To_n(largest_value,0)
        
    # This next line finds the order of magnitude of the maximum value so it can be appropriately scaled
	Scale = largest_value_0sf/float(str(largest_value_0sf)[0]) 
	Order_of_Magnitude = int(py.log10(Scale))
        
	x_scaled=[xi/Scale for xi in x]
	return (x_scaled,Order_of_Magnitude)
	

def Round_To_n(x,n):
	return round(x, -int(floor(py.log10(x)))+n)



def Plot_a_phi(time,a,phi):
	# Function to help scale the x-axis
	(scale_val,t_scaled) = sort_out_time_axis(time)	

	ax1 = py.subplot(111)

	# Plot aprime(t)
	ax1.plot(t_scaled,aprime)
	#py.axhline(90,ls="--",color="k")

	# Ploptions
	y_max = 105
	py.ylim(0,y_max)
	#py.yticks(fig2.get_yticks()[0:-2])
	py.yticks(py.arange(0,y_max,15))
	py.ylabel("$a' \;[^{\circ}]$",rotation="vertical")


	# Plot phiprime(t)
	ax2 = ax1.twinx()

	phiprime = Fix_Phi_Degrees(phiprime)
	if abs(phiprime[-1])>100:
		def Scale_Axis(axis):
			max_item = max(axis) ; min_item = min(axis)
			if abs(max_item) < abs(min_item): max_item = abs(min_item)
			scale = round_to_n(max_item,0)
			axis_scaled = [ ai/scale for ai in axis]
			return (axis_scaled , scale)

		(phiprime_scaled , scale) = Scale_Axis(phiprime)
		ax2.plot(t_scaled,phiprime_scaled)

		py.ylabel("$\phi' \; 10^{"+str(int(py.log10(scale)))+"} [^{\circ}]$",rotation="vertical")

	else : 
		ax2.plot(t_scaled,phiprime)
		py.ylabel("$\phi' \; [^{\circ}]$",rotation="vertical")
	#Ploptions
	
	py.xlabel(r"$t\;  [1\times 10^{"+str(scale_val)+"}s]$",fontsize=16)

	py.show()


