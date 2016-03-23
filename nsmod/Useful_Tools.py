#!/usr/bin/python
""" Some useful functions"""

import pylab as py
import numpy as np
from math import floor


def Texify_Float(f, n=1, power=True):
    """

    Takes a float and returns a string that looks nice in Latex, takes args

    n=int ~ Number of sig.fig
    power = Bool ~ Whether to produce a *10^{x} or just a regular number"""
    f = float(f)
    if f == 0:
        return "0"
    sign = [r"\text{-}", "", ""][1+int(np.sign(f))]
    f = abs(f)
    f = float(f)
    s = "{:1.2e}".format(f)
    f_nodecimal = float(s[:2] + "e" + s.split("e")[1])
    f_power = int(np.floor(np.log10(f_nodecimal)))
    if power and f_power not in [-1, 0, 1, 2]:
        f_SF = Round_To_n(f, n) * pow(10, -f_power)
        if f_SF < 1:
            f_SF = Round_To_n(f, n) * pow(10, -f_power + 1)
        if n == 0:
            f_SF = int(f_SF)
        return r"{{{}{}}}\times 10^{{{}}}".format(sign, f_SF, f_power)
    else:
        f_SF = Round_To_n(f, n)
        return sign + str(f_SF)

#def Texify_Float(f, n=1, power=True):
#    """
#
#    Takes a float and returns a string that looks nice in Latex, takes arguments
#
#    n=int ~ Number of sig.fig
#    power = Bool ~ Whether to produce a *10^{x} or just a regular number"""
#
#    f = float(f)
#    f_power = int(py.log10(Round_To_n(f, 0)))
#
#    if power and f_power not in  [-1, 0, 1, 2]:
#        f_SF = Round_To_n(f, n) * pow(10, -f_power)
#        return r" %s\times 10^{%s} " % (f_SF, f_power)
#    else:
#        f_SF = Round_To_n(f, n)
#        return  str(f_SF)

# This is not working at the moment, need to think some more....
def rescale_axis(ax, axis='x', label="time"):
    """

    Replacement of Sort_Out_Some_Axis

    """
    
    largest_value = ax.get_xlims()[1]
    scale = Round_To_n(largest_value, 0)

    order_of_magnitude = int(py.log10(scale))

    x_scaled = [xi / scale for xi in x]

    
    return (x_scaled, Order_of_Magnitude)


def Sort_Out_Some_Axis(x):
    """ Scales an axis appropriately returns a list of the scaled axis
        and the order of magnitude through which it has been scaled

    """

    largest_value = max(x)
    largest_value_0sf = Round_To_n(largest_value, 0)

    # This next line finds the order of magnitude of the maximum value
    Scale = largest_value_0sf / float(str(largest_value_0sf)[0])
    Order_of_Magnitude = int(py.log10(Scale))

    x_scaled = [xi / Scale for xi in x]
    return (x_scaled, Order_of_Magnitude)


def Round_To_n(x, n):
    return round(x, -int(np.floor(np.sign(x) * np.log10(abs(x)))) + n)


def Print_Parameters(file_name):
    """ Print a list of parameters about the file to the terminal """
    from lib.File_Functions import Parameter_Dictionary
    Parameter_Dictionary = Parameter_Dictionary(file_name)
    for key in Parameter_Dictionary:
        # If value is float lets make it look nice when it is printed
        try:
            value = float(Parameter_Dictionary[key])
            if value > 1000.0:
                value = Texify_Float(value, n=3, power=True)
            else:
                value = Texify_Float(value, n=3, power=False)
        except ValueError:
            value = Parameter_Dictionary[key]
        print " %s = %s " % (key, value)

# Obsolete?
#def Plot_a_phi(time, a, phi):
    ## Function to help scale the x-axis
    #(scale_val,t_scaled) = sort_out_time_axis(time)

    #ax1 = py.subplot(111)

    ## Plot aprime(t)
    #ax1.plot(t_scaled,aprime)
    ##py.axhline(90,ls="--",color="k")

    ## Ploptions
    #y_max = 105
    #py.ylim(0,y_max)
    ##py.yticks(fig2.get_yticks()[0:-2])
    #py.yticks(py.arange(0,y_max,15))
    #py.ylabel("$a' \;[^{\circ}]$",rotation="vertical")


    ## Plot phiprime(t)
    #ax2 = ax1.twinx()

    #phiprime = Fix_Phi_Degrees(phiprime)
    #if abs(phiprime[-1])>100:
        #def Scale_Axis(axis):
            #max_item = max(axis) ; min_item = min(axis)
            #if abs(max_item) < abs(min_item): max_item = abs(min_item)
            #scale = round_to_n(max_item,0)
            #axis_scaled = [ ai/scale for ai in axis]
            #return (axis_scaled , scale)

        #(phiprime_scaled , scale) = Scale_Axis(phiprime)
        #ax2.plot(t_scaled,phiprime_scaled)

        #py.ylabel("$\phi' \; 10^{"+str(int(py.log10(scale)))+"} [^{\circ}]$",rotation="vertical")

    #else :
        #ax2.plot(t_scaled,phiprime)
        #py.ylabel("$\phi' \; [^{\circ}]$",rotation="vertical")
    ##Ploptions

    #py.xlabel(r"$t\;  [1\times 10^{"+str(scale_val)+"}s]$",fontsize=16)

    #py.show()


def ThreeD_Sphere(axis, elevation, azimuth, x, y, z,
                    color="b", ls=".", lw=1, delta=1.0):
    """Function which given an axis with 3D projection, the azimuth and
       elevation plots the x,y,z shading the values which are on the opposite
       side of the unit sphere to the viewer.

    Inputs
    axis ~ an axis instance
    elevation ~
    axiumuth ~ view positions in degrees
    x,y,z 3 N length arrays to be plotted


    """

    elevation = np.radians(elevation)
    azimuth = np.radians(azimuth)

    # Init. lists for the front and back
    front_x = []
    front_y = []
    front_z = []
    back_x = []
    back_y = []
    back_z = []

    def angle_between_vectors(a, b):
        """ Returns angle between a and b """
        angle = py.arccos(sum(py.array(a) * py.array(b)) /
                     (py.norm(a) * py.norm(b))) * 180 / py.pi
        return angle

    # Define vector of the viewing position
    elevation = np.radians(90.0)  # Tempory fix
    view_pos = [py.sin(elevation) * py.cos(azimuth),
                py.sin(elevation) * py.sin(azimuth),
                py.cos(elevation)]

    for i in xrange(len(x)):
        point_vec = [x[i], y[i], z[i]]
        angle = angle_between_vectors(point_vec, view_pos)

        #print py.arctan2(y[i],x[i])*180/py.pi
        #raw_input()
        if angle < 95.0:
            front_x.append(x[i])
            front_y.append(y[i])
            front_z.append(z[i])
        else:
            back_x.append(x[i])
            back_y.append(y[i])
            back_z.append(z[i])

    # Plot the points using plot3D
    if ls in ["-", "--", ":", "-."]:
        def plot(x, y, z, delta, alpha):
            delta = float(delta)
            j = 0
            i = 0

            for i in range(1, len(x)):
                mag = py.norm([x[i] - x[i - 1],
                               y[i] - y[i - 1],
                               z[i] - z[i - 1]])
                if mag > delta:
                    axis.plot3D(x[j:i - 1],
                                y[j:i - 1],
                                z[j:i - 1],
                                ls, alpha=alpha, color=color, lw=lw)
                    j = i
            if j == 0:
                print "Bugger"
            if i != j:
                axis.plot3D(x[j:i - 1], y[j:i - 1], z[j:i - 1],
                            ls, alpha=alpha, color=color, lw=lw)

        plot(front_x, front_y, front_z, delta, alpha=1.0)
        plot(back_x, back_y, back_z, delta, alpha=0.1)
    else:
        axis.plot3D(front_x, front_y, front_z,
                    ls, alpha=1.0, color=color, lw=lw)
        axis.plot3D(back_x, back_y, back_z,
                    ls, alpha=0.3, color=color, lw=lw)


def Fit_Function(x, y, n):
    """ Fits a polynomial of degree n to y at the points x. """

    f_p = py.polyfit(x, y, n)

    y_fit = []
    x_fit = py.linspace(x[0], x[-1], 100 * len(x))
    for i in range(len(x_fit)):
        f_val = 0.0
        for j in range(n + 1):
            f_val += f_p[j] * pow(x_fit[i], n - j)
        y_fit.append(f_val)

    return x_fit, y_fit, f_p

