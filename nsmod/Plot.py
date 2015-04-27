#!/usr/bin/python

import numpy as np
import pylab as py
from matplotlib import pylab as plt
from math import pi
from matplotlib.patches import FancyArrowPatch
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
#from mpl_toolkits.mplot3d import proj3d

# Import functions external files
import File_Functions
import Physics_Functions
import Useful_Tools
from Physics_Functions import Beta_Function
from nsmod.Pulse_width_fitting import W50

from matplotlib import rc_file, ticker
rc_file("/home/greg/Neutron_star_modelling/matplotlibrc")

SCI_FORMATTER = ticker.ScalarFormatter(useOffset=False,
                                       useMathText=True)

# Plotting functions
def simple_plot(file_name, tmax=None, tmin=None, axes=None, *args, **kwargs):
    """

    Plots `file_name` in cartesian components

    Parameters
    ----------
    file_name: str
    tmax: float 
    tmin: float
    axes : 3 axes instances,
        If given plot on top of these figures, else new ones are created.

    """
    
    (time, wx, wy, wz) = File_Functions.One_Component_Import(file_name)

    if axes:
        (ax1, ax2, ax3) = axes
    else: 
        fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, sharex=True)

    ax1.plot(time, wx, *args, **kwargs)
    ax1.set_ylabel(r"$\omega_{x}$", rotation="horizontal")
    ax1.set_yticks(ax1.get_yticks()[1:-1])
    ax1.set_xlim(tmin, tmax)

    ax2.plot(time, wy, *args, **kwargs)
    ax2.set_ylabel(r"$\omega_{y}$", rotation="horizontal")
    ax2.set_yticks(ax2.get_yticks()[1:-1])
    ax2.set_xlim(tmin, tmax)

    ax3.plot(time, wz, *args, **kwargs)
    ax3.set_ylabel(r"$\omega_{z} $", rotation="horizontal")
    ax3.set_xlim(tmin, tmax)

    ax3.set_xlabel(r"$t$")

    return (ax1, ax2, ax3)

def Spherical_Plot(file_name, axes=None, tmax=None, tmin=0.0, 
                   end_val=False, save_fig=False, **kwargs):
    """

    Plot the input data after transforming to spherical polar coordinates

    Parameters
    ----------
    file_name : string
        .hdf5 file name
    fig : matplotlib figure instance
        If none a new instance is created, if exists it must contain 3 axes
    tmax : float
        Maximum time to plot over
    tmin : float
        Minimum time to plot over
    end_val : bool
        If True print the average of the last 100 points
    save_fig : bool
        If True save an appropriately named figure
    Returns
    -------
    fig : matplotlib figure instance
        
    """

    if axes:
        (ax1, ax2, ax3) = axes
    else:
        fig, (ax1, ax2, ax3) = plt.subplots(nrows=3)

    # Default settings
    labelx = -0.1  # x position of the yaxis labels

    # Handle any additional options which are in the dictionary

    (time, omega_x, omega_y, omega_z) = \
            File_Functions.One_Component_Import(file_name)

    if not tmax:
        tmax = max(time)


    # Transform to spherical polar coordinates
    (omega, a, varphi) = Physics_Functions.Cartesian_2_Spherical(
                        omega_x, omega_y, omega_z, fix_varphi=True)


    # Plot omega(t)
    ax1.set_xticklabels([])
    ax1.plot(time, omega, **kwargs)
    ax1.set_xlim(tmin, tmax)

    ax1.set_ylim(0, 1.1 * max(omega))
    #py.yticks(fig1.get_yticks()[1:-1])
    ax1.set_ylabel(r"$\omega$  [rad/s] ", rotation="vertical")
    ax1.yaxis.set_label_coords(labelx, 0.5)

    # Plot a(t)
    ax2.set_xticklabels([])
    ax2.plot(time, a, **kwargs)
    #py.axhline(90,ls="--",color="k")

    ax2.set_ylim(0, 105)
    #py.yticks(fig2.get_yticks()[0:-2])
    ax2.set_yticks(py.arange(0, 105, 15))
    ax2.set_ylabel(r"$a $ [deg]", rotation="vertical")
    ax2.yaxis.set_label_coords(labelx, 0.5)
    ax2.set_xlim(tmin, tmax)

    # Plot varphi(t)
    ax3.plot(time, varphi, **kwargs)

    #Ploptions
    #ax3.set_ylim(0,110)
    #ax3.set_yticks(py.arange(0,105,15))
    ax3.set_yticks(ax3.get_yticks()[0:-1])
    ax3.set_ylabel(r"$\varphi$ [deg]", rotation="vertical")
    ax3.yaxis.set_label_coords(labelx, 0.5)
    ax3.set_xlabel(r"time [s]")
    #ax3.xaxis.set_major_formatter(SCI_FORMATTER)
    ax3.set_xlim(tmin, tmax)
    if end_val:
        print " Data on the end value of the spherical components of omega"
        omega_end = omega[-100:-1]
        print ("Average of |omega|: {0} s^-1 \n Range of omega : {1}"
                .format(py.average(omega_end), max(omega_end) - min(omega_end)))
        a_end = a[-100:-1]
        print (" Average of a: {0} s^-1 \n Range of a : {1}"
              .format(py.average(a_end), max(a_end) - min(a_end)))
        varphi_end = omega[-100:-1]
        print (" Average of varphi: {0} s^-1 \n Range of varphi : {1}"
            .format(py.average(varphi_end), max(varphi_end) - min(varphi_end)))

    py.subplots_adjust(left=0.13, right=0.9, top=0.9, bottom=0.12, hspace=0.0)
    
    if save_fig:
        File_Functions.Save_Figure(file_name, 'Spherical_Plot')

    return (ax1, ax2, ax3)

def Alpha_Plot(file_name, Option_Dictionary={}):
    """

    Plots the alignment of the input file [t,omega_x , omega_y and omega_z]
    against the magnetic dipole

    """

    # Import the data in components x,y,z
    #(time,x,y,z) = File_Functions.Import_Data(file_name,max_n)
    f = File_Functions.Read_File(file_name)
    time = f['time'].value
    omega_x = f['w1'].value
    omega_y = f['w2'].value
    omega_z = f['w3'].value

    # Get the paramters of the run
    Parameter_Dictionary = File_Functions.Parameter_Dictionary(file_name)

    # Extract some parameters about the pulsar
    chi0 = float(Parameter_Dictionary["chi0"]) * pi / 180  # radians

    # Transform to spherical polar coordinates specifying that we want
    # the angles to be in Radians rather than degrees
    (omega, a, varphi) = Physics_Functions.Cartesian_2_Spherical(
                        omega_x, omega_y, omega_z, "Radians")

    # Function to help scale the t-axis
    (t_scaled, scale_val) = Useful_Tools.Sort_Out_Some_Axis(time)

    # Calculate the angle made with the magnetic dipole
    # assumed to lie at chi0 to the z axis in the x-z plane
    def alpha_func(a, varphi, chi0):
        """ Calculate angle between omega and magnetic dipole """
        chi0_radians = chi0 * pi / 180
        Sx = py.sin(chi0_radians)
        Cx = py.cos(chi0_radians)
        return py.arccos(Sx * py.sin(a) * py.cos(varphi) + Cx * py.cos(a))

    alpha = [alpha_func(a[i], varphi[i], chi0) * 180 / pi for i in range(len(a))]

    fig = py.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(t_scaled, alpha, lw=2)
    ax1.set_ylabel(r"$\alpha$ [deg]", rotation="vertical")
    ax1.set_xlabel(r"time  [$1\times 10^{}$ s]".format(str(scale_val)))

    if 'save_fig' in Option_Dictionary and Option_Dictionary['save_fig']:
        File_Functions.Save_Figure(file_name, "Alpha")
    else:
        py.show()

# HACK DUE TO WEIRD IMPORT ERROR
#def ThreeD_Plot_Cartesian(file_name, Option_Dictionary={}):
#    """
#
#    Plots the components of input file in 3D Option_Dictionary, takes:
#    # Integers to retrict the number of plotted points
#    start = int
#    stop = int
#
#    # Colouring commands
#    power = float
#    """
#
#    # Import the data in components x,y,z we use generic names incase the
#    # effective body frame (EBF) axis is used
#
#    f = File_Functions.Read_File(file_name)
#    time = f['time'].value
#    x = f['w1'].value
#    y = f['w2'].value
#    z = f['w3'].value
#
#    # Set the defaults and then overide if they exist in Option_Dictionary
#    start = 0
#    stop = -1
#    if "start" in Option_Dictionary:
#        start = int(Option_Dictionary["start"])
#    if "stop" in Option_Dictionary:
#        stop = int(Option_Dictionary["stop"])
#
#    # Reduce the number of points
#    time = time[start:stop]
#    x = x[start:stop]
#    y = y[start:stop]
#    z = z[start:stop]
#
#    # Check if we should use the Effective body frame axis
#    if 'EBF' in Option_Dictionary:
#        if Option_Dictionary.get('verbose'):
#            print " Using the effective body frame axis in this plot "
#
#        Parameter_Dictionary = File_Functions.Parameter_Dictionary(file_name)
#        epsI = float(Parameter_Dictionary["epsI"])
#        epsA = float(Parameter_Dictionary["epsA"])
#        chi0 = float(Parameter_Dictionary["chi0"]) * pi / 180
#        beta = Beta_Function(epsI, epsA, chi0)
#        (x, y, z) = Physics_Functions.Cartesian_2_EBF(x, y, z, beta)
#
#    # Create subplot and define view angle, this is fixed
#    ax = py.subplot(111, projection='3d')
#    elev = 15.0
#    azim = -134.0
#    ax.view_init(elev, azim)
#
#    # Create and label the primed axis
##    ax.plot(py.zeros(100),py.zeros(100),py.linspace(-max(z),max(z),100),
##                color="k")
##    ax.text(0,0,max(z)*1.1,"$z'$")
##    ax.plot(py.zeros(100),py.linspace(max(y),min(y),100),py.zeros(100),
##                color="k")
##    ax.text(0,max(y)*1.1,0,"$y'$")
##    ax.plot(py.linspace(max(x),min(x),100),py.zeros(100),py.zeros(100),
#                 #color="k")
##    ax.text(max(x)*1.1,0,0,"$x'$")
#
#    # Compute same variables used for colouring and plot the x',y' and z'
#    # transforming the colour as time changes
#    if 'power' in Option_Dictionary:
#        power = float(Option_Dictionary['power'])
##        n=len(x)
##        d = int(Useful_Tools.Round_To_n(n,0))/100
##        s=n/d
##        for i in range(1,d-1):
##            ax.plot(x[s*i:s*i+s],y[s*i:s*i+s],z[s*i:s*i+s],
##            color=(0.0,1-float(pow(i,power)*1.0/pow(d,power)),0.8),alpha=0.5)
#
#        n = len(time)
#        for i in range(0, n - 10, 10):
#            color = (0.0,
#                     1 - float(pow(i, power) * 1.0 / pow(n, power)),
#                     float(pow(i, power) * 1.0 / pow(n, power))
#                    )
#            ax.plot(x[i:i + 11], y[i:i + 11], z[i:i + 11],
#                                color=color, alpha=1.0)
#    else:
#        ax.plot(x, y, z)
#
#    if 'EBF' in Option_Dictionary:
#        ax.set_xlabel(r"$e_{1}$")
#        ax.set_ylabel(r"$e_{2}$")
#        ax.set_zlabel(r"$e_{3}$")
#        # These axis labels may prove to be wrong in the case of epsI<0
#    else:
#        ax.set_xlabel(r"$\hat{\omega_{x}}$")
#        ax.set_ylabel(r"$\hat{\omega_{y}}$")
#        ax.set_zlabel(r"$\hat{\omega_{z}}$")
#    ax.grid(False)
#
#    # Remove ticks since we do not need them
#    ax.set_xticklabels([])
#    ax.set_yticklabels([])
#    ax.set_zticklabels([])
#
#    if 'save_fig' in Option_Dictionary and Option_Dictionary['save_fig']:
#        File_FUnctions.Save_Figure(file_name, "ThreeD_Plot_Cartesian")
#    else:
#        py.show()


def Angle_Space_Plot(file_name, Option_Dictionary={}):
    """

    Used to plot the angular components against each other
      of the spin vector. Option_Dictionary takes the following as arguments
    nmax : int ~ take only the first nmax points from the file
    EBF : True ~ Rotate the Effective Body Frame axis
    2D : True ~ This will plot the angular components varphi and a in normal plot
    3D : True ~ This will plot the angular components varphi and a projected onto
                the unit sphere
    3D : elevation/azimuthal the viewing angle deliminated by a slash
        split=n1/n2/n3/n4/... For use with 3D this will segment the data
                     into chunks starting and stopping at the integers n1,n2...

    """

    # Import the data in components x,y,z we use generic x,y,z
    # so as to not confused the EFB

    if "nmax" in Option_Dictionary:
        nmax = Option_Dictionary['nmax']
    else:
        nmax = None

    (time, x, y, z) = File_Functions.One_Component_Import(file_name, nmax)

    # Check if we should use the Effective body frame axis
    if 'EBF' in Option_Dictionary:
        if Option_Dictionary.get('verbose'):
            print " Using the effective body frame axis in this plot "

        Parameter_Dictionary = File_Functions.Parameter_Dictionary(file_name)
        epsI = float(Parameter_Dictionary["epsI"])
        epsA = float(Parameter_Dictionary["epsA"])
        chi0 = float(Parameter_Dictionary["chi0"]) * pi / 180
        beta = Beta_Function(epsI, epsA, chi0)
        (x, y, z) = Physics_Functions.Cartesian_2_EBF(x, y, z, beta)

    if '2D' in Option_Dictionary:

        # Transform to spherical coordinates
        (omega, a, varphi) = Physics_Functions.Cartesian_2_Spherical(
                                            x, y, z, fix_varphi=True)

        fig = py.figure()
        ax1 = fig.add_subplot(111)
        n = len(time)

        ax1.plot(varphi, a)
        if 'beta' in Option_Dictionary:
            ax1.set_xlabel(r"$\varphi'$ [deg]")
            ax1.set_ylabel(r"$a'$ [deg]")
        else:
            ax1.set_xlabel(r"$\varphi$ [deg]")
            ax1.set_ylabel(r"$a$ [deg]")
        #ax1.set_xlim(0,180)
        #ax1.set_ylim(0,180)

        if 'arrow' in Option_Dictionary and Option_Dictionary['arrow']:
            for i in range(0, n, int(n / 20.0)):
                py.arrow(varphi[i], a[i],
                0.1 * (varphi[i + 1] - varphi[i]), 0.1 * (a[i + 1] - a[i]),
                color="k", fill=True, head_width=1.0,
                head_starts_at_zero=True, alpha=0.8
                )

            else:
                try:
                    i = int(Option_Dictionary['arrow'])
                    py.arrow(varphi[i], a[i],
                    0.1 * (varphi[i + 1] - varphi[i]), 0.1 * (a[i + 1] - a[i]),
                    color="k", fill=True, head_width=1.0,
                    head_starts_at_zero=True, alpha=0.8
                    )
                except ValueError:
                    pass

        if 'save_fig' in Option_Dictionary and Option_Dictionary['save_fig']:
            File_Functions.Save_Figure(file_name, "Angle_Space_Plot_2D")
        else:
            py.show()

    if '3D' in Option_Dictionary:

        ax = py.subplot(111, projection='3d')

        # Set the viewing position
        if Option_Dictionary.get('azimuth'):
            azimuth = float(Option_Dictionary['azimuth'])
        else:
            azimuth = -134.0
        if Option_Dictionary.get('elevation'):
            elevation = float(Option_Dictionary['elevation'])
        else:
            elevation = 15.0

        ax.view_init(elevation, azimuth)

        # Transform to spherical coordinates,
        (omega, a, varphi) = Physics_Functions.Cartesian_2_Spherical(
                                        x, y, z,
                                        fix_varphi=True,
                                        Angle_Type="Radians")

        if 'split' in Option_Dictionary:
            # Import the values to split by
            values = [int(it) for it in Option_Dictionary['split'].split("/")]

            # Define color library for the splits
            colors = ["b", "r", "k", "c"]
            for j in range(0, len(values), 2):
                # Trajectory of path
                low = values[j]
                high = values[j + 1]
                x = [py.sin(a[i]) * py.cos(varphi[i]) for i in range(low, high)]
                y = [py.sin(a[i]) * py.sin(varphi[i]) for i in range(low, high)]
                z = [py.cos(a[i]) for i in range(low, high)]
                Useful_Tools.ThreeD_Sphere(ax, elevation, azimuth, x, y, z,
                                        ls="-", lw=0.8, color=colors[j / 2])

        else:
            # Trajectory of path
            x = [py.sin(a[i]) * py.cos(varphi[i]) for i in range(len(time))]
            y = [py.sin(a[i]) * py.sin(varphi[i]) for i in range(len(time))]
            z = [py.cos(a[i]) for i in range(len(time))]

            if 'delta' in Option_Dictionary:
                delta = Option_Dictionary['delta']
            else:
                delta = 1.0
            Useful_Tools.ThreeD_Sphere(ax, elevation, azimuth, x, y, z,
                                ls="-", lw=0.8, color="k", delta=delta)

        if 'arrows' in Option_Dictionary:

            class Arrow3D(FancyArrowPatch):
                def __init__(self, xs, ys, zs, *args, **kwargs):
                    FancyArrowPatch.__init__(self, (0, 0), (0, 0),
                                            *args, **kwargs)
                    self._verts3d = xs, ys, zs

                def draw(self, renderer):
                    xs3d, ys3d, zs3d = self._verts3d
                    xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d,
                                                      renderer.M)
                    self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
                    FancyArrowPatch.draw(self, renderer)

            arrow_points = [int(arrow)
                            for arrow in Option_Dictionary['arrows'].split("/")]
            for i in arrow_points:
                print i
                artist = Arrow3D([x[i], x[i + 1]],
                                 [y[i], y[i + 1]], [z[i], z[i + 1]],
                                 mutation_scale=20, lw=1, arrowstyle="-|>",
                                 color="b"
                                )

                ax.add_artist(artist)

        # Sphere of unit radius
        u = np.linspace(0, 2 * np.pi, 100)
        v = np.linspace(0, np.pi, 100)
        x = 1.0 * np.outer(np.cos(u), np.sin(v))
        y = 1.0 * np.outer(np.sin(u), np.sin(v))
        z = 1.0 * np.outer(np.ones(np.size(u)), np.cos(v))

        ax.plot_surface(x, y, z,
                        rstride=4, cstride=9, color='white', alpha=0.9, lw=0.1)

        if 'EBF' in Option_Dictionary:
            ax.set_xlabel(r"$e_{1}$")
            ax.set_ylabel(r"$e_{2}$")
            ax.set_zlabel(r"$e_{3}$", rotation="horizontal")
            # These axis labels may prove to be wrong in the case of epsI<0
        else:
            ax.set_xlabel(r"$\hat{\omega_{x}}$")
            ax.set_ylabel(r"$\hat{\omega_{y}}$")
            ax.set_zlabel(r"$\hat{\omega_{z}}$")

        # Remove the ticks as these are not required
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_zticklabels([])

        ax.grid(False)

        if 'save_fig' in Option_Dictionary and Option_Dictionary['save_fig']:
            File_Functions.Save_Figure(file_name, "Angle_Space_Plot_3D")
        else:
            py.show()


def Simple_Plot_Transform(file_name, Option_Dictionary={}):
    """

    Same as Simple_Plot() except transform to the effective MOI tensor
    principle axis, this has limited functionality as it is generally
    used for checking data is correct and not final plots

    """

    # Import the data in components x,y,z
    f = File_Functions.Read_File(file_name)
    time = f['time'].value
    w1 = f['w1'].value
    w2 = f['w2'].value
    w3 = f['w3'].value

    # Get the paramters of the run
    Parameter_Dictionary = File_Functions.Parameter_Dictionary(file_name)

    epsI3 = float(Parameter_Dictionary["epsI3"])
    epsA = float(Parameter_Dictionary["epsA"])
    chi0 = float(Parameter_Dictionary["chi0"]) * pi / 180
    beta = Beta_Function(epsI3, epsA, chi0)

    (w1_prime, w2_prime, w3_prime) = Physics_Functions.Cartesian_2_EBF(
                                                        w1, w2, w3, beta)

    fig1 = py.subplot(3, 1, 1)
    fig1.plot(time, w1_prime)
    py.ylabel(r"$\omega_{x}' $", fontsize=20, rotation="horizontal")

    fig2 = py.subplot(312)
    fig2.plot(time, w2_prime)
    py.yticks()
    py.ylabel(r"$\omega_{y}' $", fontsize=20, rotation="horizontal")

    fig2 = py.subplot(313)
    fig2.plot(time, w3_prime)
    py.xlabel(r"$t$", fontsize=20)
    py.ylabel(r"$\omega_{z}' $", fontsize=20, rotation="horizontal")
    py.subplots_adjust(left=0.13, right=0.9, top=0.9, bottom=0.12, hspace=0.0)

    py.show()


def Spherical_Plot_Transform(file_name, Option_Dictionary={}):
    """

    Plot the input data after transforming to spherical polar coordinates in
    the primed coordinates
    options should be passed as a dictionary with the following arguments: \n

    opts : max_len ~ A maximum number of points to include
    end_val : True ~ Prints an average of the last 100 points in the plot
    save_fig : True ~ Save the figure in an appropriate way

    """

    # Default settings

    labelx = -0.1  # x position of the yaxis labels

    # Import the data in components x,y,z

    f = File_Functions.Read_File(file_name)
    time = f['time'].value
    w1 = f['w1'].value
    w2 = f['w2'].value
    w3 = f['w3'].value

    # Get the paramters of the run

    Parameter_Dictionary = \
        File_Functions.Parameter_Dictionary(file_name)

    epsI3 = float(Parameter_Dictionary['epsI3'])
    epsA = float(Parameter_Dictionary['epsA'])
    chi0 = float(Parameter_Dictionary['chi0']) * pi / 180
    beta = Beta_Function(epsI3, epsA, chi0)
    print
    print 'Beta  = {}s degrees for %s'.format(beta * 180 / pi,
            file_name)

    (w1_prime, w2_prime, w3_prime) = \
        Physics_Functions.Cartesian_2_EBF(w1, w2, w3, beta)

    # Transform to the spherical polar coordinates

    (omega_prime, a_prime, varphi_prime) = \
        Physics_Functions.Cartesian_2_Spherical(w1_prime, w2_prime,
            w3_prime)

    if 'no_omega' in Option_Dictionary:

        # Function to help scale the x-axis

        (t_scaled, scale_val) = Useful_Tools.Sort_Out_Some_Axis(time)

        fig = py.figure()
        ax1 = fig.add_subplot(111)

        # Plot a_prime(t)

        ax1.plot(t_scaled, a_prime, lw=1.0)

        # py.axhline(90,ls="--",color="k")

        # Ploptions

        y_max = 105
        py.ylim(0, y_max)

        # py.yticks(fig2.get_yticks()[0:-2])

        py.yticks(py.arange(0, y_max, 15))
        py.ylabel(r"$a' \;[^{\circ}]$", rotation='horizontal',
                  fontsize=18)

        # Plot varphi_prime(t)

        ax2 = ax1.twinx()

        varphi_prime = Physics_Functions.Fix_Varphi(varphi_prime)
        if abs(varphi_prime[-1]) > 100:

            def Scale_Axis(axis):
                """ """
                max_item = max(axis)
                min_item = min(axis)
                if abs(max_item) < abs(min_item):
                    max_item = abs(min_item)
                scale = Useful_Tools.Round_To_n(max_item, 0)
                axis_scaled = [ai / scale for ai in axis]
                return (axis_scaled, scale)

            (varphi_prime_scaled, scale) = Scale_Axis(varphi_prime)
            ax2.plot(t_scaled, varphi_prime_scaled, color='b', lw=1.0)

            py.ylabel(r"$\varphi' \; 1\times 10^{"
                      + str(int(py.log10(scale))) + "} [^{\circ}]$",
                      rotation='vertical')
        else:

            ax2.plot(t_scaled, varphi_prime, color='b', lw=1.0)
            py.ylabel(r"$\varphi' \; [^{\circ}]$", rotation='vertical')

        # Ploptions

        ax2.tick_params(axis='y', colors='blue')
        ax2.yaxis.label.set_color('blue')
        py.xlabel(r"time  [$1\times 10^{" + str(scale_val) + '}$ s]',
                  fontsize=16)
    else:

        # Function to help scale the x-axis

        (t_scaled, scale_val) = Useful_Tools.Sort_Out_Some_Axis(time)

        # Plot omega_prime(t)

        fig = py.figure()
        ax1 = fig.add_subplot(311)
        ax1.set_xticklabels([])
        ax1.plot(t_scaled, omega_prime)

        ax1.set_ylim(0, 1.1 * max(omega_prime))

        # py.yticks(fig1.get_yticks()[1:-1])

        ax1.set_ylabel(r"$\omega'$ [Hz]", rotation='vertical')
        ax1.yaxis.set_label_coords(labelx, 0.5)

        # Plot a_prime(t)

        ax2 = fig.add_subplot(312)
        ax2.set_xticklabels([])
        ax2.plot(t_scaled, a_prime)

        # py.axhline(90,ls="--",color="k")

        y_max = 120
        ax2.set_ylim(0, y_max)

        # py.yticks(fig2.get_yticks()[0:-2])

        ax2.set_yticks(py.arange(0, y_max, 15))
        ax2.set_ylabel(r"$a'$ [deg]", rotation='vertical')
        ax2.yaxis.set_label_coords(labelx, 0.5)

        # Plot varphi_prime(t)

        ax3 = fig.add_subplot(313)

        # Check and fix rotations of 2pi in varphi

        varphi_prime = Physics_Functions.Fix_Varphi(varphi_prime)

        # Often varphi becomes very large in which case we scale the axis

        if abs(varphi_prime[-1]) > 1000:

            (varphi_prime_scaled, scale) = \
                Useful_Tools.Sort_Out_Some_Axis(varphi_prime)
            ax3.plot(t_scaled, varphi_prime_scaled)
            ax3.set_ylabel(r"$\varphi' [\;1\times 10^{"
                           + str(int(py.log10(scale))) + '} $deg]',
                           rotation='vertical')
        else:

            ax3.plot(t_scaled, varphi_prime)
            ax3.set_ylabel(r"$\varphi'$  [deg]", rotation='vertical')

        # Ploptions

        ax3.set_xlabel(r"time  [$1\times 10^{}$ s]".format(str(scale_val)))
        ax3.yaxis.set_label_coords(labelx, 0.5)
        ax3.set_yticks(ax3.get_yticks()[0:-1])

    if 'end_val' in Option_Dictionary:
        print ' Data on the end value of the spherical components of omega'
        omega_end = omega_prime[-100:-1]
        print ' Average of |omega| :  %s s^-1  \n Range of omega : %s' \
            % (py.average(omega_end), max(omega_end) - min(omega_end))
        a_end = a_prime[-100:-1]
        print ' Average of a  : %s degrees \n Range of a : %s' \
            % (py.average(a_end), max(a_end) - min(a_end))
        varphi_end = varphi_prime[-100:-1]
        print ' Average of varphi :  %s  degrees \n Range of varphi : %s' \
            % (py.average(varphi_end), max(varphi_end) - min(varphi_end))

    py.subplots_adjust(left=0.13, right=0.9, top=0.9, bottom=0.12,
                       hspace=0.0)

    if 'save_fig' in Option_Dictionary:
        File_Functions.Save_Figure(file_name, 'Spherical_Plot_Transform')
    else:
        py.show()


def Observables_Plot(file_name):
    """

    Plot some physical observables

    """
    
    import NLD_Functions
    # Default settings
    labelx = -0.1  # x position of the yaxis labels

    (time, w1, w2, w3) = File_Functions.One_Component_Import(file_name)
    (t_scaled, scale_val) = Useful_Tools.Sort_Out_Some_Axis(time)
    N = len(time)

    nu = [2 * pi * py.norm([w1[i], w2[i], w3[i]]) for i in xrange(N)]
    (time_2, omega_dot) = NLD_Functions.Dotted_Variable_Triaxial(file_name)
    nu_dot = [o * 2 * pi for o in omega_dot]

    PhaseResidual = Physics_Functions.PhaseResidual(time, w1, w2, w3)

    # Plotting
    ax1 = py.subplot(311)
    ax1.plot(t_scaled, nu)
    ax1.set_ylabel(r"$\nu$")
    ax1.set_xticklabels([])
    ax1.yaxis.set_label_coords(labelx, 0.5)

    ax2 = py.subplot(312)
    ax2.plot(t_scaled, nu_dot)
    ax2.set_ylabel(r"$\dot{\nu}$")
    ax2.set_xticklabels([])
    ax2.yaxis.set_label_coords(labelx, 0.5)

    ax3 = py.subplot(313)
    ax3.plot(t_scaled, PhaseResidual)
    ax3.set_ylabel(r"$\Delta\Phi$")
    ax3.set_xlabel(r"time  [$1\times 10^{}$ s]".format(str(scale_val)))
    ax3.yaxis.set_label_coords(labelx, 0.5)

    py.subplots_adjust(hspace=0.0)
    py.show()


def Euler_Angles(file_name, axes=None, save_fig=False, analytic=False, 
                 *args, **kwargs):
    """ Plot the Euler angles in three subplots """

    if axes != None:
        (ax1, ax2, ax3)= axes
    else:
        fig, (ax1, ax2, ax3) = plt.subplots(nrows=3)

    labelx = -0.1  # x position of the yaxis labels

    (time, w1, w2, w3) = File_Functions.One_Component_Import(file_name)
    data_dictionary = File_Functions.Read_File(file_name)
    theta = np.degrees(data_dictionary['theta'].value)
    phi = np.degrees(data_dictionary['phi'].value)
    psi = np.degrees(data_dictionary['psi'].value)

    PD = File_Functions.Parameter_Dictionary(file_name)
    #epsI1 = float(PD['epsI1'])
    epsI3 = float(PD['epsI3'])

    #theta, phi, psi = Physics_Functions.Inertial_Frame(
    #                    time, np.array([w1, w2, w3]), epsI3)

    # Plot the analytic solution from Jones 2001
    if analytic:
        parameter_dictionary = File_Functions.Parameter_Dictionary(file_name)
        a0 = parameter_dictionary['a0']
        omega0 = parameter_dictionary['omega0']
        epsI3 = parameter_dictionary['epsI3']

        COLOR = "r"
        lw = 3

        thetadot = 0
        theta = time * thetadot + a0
        ax1.plot(time, theta, ls="--", color=COLOR, lw=lw)

        phidot = omega0
        phi = np.degrees(phidot * time)
        ax2.plot(time, phi, ls="--", color=COLOR, lw=lw)

        psidot = -1 * epsI3 * phidot
        psi = np.degrees(psidot * time) + 90.0 
        ax3.plot(time, psi, ls="--", color=COLOR, lw=lw)

    # Plotting
    ax1.plot(time, theta, *args, **kwargs)
    ax1.set_ylabel(r"$\theta$ [deg]")
    ax1.set_xticklabels([])
    ax1.yaxis.set_label_coords(labelx, 0.5)
    ax1.yaxis.set_major_locator(MultipleLocator(0.5))
    #ax1.set_yticks(ax1.get_yticks()[1:])
    #ax1.ticklabel_format(useOffset=False, axis='y')

    if abs(phi[-1])>10000:
        ax2.plot(time, phi, *args, **kwargs)
        ax2.set_ylabel(r"$\phi$ [deg]")
    else:
        ax2.plot(time, phi, *args, **kwargs)
        ax2.set_ylabel(r"$\phi$", rotation="horizontal")
    ax2.set_xticklabels([])
    ax2.yaxis.set_label_coords(labelx, 0.5)
    ax2.set_yticks(ax2.get_yticks()[1:])


    if abs(psi[-1])>10000:
        ax3.plot(time, psi, *args, **kwargs)
        ax3.set_ylabel(r"$\psi$ [deg]")
    else:
        ax3.plot(time, psi, *args, **kwargs)
        ax3.set_ylabel(r"$\psi$ [deg]")
    ax3.set_xlabel(r"time  [s]")
    ax3.yaxis.set_label_coords(labelx, 0.5)

    py.subplots_adjust(hspace=0.0)

    if save_fig:
        File_Functions.Save_Figure(file_name, 'Euler_Angles')


    return (ax1, ax2, ax3)

def big_phi_dot(file_name, ax=None, save_fig=False, *args, **kwargs):
    """ 
    
    Plot Phi_dot for data given in file_name 
    
    :param file_name: Name of the `hdf5` file containing the data
    :type file_name: str
    :param ax: axis instance to use in plotting, if none a new one will be created

    :save_fig: if True the file will be saved appropriately
    :type save_fig: bool

    """

    if ax is None:
        ax = py.subplot(111)

    time, w1, w2, w3, theta, phi, psi = File_Functions.Euler_Angles_Import(file_name)
    chi0 = float(File_Functions.Parameter_Dictionary(file_name)['chi0'])
    (t_scaled, scale_val) = Useful_Tools.Sort_Out_Some_Axis(time)
    
    Phi_dot_list = Physics_Functions.Phi_dot(np.array([w1, w2, w3]), theta, phi, psi, np.radians(chi0))
    ax.plot(t_scaled, Phi_dot_list, *args, **kwargs) 
    ax.set_xlabel(r"time  [$1\times 10^{}$ s]".format(str(scale_val)))
    ax.set_ylabel(r"$\dot{\Phi}$", rotation="horizontal", size=26)

    return ax


def big_theta(file_name, ax=None, save_fig=False, *args, **kwargs):
    """ 
    
    Plot Phi_dot for data given in file_name 
    
    :param file_name: Name of the `hdf5` file containing the data
    :type file_name: str
    :param ax: axis instance to use in plotting, if none a new one will be created

    :save_fig: if True the file will be saved appropriately
    :type save_fig: bool

    """

    if ax is None:
        ax = py.subplot(111)

    time, w1, w2, w3, theta, phi, psi = File_Functions.Euler_Angles_Import(file_name)
    chi0 = float(File_Functions.Parameter_Dictionary(file_name)['chi0'])
    (t_scaled, scale_val) = Useful_Tools.Sort_Out_Some_Axis(time)
    
    Theta_list = np.degrees(Physics_Functions.Theta(theta, psi, np.radians(chi0)))
    ax.plot(t_scaled, Theta_list, *args, **kwargs) 
    ax.set_xlabel(r"time  [$1\times 10^{}$ s]".format(str(scale_val)))
    ax.set_ylabel(r"$\Theta$", rotation="horizontal", size=26)

    return ax


def PhaseResidual(file_name, ax=None, save_fig=False, order=3, analytic="", 
                    tstart=None, tend=None, *args, **kwargs):
    """ 
    Plot the phase residuals for the data given in file_name. 
    
    Parameters:
    -----------
    file_name : 
        String referencing the h5py data file
    order : 
        Order of polynomial to fit, must be either 2 or 3
    ax: 
        An axis instance to plot, if None a new instance is initiated
    save_fig : 
        Option to save the figure using the default save feature
    analytic : one of, or a list of, "49", "63", "75"
        If given plot the analytic predictions of Jones 2001
    tstart, tend : float
        If given as floats the start and ends points to use
    
    Note: *args and **kwargs are passed onto the matplotlib plot function

    Returns:
    --------
    ax: the axis instance

    """

    out = File_Functions.Euler_Angles_Import(file_name)
    [time, w1, w2, w3, theta, phi, psi] = out
    PD = File_Functions.Parameter_Dictionary(file_name)
    chi0 = np.radians(PD['chi0'])

    if tstart and tend:
        idxs = (tstart < time) * (time < tend)
    elif tstart:
        idxs = time > tstart
    elif tend:
        idxs = time < tend     
    
    try:
        time = time[idxs]
        w1 = w1[idxs]
        w2 = w2[idxs]
        w3 = w3[idxs]
        theta = theta[idxs]
        phi = phi[idxs]
        psi = psi[idxs]
    except NameError:
        pass

    Pres = Physics_Functions.PhaseResidual(time, w1, w2, w3, 
                                             theta, phi, psi, chi0, 
                                             order=order)
    Cycles = Pres / (2*np.pi)
    if not ax:
        fig, ax = plt.subplots()

    ax.plot(time, Cycles, *args, **kwargs)
    ax.set_ylabel(r"Phase residual [cycles]", rotation='vertical')
    ax.set_xlabel(r"time  [s]")
    ax.axhline(0, ls="-", color="k", zorder=-100)

    if "49" in analytic:
        DeltaPhi_49 = PD['DeltaPhi_49']
        ax.axhline(DeltaPhi_49/(2*np.pi), ls="--", color="k", zorder=-100, 
                   label="$|\Delta\Phi^{49}|$")
        ax.axhline(-DeltaPhi_49/(2*np.pi), ls="--", color="k", zorder=-100)
    if "63" in analytic:
        DeltaPhi_63 = PD['DeltaPhi_63']
        ax.axhline(DeltaPhi_63/(2*np.pi), ls="--", color="r", zorder=-100, 
                   label="$|\Delta\Phi^{63}|$")
        ax.axhline(-DeltaPhi_63/(2*np.pi), ls="--", color="r", zorder=-100)
    if "75" in analytic:
        DeltaPhi_75 = PD['DeltaPhi_75']
        ax.axhline(DeltaPhi_75/(2*np.pi), ls="--", color="b", zorder=-100,
                   label="$|\Delta\Phi^{75}|$")
        ax.axhline(-DeltaPhi_75/(2*np.pi), ls="--", color="b", zorder=-100)
    if "49SD" in analytic:
        DeltaPhi_49_SpindownTorque = PD['DeltaPhi_49_SpindownTorque']
        ax.axhline(DeltaPhi_49_SpindownTorque/(2*np.pi), ls="--", color="k", 
                   zorder=-100, label="$|\Delta\Phi^{49}|$")
        ax.axhline(-DeltaPhi_49_SpindownTorque/(2*np.pi), ls="--", color="k",
                   zorder=-100) 
    return ax

def SpindownRate(file_name, ax=None, normalise=False, divisor=10, analytic="",
                 nmax=None, method='Lyne', *args, **kwargs):
    """

    Plot the spindown rate from numeric simulation

    Parameters:
    ----------
    file_name : str
        Name of the h5py data file
    ax : axis instance
        If None a new instance is initiated
    normalise : bool
        Option to normalise the plotted output, useful when comparing several
        different simulations
    method : str
        Either One of "Lyne" or "Numeric" describing how to compute nu_dot

    Note: One can also pass *args and **kwargs onto the matplotlib plot
          function

    Returns:
    --------
    ax : the axis instance

    """

    out_EA = File_Functions.Euler_Angles_Import(file_name, nmax=nmax)
    [time, w1, w2, w3, theta, phi, psi] = out_EA


    PD = File_Functions.Parameter_Dictionary(file_name)
    chi0 = np.radians(PD['chi0'])
    tauP = PD['tauP']

    if method in ['Lyne']:
        time, nu_dot = Physics_Functions.nu_dot_Lyne(time, w1, w2, w3,
                                                     theta, phi, psi, chi0,
                                                     tauP, divisor=divisor)
    elif method in ['Numeric', 'numeric']:
        time, nu_dot = Physics_Functions.nu_dot_numeric(time, theta, psi,
                                                        phi, chi0)

    if not ax:
        fig, ax = plt.subplots()

    if normalise:
        nu_dot = nu_dot / sum(nu_dot**2)**0.5
    else:
        nu_dot = nu_dot

    ax.plot(time, nu_dot, *args, **kwargs)

    ax.set_xlabel(r"time [s]")
    ax.set_ylabel(r"$\dot{\nu}$", rotation="horizontal", size=26)
    #ax.set_ylim(ax.get_ylim()[0], max(ax.get_ylim()[1], 0))
    #ax.axhline(0, ls="--", color="k", zorder=-100)

    # Plot analytic calcualtions
    if PD.has_key('nu_dot0'):
        nu_dot0 = PD['nu_dot0']
    else:
        nu_dot0 = 0

    if "EM" in analytic:
        ax.axhline(nu_dot0, label=r"$\dot{\nu}_{0}$",
                   color="r")

    if "FP" in analytic:
        Delta_nu0 = PD['delta_omega_dot0_FP'] / (2*np.pi)
        ax.axhline(nu_dot0 + Delta_nu0,
                   label=r"$|\Delta\dot{\nu}|_{\mathrm{p}}$")
        ax.axhline(nu_dot0 - Delta_nu0)
        #ax.fill_between(time, nu_dot0-Delta_nu0,
        #                 nu_dot0 + Delta_nu0, color="b",
        #                alpha=0.1, zorder=-100)

    if "58" in analytic:
        c = "green"
        Delta_nu0_58 = PD['delta_omega_dot0_FP_EM'] / (2*np.pi)
        ax.axhline(nu_dot0 + Delta_nu0_58, color=c,
                   label=r"$|\Delta\dot{\nu}|^{58}_{\mathrm{p}}$")
        ax.axhline(nu_dot0 - Delta_nu0_58, color=c)
        #ax.fill_between(time,
        #                 nu_dot0-Delta_nu0_58, nu_dot0 + Delta_nu0_58,
        #                 color=c, alpha=0.2)

    ax.set_xlim(time[0], time[-1])

    return ax


def Intensity(file_name, PhiO, ThetaO, sigmaB,
              ax=None, *args, **kwargs):
    """
    Plot the amplitude using a 2D Gaussian beam model

    Parameters:
    ----------
    Phi0, Theta0 : float
        Observers angular position in the inertial frame
    sigmaPhi, sigmaTheta : float
        Pulse shape parameters
   file_name: string
        Reference to the h5py data file
    ax : matplotlib axis instance, optional
        Axis to plot on
    save_fig : bool, optional
        Save the figure using the default save feature

    Note: One can also pass *args and **kwargs to the matplotlib plot function

    Returns:
    --------
    ax: the axis instance


    """
    if not ax:
        fig, ax = plt.subplots()
    out_EA = File_Functions.Euler_Angles_Import(file_name)
    [time, w1, w2, w3, theta, phi, psi] = out_EA

    PD = File_Functions.Parameter_Dictionary(file_name)
    chi0 = np.radians(PD['chi0'])

    Phi = Physics_Functions.Phi(theta, phi, psi, chi0, fix=True)
    Theta = Physics_Functions.Theta(theta, psi, chi0)

    Intensity = Physics_Functions.Intensity(Phi, Theta, PhiO, ThetaO,
                                            sigmaB, I0=1)
    ax.plot(time, Intensity, *args, **kwargs)

    IMax = Physics_Functions.IntensityMax(Theta, ThetaO, sigmaB, I0=1)
    ax.plot(time, IMax, "--b")

    ax.set_xlabel("time [s]")
    ax.set_ylabel("Normalised Intensity")

    return ax

def PulseWidth(file_name, Theta0, sigmaB, p=50,
               eta=0.01, ax=None, *args, **kwargs):
    """ 
    Plot the pulse width using a 2D Gaussian beam model

    Parameters:
    ----------
    Theta0 : float
        Observers angular position in the inertial frame in radians
    sigmaPhi, sigmaTheta : float
        Pulse shape parameters
    eta : float [0, 1]
        Scaling of the free precession time scale to one in which the modulation
        is negligible
    file_name: string 
        Reference to the h5py data file
    ax : matplotlib axis instance, optional
        Axis to plot on
    save_fig : bool, optional
        Save the figure using the default save feature
    
    Note: One can also pass *args and **kwargs onto the matplotlib plot function

    Returns:
    --------
    ax: the axis instance


    """    
    
    if not ax:
        fig, ax = plt.subplots()

    out_EA = File_Functions.Euler_Angles_Import(file_name)
    [time, w1, w2, w3, theta, phi, psi] = out_EA
    
    PD = File_Functions.Parameter_Dictionary(file_name)
    chi0 = np.radians(PD['chi0'])
    
    #Phi = Physics_Functions.Phi(theta, phi, psi, chi0, fix=True)
    Theta = Physics_Functions.Theta(theta, psi, chi0)
    Phi_dot = Physics_Functions.Phi_dot(np.array([w1, w2, w3]), 
                                             theta, phi, psi, chi0)

    #Amplitude = Physics_Functions.Amplitude(Phi, Theta, Phi0, Theta0, 
    #                              sigmaTheta, sigmaPhi, A0=1) 

    #tCONS = eta * PD['tauP']

    #time_list, W50_list = W50(time, Amplitude, tCONS)
    Wp_list = Physics_Functions.Wp(Phi_dot, Theta, Theta0, sigmaB, p)


    ax.plot(time, Wp_list, *args, **kwargs)

    ax.set_xlabel('time')
    ax.set_ylabel('$W_{p}$', rotation='vertical')

    return ax
