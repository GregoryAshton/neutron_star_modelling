#!/usr/bin/python

import numpy as np
import pylab as py
from math import pi
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d

# Import functions external files
import File_Functions
import Physics_Functions
import Useful_Tools
import NLD_Functions
from Physics_Functions import Beta_Function


#def Save_Figure(file_name, type_of_plot, format_type=".png"):
#    """ Saves the current figure with a relevant file name"""
#    plot_file_name = type_of_plot + "_" + \
#                    file_name.rstrip(".hdf5") + format_type
#    py.savefig(plot_file_name)
#    print "Saving figure as %s" % plot_file_name


def Additional_Code(Option_Dictionary):
    """
    Allows the user to 'append' python code in to Plot.py for example to change
    the xlim of a plot without changing the source code.


    """
    try:
        raw_str_code = Option_Dictionary['raw']
        exec raw_str_code
    except SyntaxError:
        print (" SyntaxError: Something is wrong with how you"
               "provided the raw code, I recived \n {}".format(raw_str_code)
              )
    except KeyError:
        # No additional code
        pass


def Defaults():
    """ Default plotting options """
    # Set the default font for all plots
    from matplotlib import rc
    rc('font', **{'family': 'serif',
        'serif': ['Computer Modern']})
    rc('text', usetex=True)

    # Set the defaults for axis
    py.rcParams['axes.color_cycle'] = ['k', 'b', 'r', 'g']
    py.rcParams['font.size'] = 15
    py.rcParams['axes.labelsize'] = 20
    py.rcParams['lines.linewidth'] = 2
    py.rcParams['axes.grid'] = True
    py.rcParams['figure.figsize'] = (10.0, 8.0)
    #py.subplots_adjust(left=0.13, right=0.9, top=0.9, bottom=0.12,hspace=0.0)


# Plotting functions
def Simple_Plot(file_name, Option_Dictionary={}):
    """

    Plots the given file as a function of time

    """

    # Import the data
    #(time,x,y,z) = File_Functions.Import_Data(file_name,max_n)

    (time, x, y, z) = \
            File_Functions.One_Component_Import(file_name)

    # Handle any additional options which are in the dictionary
    if 'tmax' in Option_Dictionary:
        tmax = Option_Dictionary['tmax']
    else:
        tmax = max(time)
    if 'tmin' in Option_Dictionary:
        tmin = Option_Dictionary['tmin']
    else:
        tmin = min(time)

    fig1 = py.subplot(3, 1, 1)
    fig1.set_xticklabels([])
    fig1.plot(time, x)
    py.ylabel(r"$\omega_{x}$", rotation="horizontal")
    py.yticks(fig1.get_yticks()[1:-1])
    py.xlim(tmin, tmax)

    fig2 = py.subplot(3, 1, 2)
    fig2.set_xticklabels([])
    fig2.plot(time, y)
    py.yticks()
    py.ylabel("$\omega_{y}$", rotation="horizontal")
    py.yticks(fig2.get_yticks()[1:-1])
    py.xlim(tmin, tmax)

    fig2 = py.subplot(3, 1, 3)
    fig2.plot(time, z)
    py.ylabel("$\omega_{z} $", rotation="horizontal")
    py.xlim(tmin, tmax)

    py.xlabel(r"$t$")

    Additional_Code(Option_Dictionary)

    py.show()


def Spherical_Plot(file_name, Option_Dictionary={}):
    """

    Plot the input data after transforming to spherical polar coordinates

    The opts dictionary may contain
    nmax=int limit the data from 0:nmax
    tmax=float, tmin=float ~ limit the xaxis
    end_val=True ~ print the average of the last 100 points
    save_fig=True ~ saves the figure
    """

    # Default settings
    labelx = -0.1  # x position of the yaxis labels

    # Handle any additional options which are in the dictionary

    (time, omega_x, omega_y, omega_z) = \
            File_Functions.One_Component_Import(file_name)

    if 'tmax' in Option_Dictionary:
        tmax = Option_Dictionary['tmax']
    else:
        tmax = max(time)
    if 'tmin' in Option_Dictionary:
        tmin = Option_Dictionary['tmin']
    else:
        tmin = 0.0

    # Transform to spherical polar coordinates
    (omega, a, phi) = Physics_Functions.Cartesian_2_Spherical(
                        omega_x, omega_y, omega_z, fix_phi=True)

    # Function to help scale the x-axis
    (t_scaled, scale_val) = Useful_Tools.Sort_Out_Some_Axis(time)

    # Plot omega(t)
    fig = py.figure()
    ax1 = fig.add_subplot(3, 1, 1)
    ax1.set_xticklabels([])
    ax1.plot(t_scaled, omega)
    ax1.set_xlim(tmin * pow(10, -scale_val), tmax * pow(10, -scale_val))

    ax1.set_ylim(0, 1.1 * max(omega))
    #py.yticks(fig1.get_yticks()[1:-1])
    ax1.set_ylabel(r"$\omega$  [Hz] ", rotation="vertical")
    ax1.yaxis.set_label_coords(labelx, 0.5)

    # Plot a(t)
    ax2 = fig.add_subplot(3, 1, 2)
    ax2.set_xticklabels([])
    ax2.plot(t_scaled, a)
    #py.axhline(90,ls="--",color="k")

    ax2.set_ylim(0, 105)
    #py.yticks(fig2.get_yticks()[0:-2])
    ax2.set_yticks(py.arange(0, 105, 15))
    ax2.set_ylabel("$a $ [deg]", rotation="vertical")
    ax2.yaxis.set_label_coords(labelx, 0.5)
    ax2.set_xlim(tmin * pow(10, -scale_val), tmax * pow(10, -scale_val))

    # Plot phi(t)
    ax3 = fig.add_subplot(3, 1, 3)
    ax3.plot(t_scaled, phi)

    #Ploptions
    #ax3.set_ylim(0,110)
    #ax3.set_yticks(py.arange(0,105,15))
    ax3.set_yticks(ax3.get_yticks()[0:-1])
    ax3.set_ylabel("$\phi$ [deg]", rotation="vertical")
    ax3.yaxis.set_label_coords(labelx, 0.5)
    ax3.set_xlabel(r"time  [$1\times 10^{}$ s]".format(scale_val))
    ax3.set_xlim(tmin * pow(10, -scale_val), tmax * pow(10, -scale_val))
    if 'end_val' in Option_Dictionary:
        print " Data on the end value of the spherical components of omega"
        omega_end = omega[-100:-1]
        print ("Average of |omega|: {0} s^-1 \n Range of omega : {1}"
                .format(py.average(omega_end), max(omega_end) - min(omega_end)))
        a_end = a[-100:-1]
        print (" Average of a: {0} s^-1 \n Range of a : {1}"
              .format(py.average(a_end), max(a_end) - min(a_end)))
        phi_end = omega[-100:-1]
        print (" Average of phi: {0} s^-1 \n Range of phi : {1}"
            .format(py.average(phi_end), max(phi_end) - min(phi_end)))

    py.subplots_adjust(left=0.13, right=0.9, top=0.9, bottom=0.12, hspace=0.0)

    Additional_Code(Option_Dictionary)

    if 'save_fig' in Option_Dictionary and Option_Dictionary['save_fig']:
        File_Functions.Save_Figure(file_name, "Spherical_Plot")
    else:
        py.show()


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
    chi = float(Parameter_Dictionary["chi"]) * pi / 180  # radians

    # Transform to spherical polar coordinates specifying that we want
    # the angles to be in Radians rather than degrees
    (omega, a, phi) = Physics_Functions.Cartesian_2_Spherical(
                        omega_x, omega_y, omega_z, "Radians")

    # Function to help scale the t-axis
    (t_scaled, scale_val) = Useful_Tools.Sort_Out_Some_Axis(time)

    # Calculate the angle made with the magnetic dipole
    # assumed to lie at chi to the z axis in the x-z plane
    def alpha_func(a, phi, chi):
        """ Calculate angle between omega and magnetic dipole """
        chi_radians = chi * pi / 180
        Sx = py.sin(chi_radians)
        Cx = py.cos(chi_radians)
        return py.arccos(Sx * py.sin(a) * py.cos(phi) + Cx * py.cos(a))

    alpha = [alpha_func(a[i], phi[i], chi) * 180 / pi for i in range(len(a))]

    fig = py.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(t_scaled, alpha, lw=2)
    ax1.set_ylabel(r"$\alpha$ [deg]", rotation="vertical")
    ax1.set_xlabel(r"time  [$1\times 10^{}$ s]".format(str(scale_val)))

    Additional_Code(Option_Dictionary)

    if 'save_fig' in Option_Dictionary and Option_Dictionary['save_fig']:
        File_FUnctions.Save_Figure(file_name, "Alpha")
    else:
        py.show()


def ThreeD_Plot_Cartesian(file_name, Option_Dictionary={}):
    """

    Plots the components of input file in 3D Option_Dictionary, takes:
    # Integers to retrict the number of plotted points
    start = int
    stop = int

    # Colouring commands
    power = float
    """

    # Import the data in components x,y,z we use generic names incase the
    # effective body frame (EBF) axis is used

    f = File_Functions.Read_File(file_name)
    time = f['time'].value
    x = f['w1'].value
    y = f['w2'].value
    z = f['w3'].value

    # Set the defaults and then overide if they exist in Option_Dictionary
    start = 0
    stop = -1
    if "start" in Option_Dictionary:
        start = int(Option_Dictionary["start"])
    if "stop" in Option_Dictionary:
        stop = int(Option_Dictionary["stop"])

    # Reduce the number of points
    time = time[start:stop]
    x = x[start:stop]
    y = y[start:stop]
    z = z[start:stop]

    # Check if we should use the Effective body frame axis
    if 'EBF' in Option_Dictionary:
        if Option_Dictionary.get('verbose'):
            print " Using the effective body frame axis in this plot "

        Parameter_Dictionary = File_Functions.Parameter_Dictionary(file_name)
        epsI = float(Parameter_Dictionary["epsI"])
        epsA = float(Parameter_Dictionary["epsA"])
        chi = float(Parameter_Dictionary["chi"]) * pi / 180
        beta = Beta_Function(epsI, epsA, chi)
        (x, y, z) = Physics_Functions.Cartesian_2_EBF(x, y, z, beta)

    # Create subplot and define view angle, this is fixed
    ax = py.subplot(111, projection='3d')
    elev = 15.0
    azim = -134.0
    ax.view_init(elev, azim)

    # Create and label the primed axis
#    ax.plot(py.zeros(100),py.zeros(100),py.linspace(-max(z),max(z),100),
#                color="k")
#    ax.text(0,0,max(z)*1.1,"$z'$")
#    ax.plot(py.zeros(100),py.linspace(max(y),min(y),100),py.zeros(100),
#                color="k")
#    ax.text(0,max(y)*1.1,0,"$y'$")
#    ax.plot(py.linspace(max(x),min(x),100),py.zeros(100),py.zeros(100),
                 #color="k")
#    ax.text(max(x)*1.1,0,0,"$x'$")

    # Compute same variables used for colouring and plot the x',y' and z'
    # transforming the colour as time changes
    if 'power' in Option_Dictionary:
        power = float(Option_Dictionary['power'])
#        n=len(x)
#        d = int(Useful_Tools.Round_To_n(n,0))/100
#        s=n/d
#        for i in range(1,d-1):
#            ax.plot(x[s*i:s*i+s],y[s*i:s*i+s],z[s*i:s*i+s],
#            color=(0.0,1-float(pow(i,power)*1.0/pow(d,power)),0.8),alpha=0.5)

        n = len(time)
        for i in range(0, n - 10, 10):
            color = (0.0,
                     1 - float(pow(i, power) * 1.0 / pow(n, power)),
                     float(pow(i, power) * 1.0 / pow(n, power))
                    )
            ax.plot(x[i:i + 11], y[i:i + 11], z[i:i + 11],
                                color=color, alpha=1.0)
    else:
        ax.plot(x, y, z)

    if 'EBF' in Option_Dictionary:
        ax.set_xlabel(r"$e_{1}$")
        ax.set_ylabel(r"$e_{2}$")
        ax.set_zlabel(r"$e_{3}$")
        # These axis labels may prove to be wrong in the case of epsI<0
    else:
        ax.set_xlabel(r"$\hat{\omega_{x}}$")
        ax.set_ylabel(r"$\hat{\omega_{y}}$")
        ax.set_zlabel(r"$\hat{\omega_{z}}$")
    ax.grid(False)

    # Remove ticks since we do not need them
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_zticklabels([])

    Additional_Code(Option_Dictionary)

    if 'save_fig' in Option_Dictionary and Option_Dictionary['save_fig']:
        File_FUnctions.Save_Figure(file_name, "ThreeD_Plot_Cartesian")
    else:
        py.show()


def Angle_Space_Plot(file_name, Option_Dictionary={}):
    """

    Used to plot the angular components against each other
      of the spin vector. Option_Dictionary takes the following as arguments
    nmax : int ~ take only the first nmax points from the file
    EBF : True ~ Rotate the Effective Body Frame axis
    2D : True ~ This will plot the angular components phi and a in normal plot
    3D : True ~ This will plot the angular components phi and a projected onto
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
        chi = float(Parameter_Dictionary["chi"]) * pi / 180
        beta = Beta_Function(epsI, epsA, chi)
        (x, y, z) = Physics_Functions.Cartesian_2_EBF(x, y, z, beta)

    if '2D' in Option_Dictionary:

        # Transform to spherical coordinates
        (omega, a, phi) = Physics_Functions.Cartesian_2_Spherical(
                                            x, y, z, fix_phi=True)

        fig = py.figure()
        ax1 = fig.add_subplot(111)
        n = len(time)

        ax1.plot(phi, a)
        if 'beta' in Option_Dictionary:
            ax1.set_xlabel("$\phi'$ [deg]")
            ax1.set_ylabel("$a'$ [deg]")
        else:
            ax1.set_xlabel("$\phi$ [deg]")
            ax1.set_ylabel("$a$ [deg]")
        #ax1.set_xlim(0,180)
        #ax1.set_ylim(0,180)

        if 'arrow' in Option_Dictionary and Option_Dictionary['arrow']:
            for i in range(0, n, int(n / 20.0)):
                py.arrow(phi[i], a[i],
                0.1 * (phi[i + 1] - phi[i]), 0.1 * (a[i + 1] - a[i]),
                color="k", fill=True, head_width=1.0,
                head_starts_at_zero=True, alpha=0.8
                )

            else:
                try:
                    i = int(Option_Dictionary['arrow'])
                    py.arrow(phi[i], a[i],
                    0.1 * (phi[i + 1] - phi[i]), 0.1 * (a[i + 1] - a[i]),
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
        (omega, a, phi) = Physics_Functions.Cartesian_2_Spherical(
                                        x, y, z,
                                        fix_phi=True,
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
                x = [py.sin(a[i]) * py.cos(phi[i]) for i in range(low, high)]
                y = [py.sin(a[i]) * py.sin(phi[i]) for i in range(low, high)]
                z = [py.cos(a[i]) for i in range(low, high)]
                Useful_Tools.ThreeD_Sphere(ax, elevation, azimuth, x, y, z,
                                        ls="-", lw=0.8, color=colors[j / 2])

        else:
            # Trajectory of path
            x = [py.sin(a[i]) * py.cos(phi[i]) for i in range(len(time))]
            y = [py.sin(a[i]) * py.sin(phi[i]) for i in range(len(time))]
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

        Additional_Code(Option_Dictionary)

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

    epsI = float(Parameter_Dictionary["epsI"])
    epsA = float(Parameter_Dictionary["epsA"])
    chi = float(Parameter_Dictionary["chi"]) * pi / 180
    beta = Beta_Function(epsI, epsA, chi)

    (w1_prime, w2_prime, w3_prime) = Physics_Functions.Cartesian_2_EBF(
                                                        w1, w2, w3, beta)

    fig1 = py.subplot(3, 1, 1)
    fig1.plot(time, w1_prime)
    py.ylabel(r"$\omega_{x}' $", fontsize=20, rotation="horizontal")

    fig2 = py.subplot(312)
    fig2.plot(time, w2_prime)
    py.yticks()
    py.ylabel("$\omega_{y}' $", fontsize=20, rotation="horizontal")

    fig2 = py.subplot(313)
    fig2.plot(time, w3_prime)
    py.xlabel(r"$t$", fontsize=20)
    py.ylabel("$\omega_{z}' $", fontsize=20, rotation="horizontal")
    py.subplots_adjust(left=0.13, right=0.9, top=0.9, bottom=0.12, hspace=0.0)

    Additional_Code(Option_Dictionary)

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

    epsI = float(Parameter_Dictionary['epsI'])
    epsA = float(Parameter_Dictionary['epsA'])
    chi = float(Parameter_Dictionary['chi']) * pi / 180
    beta = Beta_Function(epsI, epsA, chi)
    print
    print 'Beta  = {}s degrees for %s'.format(beta * 180 / pi,
            file_name)

    (w1_prime, w2_prime, w3_prime) = \
        Physics_Functions.Cartesian_2_EBF(w1, w2, w3, beta)

    # Transform to the spherical polar coordinates

    (omega_prime, a_prime, phi_prime) = \
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
        py.ylabel("$a' \;[^{\circ}]$", rotation='horizontal',
                  fontsize=18)

        # Plot phi_prime(t)

        ax2 = ax1.twinx()

        phi_prime = Physics_Functions.Fix_Phi(phi_prime)
        if abs(phi_prime[-1]) > 100:

            def Scale_Axis(axis):
                """ """
                max_item = max(axis)
                min_item = min(axis)
                if abs(max_item) < abs(min_item):
                    max_item = abs(min_item)
                scale = Useful_Tools.Round_To_n(max_item, 0)
                axis_scaled = [ai / scale for ai in axis]
                return (axis_scaled, scale)

            (phi_prime_scaled, scale) = Scale_Axis(phi_prime)
            ax2.plot(t_scaled, phi_prime_scaled, color='b', lw=1.0)

            py.ylabel(r"$\phi' \; 1\times 10^{"
                      + str(int(py.log10(scale))) + "} [^{\circ}]$",
                      rotation='vertical')
        else:

            ax2.plot(t_scaled, phi_prime, color='b', lw=1.0)
            py.ylabel("$\phi' \; [^{\circ}]$", rotation='vertical')

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
        ax2.set_ylabel("$a'$ [deg]", rotation='vertical')
        ax2.yaxis.set_label_coords(labelx, 0.5)

        # Plot phi_prime(t)

        ax3 = fig.add_subplot(313)

        # Check and fix rotations of 2pi in phi

        phi_prime = Physics_Functions.Fix_Phi(phi_prime)

        # Often phi becomes very large in which case we scale the axis

        if abs(phi_prime[-1]) > 1000:

            (phi_prime_scaled, scale) = \
                Useful_Tools.Sort_Out_Some_Axis(phi_prime)
            ax3.plot(t_scaled, phi_prime_scaled)
            ax3.set_ylabel(r"$\phi' [\;1\times 10^{"
                           + str(int(py.log10(scale))) + '} $deg]',
                           rotation='vertical')
        else:

            ax3.plot(t_scaled, phi_prime)
            ax3.set_ylabel("$\phi'$  [deg]", rotation='vertical')

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
        phi_end = phi_prime[-100:-1]
        print ' Average of phi :  %s  degrees \n Range of phi : %s' \
            % (py.average(phi_end), max(phi_end) - min(phi_end))

    py.subplots_adjust(left=0.13, right=0.9, top=0.9, bottom=0.12,
                       hspace=0.0)

    Additional_Code(Option_Dictionary)

    if 'save_fig' in Option_Dictionary:
        File_Functions.Save_Figure(file_name, 'Spherical_Plot_Transform')
    else:
        py.show()


def Observables_Plot(file_name):
    """

    Plot the some physical observables

    """

    # Default settings
    labelx = -0.1  # x position of the yaxis labels

    (time, w1, w2, w3) = File_Functions.One_Component_Import(file_name)
    (t_scaled, scale_val) = Useful_Tools.Sort_Out_Some_Axis(time)
    N = len(time)

    nu = [2 * pi * py.norm([w1[i], w2[i], w3[i]]) for i in xrange(N)]
    (time_2, omega_dot) = NLD_Functions.Dotted_Variable_Triaxial(file_name)
    nu_dot = [o * 2 * pi for o in omega_dot]

    T_res = Physics_Functions.T_residual(time, w1, w2, w3)

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
    ax3.plot(t_scaled, T_res)
    ax3.set_ylabel(r"$T_{\textrm{res}}$")
    ax1.set_xlabel(r"time  [$1\times 10^{}$ s]".format(str(scale_val)))
    ax3.yaxis.set_label_coords(labelx, 0.5)

    py.subplots_adjust(hspace=0.0)
