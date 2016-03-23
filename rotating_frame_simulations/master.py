import sys
import os
import matplotlib.pyplot as py
import shutil
from nsmod import Plot, Model, NLD_Functions
from nsmod.File_Functions import Parameter_Dictionary
from nsmod.Useful_Tools import Texify_Float


def get_file_names():
    file_names = ["data/"+file_name for file_name in os.listdir("data") if
                  file_name.endswith(".hdf5")]
    return file_names


def get_low_high():
    file_names = get_file_names()
    for f in file_names:
        if low in f:
            file_low = f
        if high in f:
            file_high = f
    return file_low, file_high

home = str(os.getcwd())
# No anomalous torque

dr_A_NA = "Pulsar_A_No_Anomalous_Torque"
dr_B_NA = "Pulsar_B_No_Anomalous_Torque"
dr_C_NA = "Pulsar_C_No_Anomalous_Torque"
dr_A = "Pulsar_A"
dr_B = "Pulsar_B"
dr_C = "Pulsar_C"

dr_S_NA = "Pulsar_Spherical_No_Anomalous_Torque"

# Imput the known data into dictionary
pulsar_a = {'epsI3': 1.0e-9,
            'epsI1': 0,
            'chi0': 30,
            'a0': 50,
            'epsA': 5.0e-11,
            'omega0': 1.0e4,
            'n': None,
            'T': 1e10,
            'error': 1e-14}
pulsar_b = {'epsI3': 4.0e-11,
            'epsI1': 0,
            'chi0': 0,
            'a0': 50,
            'epsA': 5.0e-11,
            'omega0': 1.0e4,
            'T': 2.0e8,
            'error': 1e-14}
pulsar_c = {'epsI3': 1.0e-15,
            'epsI1': 0,
            'chi0': 0,
            'a0': 50,
            'epsA': 5.0e-11,
            'omega0': 1.0e4,
            'T': 1e8}

pulsar_spherical = {'epsI3': 0.0,
                    'epsI1': 0,
                    'chi0': 30,
                    'a0': 50,
                    'epsA': 5.0e-11,
                    'omega0': 1.0e4,
                    'T': 1e8}

chi0_list = [30.0, 75.0]
low = "3.0"
high = "7.5"

if "table" in sys.argv:
    # Use a function to add some more data to the standard dictionaries
    pulsar_a = Parameter_Dictionary(pulsar_a)
    pulsar_b = Parameter_Dictionary(pulsar_b)
    pulsar_c = Parameter_Dictionary(pulsar_c)

    table1 = r"""
\begin{table}[h]
\centering
\begin{tabular}[h]{|l|c|c|c|c|c|}\hline
&$\epsilon_{I}$ & $\epsilon_{A} $ & $B_{S} $[Gauss] & $\tau_{S}$ [s] & $ \tau_{P}$  [s] \\ \hline
Pulsar A & $ KEY_A_eI $ & $ KEY_A_eA $ & $ KEY_A_Bs $ & $ KEY_A_tauS $ & $ KEY_A_tauP $\\
Pulsar B & $ KEY_B_eI $ & $ KEY_B_eA $ & $ KEY_B_Bs $ & $ KEY_B_tauS $ & $ KEY_B_tauP $\\
Pulsar C & $ KEY_C_eI $ & $ KEY_C_eA $ & $ KEY_C_Bs $ & $ KEY_C_tauS $ & $ KEY_C_tauP $\\ \hline
\end{tabular}
\caption{Table of relevant values for the selected points. In calculating the Magnetic fields we have assumed canonical value of $R=1\times10^{6}$cm for the radius, and a value $\omega_{0} = 1\times10^{4}$ rads $s^{-1}$}
\label{A B C params}
\end{table}

"""

    table1 = table1.replace("KEY_A_eI", Texify_Float(pulsar_a['epsI3']))
    table1 = table1.replace("KEY_B_eI", Texify_Float(pulsar_b['epsI3']))
    table1 = table1.replace("KEY_C_eI", Texify_Float(pulsar_c['epsI3']))

    table1 = table1.replace("KEY_A_eA", Texify_Float(pulsar_a['epsA']))
    table1 = table1.replace("KEY_B_eA", Texify_Float(pulsar_b['epsA']))
    table1 = table1.replace("KEY_C_eA", Texify_Float(pulsar_c['epsA']))

    table1 = table1.replace("KEY_A_Bs", Texify_Float(pulsar_a['Bs']))
    table1 = table1.replace("KEY_B_Bs", Texify_Float(pulsar_b['Bs']))
    table1 = table1.replace("KEY_C_Bs", Texify_Float(pulsar_c['Bs']))

    table1 = table1.replace("KEY_A_tauS", Texify_Float(pulsar_a['tauS']))
    table1 = table1.replace("KEY_B_tauS", Texify_Float(pulsar_b['tauS']))
    table1 = table1.replace("KEY_C_tauS", Texify_Float(pulsar_c['tauS']))

    table1 = table1.replace("KEY_A_tauP", Texify_Float(pulsar_a['tauP']))
    table1 = table1.replace("KEY_B_tauP", Texify_Float(pulsar_b['tauP']))
    table1 = table1.replace("KEY_C_tauP", Texify_Float(pulsar_c['tauP']))

    print table1

    table2 = r"""
\begin{table}[h]
\centering
\begin{tabular}[h]{|l|c|c|c|c|c|c|c|c|}\hline
$\epsilon_{I}$  & $\epsilon_{A} $ & $B_{S} $[Gauss] & $\tau_{S}$ [s] & $\tau_{A}$ [s] & $ \tau_{P}$  [s] & $\beta(\chi0=30^{\circ})$ & $\beta(\chi0=75^{\circ})$ \\ \hline
Y_A_eI $ & $ KEY_A_eA $ & $ KEY_A_Bs $ & $ KEY_A_tauS $ & $ KEY_A_tauA $ & $ KEY_A_tauP $ & $ KEY_A_beta30 $& $ KEY_A_beta75 $ \\
Y_B_eI $ & $ KEY_B_eA $ & $ KEY_B_Bs $ & $ KEY_B_tauS $ & $ KEY_B_tauA $ & $ KEY_B_tauP $ & $ KEY_B_beta30 $& $ KEY_B_beta75 $ \\
Y_C_eI $ & $ KEY_C_eA $ & $ KEY_C_Bs $ & $ KEY_C_tauS $ & $ KEY_C_tauA $ & $ KEY_C_tauP $ & $ KEY_C_beta30 $& $ KEY_C_beta75 $ \\ \hline
\end{tabular}
\caption{Table of relevant values for the selected points. In calculating the Magnetic fields we have assumeda canonical value of $R=1\times10^{6}$cm for the radius, and a value $\omega_{0} = 1\times10^{4}$ rads $s^{-1}$}
\label{A B C params}
\end{table}

"""

    table2 = table2.replace("KEY_A_eI", Texify_Float(pulsar_a['epsI3']))
    table2 = table2.replace("KEY_B_eI", Texify_Float(pulsar_b['epsI3']))
    table2 = table2.replace("KEY_C_eI", Texify_Float(pulsar_c['epsI3']))

    table2 = table2.replace("KEY_A_eA", Texify_Float(pulsar_a['epsA']))
    table2 = table2.replace("KEY_B_eA", Texify_Float(pulsar_b['epsA']))
    table2 = table2.replace("KEY_C_eA", Texify_Float(pulsar_c['epsA']))

    table2 = table2.replace("KEY_A_Bs", Texify_Float(pulsar_a['Bs']))
    table2 = table2.replace("KEY_B_Bs", Texify_Float(pulsar_b['Bs']))
    table2 = table2.replace("KEY_C_Bs", Texify_Float(pulsar_c['Bs']))

    table2 = table2.replace("KEY_A_tauS", Texify_Float(pulsar_a['tauS']))
    table2 = table2.replace("KEY_B_tauS", Texify_Float(pulsar_b['tauS']))
    table2 = table2.replace("KEY_C_tauS", Texify_Float(pulsar_c['tauS']))

    table2 = table2.replace("KEY_A_tauA", Texify_Float(pulsar_a['tauA']))
    table2 = table2.replace("KEY_B_tauA", Texify_Float(pulsar_b['tauA']))
    table2 = table2.replace("KEY_C_tauA", Texify_Float(pulsar_c['tauA']))

    table2 = table2.replace("KEY_A_tauP", Texify_Float(pulsar_a['tauP']))
    table2 = table2.replace("KEY_B_tauP", Texify_Float(pulsar_b['tauP']))
    table2 = table2.replace("KEY_C_tauP", Texify_Float(pulsar_c['tauP']))

    table2 = table2.replace("KEY_A_beta30", Texify_Float(pulsar_a['beta30'],
                            power=False))
    table2 = table2.replace("KEY_B_beta30", Texify_Float(pulsar_b['beta30'],
                            power=False))
    table2 = table2.replace("KEY_C_beta30", Texify_Float(pulsar_c['beta30'],
                            power=False))

    table2 = table2.replace("KEY_A_beta75", Texify_Float(pulsar_a['beta75'],
                            power=False))
    table2 = table2.replace("KEY_B_beta75", Texify_Float(pulsar_b['beta75'],
                            power=False))
    table2 = table2.replace("KEY_C_beta75", Texify_Float(pulsar_c['beta75'],
                            power=False))

    print table2


if "data" in sys.argv:
    if "spherical" in sys.argv:
        pulsar_spherical['AnomTorque'] = False
        if dr_S_NA not in os.listdir("."):
            os.mkdir(dr_S_NA)
        os.chdir(dr_S_NA)
        try:
            shutil.rmtree("data")
        except OSError as e:
            pass

        pulsar_spherical['chi0'] = 30.0
        Model.Run(**pulsar_spherical)

        os.chdir(home)

    if "A_NA" in sys.argv:
        pulsar_a['AnomTorque'] = False
        if dr_A_NA not in os.listdir("."):
            os.mkdir(dr_A_NA)
        os.chdir(dr_A_NA)
        try:
            shutil.rmtree("data")
        except OSError:
            pass

        for chi0 in chi0_list:
            pulsar_a['chi0'] = chi0
            Model.Run(**pulsar_a)

        os.chdir(home)

    if "B_NA" in sys.argv:
        pulsar_b['AnomTorque'] = False

        if dr_B_NA not in os.listdir("."):
            os.mkdir(dr_B_NA)
        os.chdir(dr_B_NA)
        try:
            shutil.rmtree("data")
        except OSError:
            pass

        for chi0 in chi0_list:
                pulsar_b['chi0'] = chi0
                Model.Run(**pulsar_b)

        os.chdir(home)

    if "C_NA" in sys.argv:
        pulsar_c['AnomTorque'] = False

        if dr_C_NA not in os.listdir("."):
            os.mkdir(dr_C_NA)
        os.chdir(dr_C_NA)
        try:
            shutil.rmtree("data")
        except OSError:
            pass

        for chi0 in chi0_list:
                pulsar_c['chi0'] = chi0
                Model.Run(**pulsar_c)

        os.chdir(home)

    if "A" in sys.argv:
        pulsar_a['AnomTorque'] = True
        if dr_A not in os.listdir("."):
            os.mkdir(dr_A)
        os.chdir(dr_A)
        try:
            shutil.rmtree("data")
        except OSError:
            pass

        for chi0 in chi0_list:
            pulsar_a['chi0'] = chi0
            Model.Run(**pulsar_a)

        os.chdir(home)

    if "B" in sys.argv:
        pulsar_b['AnomTorque'] = True
        if dr_B not in os.listdir("."):
            os.mkdir(dr_B)
        os.chdir(dr_B)
        try:
            shutil.rmtree("data")
        except OSError:
            pass

        for chi0 in chi0_list:
            pulsar_b['chi0'] = chi0
            Model.Run(**pulsar_b)

        os.chdir(home)

    if "C" in sys.argv:
        pulsar_c['AnomTorque'] = True
        if dr_C not in os.listdir("."):
            os.mkdir(dr_C)
        os.chdir(dr_C)
        try:
            shutil.rmtree("data")
        except OSError:
            pass

        for chi0 in chi0_list:
            pulsar_c['chi0'] = chi0
            Model.Run(**pulsar_c)

        os.chdir(home)

if "plot" in sys.argv:
    if "spherical" in sys.argv:
        os.chdir(home)
        os.chdir(dr_S_NA)

        file_names = get_file_names()

        Plot.Spherical_Plot(file_names[0], save_fig=True)

    if "A_NA" in sys.argv:
        os.chdir(home)
        os.chdir(dr_A_NA)

        file_low, file_high = get_low_high()

        Plot.Spherical_Plot(file_low, tmax=1e9, save_fig=True)
        Plot.Spherical_Plot(file_high, tmax=1e10, save_fig=True)
        Plot.Angle_Space_Plot(file_low, {'save_fig': False, 'nmax': 1000,
                                         "2D": True, 'delta': 0.4,
                                         'arrows': "20/25/30/195"})

    if "B_NA" in sys.argv:
        os.chdir(home)

        os.chdir(dr_B_NA)
        file_low, file_high = get_low_high()

        Plot.Spherical_Plot(file_low, tmax=1e8)
        Plot.Spherical_Plot(file_high, tmax=1e8)

    if "C_NA" in sys.argv:
        os.chdir(home)
        os.chdir(dr_C_NA)

        file_low, file_high = get_low_high()

        Plot.Spherical_Plot(file_low)
        Plot.Spherical_Plot(file_high)
        Plot.Alpha_Plot(file_low)
        Plot.Alpha_Plot(file_high)

    if "A" in sys.argv:
        os.chdir(home)
        os.chdir(dr_A)
        file_low, file_high = get_low_high()
        Plot.Spherical_Plot_Transform(file_low, {'save_fig': True})
        Plot.Spherical_Plot_Transform(file_high, {'save_fig': True})
        Plot.Angle_Space_Plot(file_high, {"EBF": True, "3D": True,
                                          'split': '0/1000/4000/6000/9000/10000',
                                          'save_fig': True})
        Plot.ThreeD_Plot_Cartesian(file_high, {'power': 0.5, 'start': 4000,
                                               'stop': 6000, 'EBF': True,
                                               'save_fig': True})

    if "B" in sys.argv:
        os.chdir(home)
        os.chdir(dr_B)
        file_low, file_high = get_low_high()

        Plot.Spherical_Plot_Transform(file_low, {'nmax': 1000,
                                                 'save_fig': True})
        Plot.Spherical_Plot_Transform(file_high, {'nmax': 1000,
                                                  'save_fig': True})
        Plot.Spherical_Plot(file_low, save_fig=True)
        Plot.Spherical_Plot(file_high, save_fig=True)

    if "C" in sys.argv:
        os.chdir(home)
        os.chdir(dr_C)
        file_low, file_high = get_low_high()
        Plot.Spherical_Plot_Transform(file_low, {'save_fig': True})
        Plot.Spherical_Plot_Transform(file_high, {'save_fig': True})

if "copy" in sys.argv:
    for dr in [dr_A, dr_B, dr_C, dr_A_NA, dr_B_NA, dr_C_NA]:
        os.system("cp %s/*.png Images/" % dr)
