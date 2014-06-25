from nsmod.switching_torque_with_Euler import main
from nsmod import Plot
import matplotlib.pyplot as plt

switching = main(epsI1=0.0, epsI3=3.0e-3, epsA=5.0e-4 , omega0=10,
                error=1e-12, T=1.0e3 , chi0 = 70.0, AnomTorque=True ,
                a0=20.0, upsilon=0.8, n=20000)

noswitch = main(epsI1=0.0, epsI3=3.0e-3, epsA=5.0e-4 , omega0=10,
                error=1e-12, T=1.0e3 , chi0 = 70.0, AnomTorque=True ,
                a0=20.0, upsilon=0.0, n=20000)

fig = Plot.Spherical_Plot(switching)
fig = Plot.Spherical_Plot(noswitch, fig=fig)
plt.show()


ax1 = plt.subplot(311)
ax2 = plt.subplot(312)
ax3 = plt.subplot(313) 

(ax1, ax2, ax3) = Plot.Euler_Angles(switching, ax_tup=(ax1, ax2, ax3), 
                                    label=r"$\upsilon=0.8$")
(ax1, ax2, ax3) = Plot.Euler_Angles(noswitch, ax_tup=(ax1, ax2, ax3), 
                                    label=r"$\upsilon=0$")

ax2.legend(loc=4)

plt.show()
