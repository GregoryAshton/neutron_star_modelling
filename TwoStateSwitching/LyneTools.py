import matplotlib.pylab as plt
import numpy as np
from scipy.integrate import cumtrapz

def GenerateModel(R=3.0):
    # Model parameters
    t1 = 100 * 24 * 3600
    D = 1e-15
    Frac = 0.71
    nu_dot_1 = 1e-10

    T = 50 * t1
    N = 1000 # Upper limit of flips 

    nu_dot_2 = Frac * nu_dot_1
    t2 = R * t1

    t1_list = np.random.normal(t1, D*t1, size=N)
    t2_list = np.random.normal(t2, D*t2, size=N)

    flip_markers = [sum(np.dstack((t1_list,t2_list)).flat[0:i]) for i in range(2*N)]

    t_list = np.linspace(0, T, 1000)

    nu_dot_list=[]
    current_nu_dot, other_nu_dot = nu_dot_1, nu_dot_2
    j=0
    for t in t_list:
        if t < flip_markers[j]:
            nu_dot_list.append(current_nu_dot)
        if t >= flip_markers[j]:
            current_nu_dot, other_nu_dot = other_nu_dot, current_nu_dot
            j+=1
            nu_dot_list.append(current_nu_dot)

    # Plotting 
    #fig, (ax1, ax2) = plt.subplots(nrows=2, sharex=True)
    #ax1.plot(t_list, nu_dot_list, "-")
    #ax1.set_ylabel(r"$\dot{\nu}$", rotation="horizontal")
    #ax1.set_xticklabels([])
    #ax1.set_yticks([nu_dot_1, nu_dot_2])
    #ax1.set_yticklabels([r"$\dot{\nu}_{1}$", r"$\dot{\nu}_{2}$"])
      
        
    nu = 2 + cumtrapz(y=nu_dot_list, x=t_list, initial=0)
    phase = cumtrapz(y=nu, x=t_list, initial=0)

             
    # Fit polynomial to phase of order order
    coefs = np.polyfit(t_list, phase, 2)
    # poly1d returns the polynomial we then evaluate this at time giving the fitted phi
    phase_fit = np.poly1d(coefs)(t_list)

    # Subtract the two to get a residual
    T_res = phase - phase_fit

    #ax2.plot(t_list, T_res)
    #ax2.set_ylabel("Phase residual \n [cycles]")
    #ax2.set_xlabel(r"time")
    #ax2.set_yticks(ax2.get_yticks()[0:-1])

    #ax1.set_xlim(0, t_list[-100])
    #ax2.set_xlim(0, t_list[-100])

    #if D < 1e-9:
    #    D = 0
    #ax1.text(0.08*t_list[-1], 1.3 * ax1.get_ylim()[-1], 
    #         "R={} \n D={}".format(R, D), 
    #         size=14, 
    #         bbox=dict(facecolor='white', edgecolor="white", alpha=0.9))


    #plt.subplots_adjust(hspace=0.0, wspace=0.5)
    #plt.show()

    return T_res

def MagnitudeOfPeaks(Tres):
    """ Calculate the magnitude of the peaks in the residual """
    Pmax = np.max(Tres)
    Pmin = np.min(Tres)
    return abs(Pmax/Pmin)

R = 3.0

x = []
y = []
for R in np.linspace(1.1, 10, 100):
   x.append(R/(R-1))
   y.append(MagnitudeOfPeaks(GenerateModel(R)))

plt.plot(x, y)
plt.show()

