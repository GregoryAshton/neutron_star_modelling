import TNtools as TN
TN.PlotDefaults()

import matplotlib.pyplot as plt
import numpy as np
import sys
from scipy.integrate import cumtrapz

def get_F1(tA, T, N, F1A, F1B):
    tA_idx = int(tA * N / float(T))
    time = np.linspace(0, T, N)
    F1 = np.zeros(len(time))
    F1[:tA_idx] = F1A
    F1[tA_idx:] = F1B
    return time, F1

def shade_regions(ax, tA):
    ylim = ax.get_ylim()

    ax.fill_between([0, tA], ylim[0], ylim[1],
                    color="k", alpha=0.2, zorder=-3)
    return ax

def F1Plot(R=0.4, T=10, N=1000, F1A=-1, F1B=-2):

    tA = R * T
    time, F1 = get_F1(tA, T, N, F1A, F1B)
    ax = plt.subplot(111)
    ax.plot(time, F1, ls="-", color="k", lw=2)
    ax.set_ylim(-3, 0)
    ax.axhline(-1.6, color="k", ls="--")

    ax.set_xlabel("time")
    trefs = [tA/2.0, T*(R+1)/2.0]
    for t in trefs:
        ax.axvline(t, lw=0.5, color="k")
    ax.set_xticks([0, trefs[0], trefs[1], T])
    ax.set_xticklabels(['0', 'RT/2', 'T(R+1)/2', 'T'])
    ax.set_yticks([0, -1.0, -1.6, -2.0])
    ax.set_yticklabels(['0', '$\dot{f}_{G} + \Delta \dot{f}_{A}$', 
                        '$\dot{f}_{G}$', '$\dot{f}_{G} + \Delta \dot{f}_{B}$'])

    offset = 0.1
    text_offset = 0.15
    x_offset = 0.1
    ax.annotate("", xy=(0, F1A+offset), xytext=(tA, F1A+offset),
            arrowprops=dict(arrowstyle="<|-|>", fc='k' ))
    ax.text(x_offset + tA/2., F1A+text_offset, r'$RT$')
    
    offset = -0.1
    text_offset = -0.25
    ax.annotate("", xy=(tA, F1B+offset), xytext=(T, F1B+offset),
            arrowprops=dict(arrowstyle="<|-|>", fc='k' ))
    ax.text(x_offset + (T+tA)/2., F1B+text_offset, r'$(1-R)T$')

    ax = shade_regions(ax, tA) 
    plt.savefig("img/F1.pdf")

def F0Plot(tA=4, T=10, N=1000, F1A=-1, F1B=-2, F0_init=30):

    time, F1 = get_F1(tA, T, N, F1A, F1B)
    F0 = cumtrapz(F1, time, initial=0)
    F0 += F0_init

    ax = plt.subplot(111)
    
    ax.plot(time, F0)

    ax = shade_regions(ax, tA)
    
def P0Plot(tA=4, T=10, N=1000, F1A=-1, F1B=-2, F0_init=30):

    time, F1 = get_F1(tA, T, N, F1A, F1B)
    F0 = cumtrapz(F1, time, initial=0)
    F0 += F0_init
    
    P0 = cumtrapz(F0, time, initial=0)
    ax = plt.subplot(111)

    ax.plot(time, P0)

def DeltaPhi(t, dST, T, R):
    t_tilde = np.mod(t, T) 
    R = float(R)
    if t_tilde <= T*R:
        post = (1-R) * t_tilde * (t_tilde-R*T) 
    elif t_tilde > T*R:
        post =  -R * (t_tilde-T) * (t_tilde-R*T)

    return np.pi * dST * post

def DPPlot(dST=10, T=20, R=0.8, N=1000):
    " Plot the phase residual "

    time = np.linspace(0, T, N)
    DP = [DeltaPhi(t, dST, T, R) for t in time]

    ax = plt.subplot(111)
    ax.plot(time, DP, lw=2)

    #\ax.axhline(0, ls="-", color="k")
    #yticks = [DeltaPhi(t, dST, T, R) for t in [.5*T*R, T, T*(R+1)/2.]]
    yticks = [-np.pi * (R*T/2)**2 *(1-R) * dST, 0, np.pi * R*T**2/4.0 * (R-1)**2 * dST]
    ax.set_yticks(yticks)
    ax.set_yticklabels([r"$-\pi\left(\frac{RT}{2}\right)^{2}(1-R)\Delta\dot{f}_{T}$",
                        r"0",
                        r"$\frac{\pi}{4} RT^{2}(R-1)^{2}\Delta\dot{f}_{T}$"],
                        fontsize=16)

    xticks = [0, R*T/2, R*T, T*(R+1)/2., T]
    ax.set_xticks(xticks)
    ax.set_xticklabels(["0", r"$RT/2$", "$RT$", 
                        r"$T(R+1)/2$", "$T$"],
                        fontsize=18)

    ax.grid()
    ax.set_ylabel("$\Delta\Phi(t)$", labelpad=-100, rotation='horizontal') 
    ax = shade_regions(ax, R*T)
    ax.text(0.2*T, .8*yticks[2], "R={}".format(R), size=22)
    plt.subplots_adjust(left=0.25)
    plt.savefig("img/DP.pdf")
   
def Mismatch():
    dST = 1
    T = 10
    R = np.linspace(0, 1, 100)
    m = (np.pi**2 * dST * T**4  * (1-R)**2 * R**2 * 
                   (-2 * R**2 + 2 * R + 1)/180.)

    ax = plt.subplot(111)
    ax.plot(R, m)
    
    xticks = [0, .5, 1]
    ax.set_xticks(xticks)
    ax.set_xticklabels(["0", r"$1/2$", "$1$"], fontsize=18)
    ax.set_xlabel("R")

    #ax.set_yticks([0, np.pi**2 * dST * T**4  * (3/32.) / 180.])
    #ax.set_yticklabels(['0', r'$\frac{\pi^{2}}{180} \Delta \dot{f}_{T} T^{4} \frac{3}{32}$'])
    ax.set_yticks([])
    ax.set_ylabel("Mismatch")
    plt.savefig("img/mismatch.pdf")
    plt.show()


if __name__ == "__main__":
    if "F1" in sys.argv:
        F1Plot()
    if "F0" in sys.argv:
        F0Plot()
    if "P0" in sys.argv:
        P0Plot(tA=4, T=10, N=1000, F1A=0, F1B=-5, F0_init=30)
    if "DP" in sys.argv:
        DPPlot(dST=1e-3, T=10, R=0.4, N=100)
    if "m" in sys.argv:
        Mismatch()


    plt.show()
