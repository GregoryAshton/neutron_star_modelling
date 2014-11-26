import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import cumtrapz
from detect_peaks import detect_peaks
import TNtools as TN

TN.PlotDefaults()

def SimulateBothMags(F1A=-7.4e-14, F1B=-5.0e-14, D=1e-10,
                     T=100*86400, R=0.6, F0_init=2.0, N=1000, plot=True):

    # Derived quantities
    tA = R * T
    tB = T - tA
    Tobs = 7 * T

    tA_list = np.random.normal(tA, D*tA, size=N)
    tB_list = np.random.normal(tB, D*tB, size=N)

    flip_markers = [sum(np.dstack((tA_list,tB_list)).flat[0:i]) for i in range(2*N)]

    time = np.linspace(0, Tobs, 10000)

    F1_list=[]
    current_F1, other_F1 = F1A, F1B
    j=0
    for t in time:
        if t < flip_markers[j]:
            F1_list.append(current_F1)
        if t >= flip_markers[j]:
            current_F1, other_F1 = other_F1, current_F1
            j+=1
            F1_list.append(current_F1)
   
    F0_list = 2 + cumtrapz(y=F1_list, x=time, initial=0)
    phase = 2 * np.pi * cumtrapz(y=F0_list, x=time, initial=0)

             
    # Fit polynomial to phase of order order
    coefs = np.polyfit(time, phase, 2)
    # poly1d returns the polynomial we then evaluate this at time giving the fitted phi
    phase_fit = np.poly1d(coefs)(time)

    # Subtract the two to get a residual
    phase_res_num = phase - phase_fit

    # Analytic calculation
    def DeltaPhi(t, dFdot_total, T, R):
        t = np.mod(t, T)
        if t < R * T:
            x1 = (1-R) * t * (t - R*T)
        if t >= R*T:
            x1 = R * (T - t) * (t - R*T)
            
        return np.pi * dFdot_total * x1 

    dFdot_total = abs(F1A - F1B)
    phase_res_anal = [DeltaPhi(t, dFdot_total, T, R) for t in time]

    # Get variability
    num_range = GetRange(phase_res_num, plot)
    anal_range = GetRange(phase_res_anal, plot)

    return num_range, anal_range


def GetRange(phase_res, show=False):
    phase_res = np.array(phase_res)
    maximums = np.array(detect_peaks(phase_res, mph=0, mpd=20, show=show))
    minimums = np.array(detect_peaks(phase_res, mph=0, mpd=20, valley=True, show=show))

    max_residuals = phase_res[maximums]
    min_residuals = phase_res[minimums]

    return np.average(max_residuals) - np.average(min_residuals)

if __name__ == "__main__":
    num_ranges = []
    anal_ranges = []
    Tvals = np.linspace(1*86400, 1000*86400, 100)
    for T in Tvals:
        (num_range, anal_range) =  SimulateBothMags(F1A=-7.4e-14, F1B=-5.0e-14, D=1e-10,
                    T=T, R=0.6, F0_init=2.0, N=1000, plot=False)
        num_ranges.append(num_range)
        anal_ranges.append(anal_range)

    ax = plt.subplot(111)
    ax.plot(Tvals, num_ranges, lw=2, label="Numeric", color="r")
    ax.plot(Tvals, anal_ranges, label="Analytic")
    ax.legend(loc=2)
    ax.set_xlabel("T")
    ax.set_ylabel("$\Delta\ddot{\Phi}$", rotation="horizontal")
    plt.savefig("img/TestingMagnitude.pdf")
    plt.show()


