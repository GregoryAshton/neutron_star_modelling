import numpy as np
from scipy.integrate import cumtrapz, simps
from scipy.signal import find_peaks_cwt
import matplotlib.pyplot as plt

def CalculateC(time, amp, plot=False):
    """ Measure C from a number of gaussian peaks """

    # Measure the number of peaks
    peakind = find_peaks_cwt(amp, np.array([1]))
    N = len(peakind)

    # Checks
    if N == 0:
        print "WARNING: No peaks found"
        return 0
    elif len(time) < N * 10:
        print "WARNING: C calculation poorly resolved"

    # Measure A
    A = np.max(amp)

    # Integrate 
    #integral = simps(amp, time)
    integral = cumtrapz(amp, time, initial=0)[-1]

    # Calculate C
    C = integral / (N * A * np.sqrt(np.pi))

    return C

def W50(time, amp, tCONS):
    """ Calculate the FWH or W50 for time and amp

    Parameters
    ----------
    time, amp : array_like
        Arrays of the time and amplitude of pulses
    tCONS : float
        A float of a typical timescale over which the pulse properties should
        be constant
    """
    n = np.argmin(np.abs(time - tCONS)) - 10
    
    C_list = []
    #av_times = time[int(n/2.)::n]
    av_times = []
    for i in xrange(0, len(time)-n, n):
        C_list.append(CalculateC(time[i:i+n], amp[i:i+n]))
        av_times.append(time[i+int(n/2.)])


    return av_times, np.array(C_list)


