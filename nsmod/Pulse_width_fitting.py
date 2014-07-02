import numpy as np
from scipy.integrate import cumtrapz, simps
from scipy.signal import find_peaks_cwt
import matplotlib.pyplot as plt

def CleanAmp(time, amp):
    """ Given random cuttings of amplitude return a cleaned version which 
    does not include half peaks
    """

    threshold = 1e-6 * np.max(amp)
    idx_first = np.where(amp < threshold )[0][0]
    idx_last = len(time) - np.where(amp[::-1] < threshold)[0][0]
    return time[idx_first:idx_last], amp[idx_first:idx_last]


def CalculateC(time, amp, plot=False):
    """ Measure C from a number of gaussian peaks """

    time, amp = CleanAmp(time, amp)

    # Measure the number of peaks
    peakind = find_peaks_cwt(amp, np.array([90]))
  

    N = 0
    max_amp = np.max(amp)
    for p in peakind:
        if amp[p] > 0.1 * max_amp:
            N +=1
    #print N
    #print len(peakind)
     
    #plt.plot(time, amp)
    #for p in peakind:
    #    plt.axvline(time[p], color="r")
    #plt.show()
    
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

import numpy as np
from itertools import chain, izip


def find(a, predicate, chunk_size=1024):
    """

    https://github.com/numpy/numpy/issues/2269

    Find the indices of array elements that match the predicate.

    Parameters
    ----------
    a : array_like
        Input data, must be 1D.

    predicate : function
        A function which operates on sections of the given array, returning
        element-wise True or False for each data value.

    chunk_size : integer
        The length of the chunks to use when searching for matching indices.
        For high probability predicates, a smaller number will make this
        function quicker, similarly choose a larger number for low
        probabilities.

    Returns
    -------
    index_generator : generator
        A generator of (indices, data value) tuples which make the predicate
        True.

    See Also
    --------
    where, nonzero

    Notes
    -----
    This function is best used for finding the first, or first few, data values
    which match the predicate.

    Examples
    --------
    >>> a = np.sin(np.linspace(0, np.pi, 200))
    >>> result = find(a, lambda arr: arr > 0.9)
    >>> next(result)
    ((71, ), 0.900479032457)
    >>> np.where(a > 0.9)[0][0]
    71


    """
    if a.ndim != 1:
        raise ValueError('The array must be 1D, not {}.'.format(a.ndim))

    i0 = 0
    chunk_inds = chain(xrange(chunk_size, a.size, chunk_size), 
                 [None])

    for i1 in chunk_inds:
        chunk = a[i0:i1]
        for inds in izip(*predicate(chunk).nonzero()):
            yield (inds[0] + i0, ), chunk[inds]
        i0 = i1
