#!/usr/bin/python
"""

Generic frequency space tools. The Numpy standard library is used to
transform into the frequency domain.

"""
import pylab as py
import numpy as np


def FFT(time, data):
    """ Returns a list of magnitudes and phase at the given frequencies """

    n = len(time)
    frequency = np.fft.fftfreq(n, d=time[-1] / n)

    complex_amplitude = np.fft.fft(data)
    magnitude = [np.sqrt(a.real**2 + a.imag**2) for a in complex_amplitude]
    phase = [np.arctan(a.imag / a.real) for a in complex_amplitude]

    return (frequency[0:n / 2], magnitude[0:n / 2], phase[0:n / 2])


def Peak_Find(magnitude, box_size=10):
    """ Function to find true maximums in FFT magnitudes"""

    rough_list = []
    index_list = []
    n = len(magnitude)

    for i in xrange(box_size, n - box_size, box_size):
        fprime_left = magnitude[i] - magnitude[i - box_size]
        fprime_right = magnitude[i + box_size] - magnitude[i]

        if fprime_left * fprime_right < 0.0 and fprime_left > 0.0:
            rough_list.append(i)

    for i in rough_list:
        max_value = max(magnitude[i - box_size:i + box_size])
        max_index = magnitude.index(max_value)
        index_list.append(max_index)

    return index_list
