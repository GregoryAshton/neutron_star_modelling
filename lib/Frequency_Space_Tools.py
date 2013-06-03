#!/usr/bin/python

import pylab as py
import numpy as np


def FFT(time, data):
    """ Returns a list of amplitudes and phase at the given frequencies """

    n = len(time)
    frequency = np.fft.fftfreq(n, d=time[-1] / n)

    complex_amplitude = np.fft.fft(data)
    magnitude = [py.norm(a) for a in complex_amplitude]
    phase = [py.arctan(a.imag / a.real) for a in complex_amplitude]

    return (frequency[0:n / 2], magnitude[0:n / 2], phase[0:n / 2])


def Peak_Find(data, box_size=10):
    """ Function to find true maximums in FFT amplitudes"""

    rough_list = []
    index_list = []
    n = len(data)

    for i in xrange(box_size, n - box_size, box_size):
        fprime_left = data[i] - data[i - box_size]
        fprime_right = data[i + box_size] - data[i]

        if fprime_left * fprime_right < 0.0 and fprime_left > 0.0:
            rough_list.append(i)

    for i in rough_list:
        max_value = max(data[i - box_size:i + box_size])
        max_index = data.index(max_value)
        index_list.append(max_index)

    return index_list