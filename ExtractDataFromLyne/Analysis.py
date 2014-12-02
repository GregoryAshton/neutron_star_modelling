""" Tools to play around with the extracted data from Lyne 2010 

Data extracted using the Grabit Matlab tools by Jiro Doke 

http://www.mathworks.com/matlabcentral/fileexchange/7173-grabit

"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
import TNtools as TN
import sys

TN.PlotDefaults()

file_name = "B1828_W10_01.txt"

def GetData(file_name):
    df = pd.read_csv(file_name, sep=" ", header=None, names=['MJD', 'W10_ms'],
                     dtype=None, skipinitialspace=True, index_col=False)

    # The extraction process is not always in the increasing MJD order
    df = df.sort('MJD')

    return df

def Plot(data, ax=None, **kwargs):
    """ Plot the data in file_name on this axis """
   
    if not ax:
        ax = plt.subplot(111)

    if type(data) == str:
        df = GetData(data)
    elif isinstance(data, pd.DataFrame):
        df = data

    ax.plot(df.MJD.values, df.W10_ms.values, "-o", **kwargs)
    ax.set_xlabel("MJD")
    ax.set_ylabel("$W_{10}$ (ms)")

    ax.set_xlim(df.MJD.values[0], df.MJD.values[-1])

    return ax

def LeastSquaresCurveFit(f, xdata, ydata, p0=None, maxfev=10000):
    """ Use Non-linear least sqaures curve fit the data using function f """

    popt, pcov = curve_fit(f, xdata, ydata, p0=p0, maxfev=maxfev)

    perr = np.sqrt(np.diag(pcov))
   
    return popt, perr

def CompareSinSquare(file_name):
    """ Compare fitting a square wavw with a sinusoid for data in file_name """

    df = GetData(file_name)

    ax = plt.subplot(111)
    
    # Subtract off the mean to help fitting
    xdata = df.MJD.values
    ydata = df.W10_ms.values
    ydata = ydata - np.mean(ydata)
    ax.plot(xdata, ydata, "-.", label="Data")

    def Sin(x, fA, A, fB, B):
        return A * np.sin(fA * x) + B * np.cos(fB *x)

    p0 = [0.1, 3, 0.1, 3]
    popt_sin, perr_sin = LeastSquaresCurveFit(Sin, xdata, ydata, p0=p0)
    print popt_sin, perr_sin

    MJD_fine = np.linspace(df.MJD.values[0], df.MJD.values[-1], 1000)
    ax.plot(MJD_fine, Sin(MJD_fine, *popt_sin))

    plt.show()

def LyneFit():
    """ Plot the B1828-11 fit from Lyne """

    file_name = "B1828_W10_01.txt"
    df = GetData(file_name)
    
    time = df.MJD.values 
    #F = 0.73
    #T = 1./F * 86400 * 365.25
    T = 450
    print T

    def W10LyneModel(time, T, R):

        NumberFlip = 100 * int((time[-1] - time[0]) / float(T))
        N = len(time)

        wA = 5.9
        wB = 10
        tA = R * T
        tB = (1-R) * T

        tA_list = np.zeros(NumberFlip) + tA
        tB_list = np.zeros(NumberFlip) + tB

        flip_markers = np.array([
                          int(sum(np.dstack((tA_list,tB_list)).flat[:i]))
                                                    for i in range(2*N)])
        flip_markers += time[0]

        w_list = []
        current_w, other_w = wA, wB
        j=0
        for t in time:
            if t < flip_markers[j]:
                w_list.append(current_w)
            if t >= flip_markers[j]:
                current_w, other_w = other_w, current_w
                j+=1
                w_list.append(current_w)

        return time, w_list

    ax = plt.subplot(111)

    ax.plot(time, df.W10_ms.values, "-o")
    
    time_fine = np.linspace(time[0], time[-1], 100)
    R = 0.5
    time, w_list = W10LyneModel(time_fine, T, R) 
    ax.plot(time_fine, w_list)

    plt.show()
   
def BiningStates(file_name, bins=50):
    """ Look at the population of states """

    df = GetData(file_name)
    W10 = df.W10_ms.values

    hist, bin_edges = np.histogram(W10, bins)

    ax = plt.subplot(111)
    width = bin_edges[1] - bin_edges[0]
    ax.bar(bin_edges[:-1], hist, width=width)


    ax.set_xlabel("Binned $W_{10}$[ms] values")
    ax.set_ylabel("Count")
    plt.savefig("img/Histogram"+file_name.replace("txt", "pdf"))

    plt.show()


if __name__ == "__main__":
    if "plt" in sys.argv:
        fig, ax = plt.subplots(figsize=(15, 6))
        ax = Plot(file_name, markersize=3.0, ax=ax)
        plt.savefig("img/ExtractedData_"+file_name.replace("txt", "pdf"))
        plt.show()
    if "fit" in sys.argv:
        CompareSinSquare(file_name) 
    if "bin" in sys.argv:
        BiningStates(file_name)
    if "lyne" in sys.argv:
        LyneFit()
