""" Tools to play around with the extracted data from Lyne 2010 

Data extracted using the Grabit Matlab tools by Jiro Doke 

http://www.mathworks.com/matlabcentral/fileexchange/7173-grabit

"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
import TNtools as TN

TN.PlotDefaults()

file_name = "B1828_W10_01.txt"

def GetData(file_name):
    df = pd.read_csv(file_name, sep=" ", header=None, names=['MJD', 'W10_ms'],
                     dtype=None, skipinitialspace=True, index_col=False)

    # The extraction process is not always in the increasing MJD order
    df = df.sort('MJD')

    return df

def Plot(df, ax=plt.subplot(111)):
    """ Plot the data in file_name on this axis """
    
    #if not ax:
    #    ax = plt.subplot(111)

    ax.plot(df.MJD.values, df.W10_ms.values, "-o")
    ax.set_xlabel("MJD")
    ax.set_ylabel("$W_{10}$ (ms)")

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

    MJD_fine = np.linspace(df.MJD.values[0], df.MJD.values[-1], 10000)
    ax.plot(MJD_fine, Sin(MJD_fine, *popt_sin))

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
    plt.savefig("img/"+file_name.replace("txt", "pdf"))

    plt.show()


if __name__ == "__main__":
    #CompareSinSquare(file_name) 
    BiningStates(file_name)
