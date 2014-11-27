""" Tools to play around with the extracted data from Lyne 2010 

Data extracted using the Grabit Matlab tools by Jiro Doke 

http://www.mathworks.com/matlabcentral/fileexchange/7173-grabit

"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

file_name = "Fig5_EF_01.txt"

def GetData(file_name):
    df = pd.read_csv(file_name, sep=" ", header=None, names=['MJD', 'W10_ms'],
                     dtype=None, skipinitialspace=True, index_col=False)

    # The extraction process is not always in the increasing MJD order
    df = df.sort('MJD')

    return df

def Plot(file_name):
    df = GetData(file_name)

    plt.plot(df.MJD.values, df.W10_ms.values, "-o")
    plt.xlabel("MJD")
    plt.ylabel("$W_{10}$ (ms)")
    plt.show()

if __name__ == "__main__":
    Plot(file_name) 
