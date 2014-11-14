#!/usr/bin/python

"""
Script to analyse the raw data extracted from the Hobbs 2010 paper

"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import TNtools as TN
import sys

TN.PlotDefaults()

df = pd.read_csv("HobbsLyneKramer2010Figure3Data.dat", sep=" ")
data_sets = ["001", "002", "003", "004", "005", "006"]
data_files = ["stripped-cleaned-images-{}.txt".format(s) for s in data_sets]
quad_files = data_files[:4]
cubic_files = data_files[4:]


if "check" in sys.argv:
    ax = plt.subplot(111)
    colors = [np.random.uniform(0, 1, size=3) for i in data_sets]
    for i in range(4):
        colors[i][0] = 0
    for i in range(4, 6):
        colors[i][1:3] = [0.2, 0.5]



    for data_set, color in zip(data_sets, colors):
                         
        file_name = "stripped-cleaned-images-{}.txt".format(data_set)
        
        Log_peaks_s = np.log10(df[df.source_file == file_name].peaks_s.values)
        hists,binEdges=np.histogram(Log_peaks_s,bins=20)
        bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
        ax.plot(bincenters, hists, label="Data from set {}".format(data_set),
                color=color, lw=3)
        ax.fill_between(bincenters, 0, hists, color=color, alpha=0.2)

    ax.legend(fontsize=12, frameon=False)
    ax.set_xlabel("$\log_{10}$ of peak to peak residual [s]")
    ax.set_ylabel("Count")
    #ax.set_xscale("log")
    plt.tight_layout()
    plt.savefig("img/check.pdf")
    plt.show()

if "quad" in sys.argv:
    df = df[df.source_file.isin(quad_files)]

    ax = plt.subplot(111)
    Log_peaks_s = np.log10(df.peaks_s.values)
    ax.hist(Log_peaks_s, bins=50)

    ax.legend(fontsize=8)
    ax.set_xlabel("$\log_{10}$ of peak to peak residual [s]")
    ax.set_ylabel("Count")
    #ax.set_xscale("log")

    plt.show()

if "quad-lyne" in sys.argv:
    df = df[df.source_file.isin(quad_files)]

    lyne_vals = ["1931+24", "2035+36", "1903+07", 
                 "1822-09", "1642-03", "1839+09", 
                 "1540-06", "2148+63", "1818-04",
                 "0950+08", "1714-34", "1907+00", 
                 "1828-11", "1826-17", "0919+06",
                 "0740-28", "1929+20"]

    #df_lyne = df[df.names.isin(lyne_vals)]

    ax = plt.subplot(111)
    Log_peaks_s = np.log10(df.peaks_s.values)
    ax.hist(Log_peaks_s, bins=50)

    ax.legend(fontsize=8)
    ax.set_xlabel("$\log_{10}$ of peak to peak residual [s]")
    ax.set_ylabel("Count")
    #ax.set_xscale("log")

    for name in lyne_vals:
        entry = df[df.names == name]
        if len(entry) == 0:
            entry = df[df.names == name.replace("+", "_").replace("-", "_")]

        if len(entry) ==0:
            print name, "not found"
        else:
            lyne_pulsars = ax.axvline(np.log10(entry.peaks_s), label="Lyne 2010 pulsars")
    plt.legend(handles=[lyne_pulsars]) 
    plt.title("Data taken from figure 3")
    plt.savefig("img/Figure3Histogram.pdf")
    plt.show()

    print "Mean of data set: {}".format(np.mean(Log_peaks_s))

if "cubic-lyne" in sys.argv:
    df = df[df.source_file.isin(cubic_files)]

    lyne_vals = ["1931+24", "2035+36", "1903+07", 
                 "1822-09", "1642-03", "1839+09", 
                 "1540-06", "2148+63", "1818-04",
                 "0950+08", "1714-34", "1907+00", 
                 "1828-11", "1826-17", "0919+06",
                 "0740-28", "1929+20"]

    #df_lyne = df[df.names.isin(lyne_vals)]

    ax = plt.subplot(111)
    Log_peaks_s = np.log10(df.peaks_s.values)
    ax.hist(Log_peaks_s, bins=50)

    ax.legend(fontsize=8)
    ax.set_xlabel("$\log_{10}$ of peak to peak residual [s]")
    ax.set_ylabel("Count")
    #ax.set_xscale("log")

    for name in lyne_vals:
        entry = df[df.names == name]
        if len(entry) == 0:
            entry = df[df.names == name.replace("+", "_").replace("-", "_")]

        if len(entry) ==0:
            print name, "not found"
        else:
            lyne_pulsars = ax.axvline(np.log10(entry.peaks_s), label="Lyne 2010 pulsars")
    plt.legend(handles=[lyne_pulsars]) 
    plt.title("Data taken from figure 13")

    plt.savefig("img/Figure13Histogram.pdf")
    plt.show()

    print "Mean of data set: {}".format(np.mean(Log_peaks_s))

