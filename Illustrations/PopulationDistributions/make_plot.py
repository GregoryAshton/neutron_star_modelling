import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc_file
import pandas as pd
import sys
from scipy.stats import norm
import BDATools as BDA


def AddNormal(ax, data):
    mu = np.mean(data)
    std = np.std(data)
    quants = np.linspace(data.min(), data.max(), 100)
    ax.plot(quants, norm.pdf(quants, loc=mu, scale=std), lw=4, color="r")
    s = "Properties\n$\mu$={}\n$\sigma$={}".format(
        BDA.Texify_Float(mu), BDA.Texify_Float(std))
    ax.annotate(s, xy=(0.7, 0.8), xycoords="figure fraction", size=22,
                #bbox=dict(boxstyle="square, pad=0.3", fc="w", ec="k", lw=1)
                )
    return ax

rc_file("/home/greg/Neutron_star_modelling/matplotlibrc")

DATA_FILE = "ATNF_data_file.txt"

data = np.genfromtxt(DATA_FILE, skip_header=4, skip_footer=1, dtype=None)
name = data[:, 0]
F0 = np.genfromtxt(data[:, 1])
F1 = np.genfromtxt(data[:, 2])
F2 = np.genfromtxt(data[:, 3])
Binary = data[:, 4] # If not "*" then it is a binary
Type = data[:, 5]
W10 = data[:, 6]

df = pd.DataFrame({'name': name,
                   'F0': F0,
                   'F1': F1,
                   'F2': F2,
                   'Binary': Binary,
                   'Type': Type,
                   'W10': W10
                   })

# Clean data
df = df[df.Binary == "*"]
df = df[df.Type == "*"]
df = df[df.F0 < 5e1 ]

if "F0" in sys.argv or len(sys.argv) == 1:
    fig = plt.figure(figsize=(6, 8))
    ax = fig.add_subplot(111)
    logF0 = np.log10(df.F0[np.isfinite(df.F0.values)])
    ax.hist(logF0, bins=100, normed=True)
    ax = AddNormal(ax, logF0)
    ax.set_xlabel("$\log_{10}(f)$")
    ax.set_ylabel("Normalised count")
    fig.tight_layout()
    plt.savefig("F0_distribution.pdf")
    plt.show()


if "F1" in sys.argv or len(sys.argv) == 1:
    fig = plt.figure(figsize=(6, 8))
    ax = fig.add_subplot(111)
    F1 = np.abs(df.F1[np.isfinite(df.F1.values)]).values
    logF1 = np.log10(F1[F1 > 0])
    ax.hist(logF1, bins=100, normed=True)
    ax = AddNormal(ax, logF1)
    ax.set_xlabel("$\log_{10}(\dot{f})$")
    ax.set_ylabel("Normalised count")
    fig.tight_layout()
    plt.savefig("F1_distribution.pdf")
    plt.show()

if "F2" in sys.argv or len(sys.argv) == 1:
    fig = plt.figure(figsize=(6, 8))
    ax = fig.add_subplot(111)
    logF2 = np.log10(np.abs(df.F2[np.isfinite(df.F2.values)]).values)
    ax.hist(logF2, bins=100, normed=True)
    ax = AddNormal(ax, logF2)
    ax.set_xlabel("$\log_{10}(\ddot{f})$")
    ax.set_ylabel("Normalised count")
    fig.tight_layout()
    plt.savefig("F2_distribution.pdf")
    plt.show()

if "W10" in sys.argv or len(sys.argv) == 1:
    fig = plt.figure(figsize=(6, 8))
    ax = fig.add_subplot(111)
    log10W10 = np.log10(df.W10[df.W10.values != "*"].values.astype(np.float))
    ax.hist(log10W10, bins=100, normed=True)
    ax = AddNormal(ax, log10W10)
    ax.set_xlabel("$\log_{10}(W_{10})$")
    ax.set_ylabel("Normalised count")
    fig.tight_layout()
    plt.savefig("W10_distribution.pdf")
    plt.show()
