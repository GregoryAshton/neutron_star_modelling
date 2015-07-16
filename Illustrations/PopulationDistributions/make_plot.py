import matplotlib.pyplot as plt
import numpy as np
import sys
from scipy.stats import norm
import BDATools as BDA
from get_data import df


def AddNormal(ax, data):
    mu = np.mean(data)
    std = np.std(data)
    quants = np.linspace(data.min(), data.max(), 50)
    ax.plot(quants, norm.pdf(quants, loc=mu, scale=std), lw=4, color="r")
    s = "Properties\n$\mu$={}\n$\sigma$={}".format(
        BDA.Texify_Float(mu, 4), BDA.Texify_Float(std, 4))
    ax.annotate(s, xy=(0.7, 0.8), xycoords="figure fraction", size=22,
                # bbox=dict(boxstyle="square, pad=0.3", fc="w", ec="k", lw=1)
                )
    return ax


def setup_hist(vals, ax, nbins=30):
    ax.hist(vals, nbins, normed=True, histtype="step")
    ax.hist(vals, nbins, normed=True, color="b", histtype="stepfilled",
            alpha=0.2)
    return ax

if "F0" in sys.argv or len(sys.argv) == 1:
    fig = plt.figure(figsize=(6, 8))
    ax = fig.add_subplot(111)
    logF0 = np.log10(df.F0[np.isfinite(df.F0.values)])
    ax = setup_hist(logF0, ax, 100)
    ax = AddNormal(ax, logF0)
    ax.set_xlabel("$\log_{10}(f)$")
    ax.set_ylabel("Normalised count")
    fig.tight_layout()
    plt.savefig("F0_distribution.pdf")


if "F1" in sys.argv or len(sys.argv) == 1:
    fig = plt.figure(figsize=(6, 8))
    ax = fig.add_subplot(111)
    F1 = np.abs(df.F1[np.isfinite(df.F1.values)]).values
    logF1 = np.log10(F1[F1 > 0])
    ax = setup_hist(logF1, ax, 100)
    ax = AddNormal(ax, logF1)
    ax.set_xlabel("$\log_{10}(\dot{f})$")
    ax.set_ylabel("Normalised count")
    fig.tight_layout()
    plt.savefig("F1_distribution.pdf")

if "F2" in sys.argv or len(sys.argv) == 1:
    fig = plt.figure(figsize=(6, 8))
    ax = fig.add_subplot(111)
    logF2 = np.log10(np.abs(df.F2[np.isfinite(df.F2.values)]).values)
    ax = setup_hist(logF2, ax, 100)
    ax = AddNormal(ax, logF2)
    ax.set_xlabel("$\log_{10}(\ddot{f})$")
    ax.set_ylabel("Normalised count")
    fig.tight_layout()
    plt.savefig("F2_distribution.pdf")

if "log10W10" in sys.argv or len(sys.argv) == 1:
    fig = plt.figure(figsize=(6, 8))
    ax = fig.add_subplot(111)
    log10W10 = np.log10(df.W10[np.isfinite(df.W10)])
    ax = setup_hist(log10W10, ax, 50)
    ax = AddNormal(ax, log10W10)
    ax.set_xlabel("$\log_{10}(W_{10})$")
    ax.set_ylabel("Normalised count")
    fig.tight_layout()
    plt.savefig("log10W10_distribution.pdf")

if "lnW10" in sys.argv or len(sys.argv) == 1:
    fig = plt.figure(figsize=(6, 8))
    ax = fig.add_subplot(111)
    lnW10 = np.log(df.W10[np.isfinite(df.W10)])
    ax = setup_hist(lnW10, ax, 30)
    ax = AddNormal(ax, lnW10)
    ax.set_xlabel("$\ln(W_{10})$")
    ax.set_ylabel("Normalised count")
    fig.tight_layout()
    plt.savefig("lnW10_distribution.pdf")
