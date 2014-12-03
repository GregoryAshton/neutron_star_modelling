import matplotlib.pyplot as plt 
import numpy as np
import TNtools as TN
import sys
from scipy.integrate import simps

TN.PlotDefaults()

def P(y, y0, A):
    if y0 - A < y < y0 + A:
        return pow(1 - pow((y-y0)/A, 2), -0.5) / (2 * np.pi * A)
    else :
        return 0

if "sin" in sys.argv:
    " Plot just the sin distribution"
    y = np.linspace(0, 10, 600)
    y0 = 5
    A = 4

    Py = [P(yi, y0, A) for yi in y]

    ax = plt.subplot(111)
    ax.plot(y, Py)
    ax.fill_between(y, 0, Py, color="k", alpha=0.2)

    ax.set_xlabel("$y$")
    ax.set_xticks([y0 - A, y0, y0+A])
    ax.set_xticklabels(["$y_{0} - A$", "$y_{0}$", "$y_{0}+A$"], size=18)

    ax.set_ylabel("$f(y; y_{0}, A)$")
    ax.set_yticks([])
    ax.set_ylim(0, 0.3)


    plt.savefig("img/SinDistribution.pdf")

    plt.show()

if "noise" in sys.argv:
    " Plot the sin and noise "

    y0 = 5
    A = 4
    sigma = 0.3
    y = np.linspace(0, 10, 600)

    def N(eps, sigma):
        return np.exp(-eps**2/(2*sigma**2)) / (sigma * np.sqrt(2 * np.pi))
   
    Py_noise = []

    for yi in y:
        ytilde = y
        Integrand = [P(ytildei, y0, A) * N(yi - ytildei, sigma) 
                                                  for ytildei in ytilde]
        I = simps(Integrand, x=y)

        Py_noise.append(I)

    ax = plt.subplot(111)

    ax.plot(y, Py_noise)
    ax.fill_between(y, 0, Py_noise, color="k", alpha=0.2)

    ax.set_xlabel("$y$")
    ax.set_xticks([y0 - A, y0, y0+A])
    ax.set_xticklabels(["$y_{0} - A$", "$y_{0}$", "$y_{0}+A$"], size=18)

    ax.set_ylabel("$P(y^{M}; y_{0}, A, \sigma)$")
    ax.set_yticks([])

    plt.savefig("img/SinDistribution_withnoise.pdf")
    plt.show()
