import matplotlib.pyplot as plt 
import numpy as np
import emcee
import scipy.optimize as op
import triangle
import GetLyneData
import TNtools as TN
import BDATools as BDA

TN.PlotDefaults()

def PlotWithData(x, y, sampler, SignalModel, nsamples=10, c="#4682b4", 
                 x_sim=None, ax=None):
    """ Plot the data along with the Signal Model """

    ax.plot(x, y, "o", markersize=3, label="data")

    if x_sim == None:
        x_sim = x
    samples = sampler.flatchain
    for s in samples[np.random.randint(len(samples), size=nsamples)]:
        ysignal = SignalModel(s[:-1], x_sim)
        ax.plot(x_sim, ysignal, color=c, alpha=0.7)
        noise_low = ysignal - s[-1]
        noise_high = ysignal + s[-1]
        ax.fill_between(x_sim, noise_low, noise_high, color=c, alpha=0.05) 

    return ax

fig, (ax2, ax1) = plt.subplots(nrows=2, sharex=True, figsize=(6, 6))

x0, y = GetLyneData.GetData(D=None,
file_name="/home/greg/timing-noise/ExtractDataFromLyne/B1828_W10_01.txt")
x = x0 - x0[0]

# Square wave
from scipy.signal import square

model_name = "Square"

y0_bg = 8.0
A_bg = 5.1
f0_bg = 0.73 / 365.25
offset_bg = 0.0
sigma_bg = 0.5
theta_bestguess = [y0_bg, A_bg, f0_bg, offset_bg, sigma_bg]
symbols = ["$y_{0}$", "$A$", "f", "$\phi_{0}$", "$\sigma$"]

def SignalModel(theta, x):
    y0, A, f0, offset = theta
    R = 0.5
    return y0 + A * square(2*np.pi*f0*x + offset, duty=R)


# Uniform prior
theta_lims = [[7.0, 9.0], [0, 10], [0.002, 0.0022], [1.0, 2*np.pi], [0.0, 10.0]]

sampler = BDA.RunEmCee(x, y, SignalModel, theta_lims=theta_lims, nburn0=250)
lnprobs = sampler.lnprobability[:, :].reshape((-1))
value_square = np.log10(BDA.HarmonicMeanApprox(lnprobs))

x_sim = np.linspace(x[0], x[-1], 1000)
ax1 = PlotWithData(x, y, sampler, SignalModel, 10, x_sim=x_sim, ax=ax1)

# Sine model
model_name = "Sine"
y0_bg = 8.0
A_bg = 6.0
offset_bg = np.pi
f0_bg = 0.73 / 365.25
sigma_bg = 0.1
theta_bestguess = [y0_bg, A_bg, offset_bg, f0_bg, sigma_bg]
symbols = ["$y_{0}$", "$A$", "$\phi_{0}$", "f", "$\sigma$"]

def SignalModel(theta, x):
    y0, A, offset, f0 = theta
    return y0 + A*np.sin(2*np.pi*f0*x + offset) 

# Uniform prior
theta_lims = [[7.0, 9.0], [0, 10], [0, 2*np.pi], [0.0018, 0.0025] , [0.0, 10.0]]

sampler = BDA.RunEmCee(x, y, SignalModel, theta_lims=theta_lims, nburn0=250)
lnprobs = sampler.lnprobability[:, :].reshape((-1))
value_sine = np.log10(BDA.HarmonicMeanApprox(lnprobs))

ax2 = PlotWithData(x, y, sampler, SignalModel, 10, x_sim=x_sim, ax=ax2)


ax1.set_yticks(ax1.get_yticks()[:-1])
ax2.set_yticks(ax2.get_yticks()[1:-1])

ax2.set_xticklabels([int(i) for i in (53000 + np.array(ax2.get_xticks()))])
ax2.set_xlabel("Modified Julian Date [days]", size=14)
ax2.set_ylabel("$W_{10}$ [ms]", size=14)
ax1.set_ylabel("$W_{10}$ [ms]", size=14)

plt.tight_layout()
plt.subplots_adjust(hspace=0)

ax1.text(.5,.9,'Model B ',
        horizontalalignment='center',
        transform=ax1.transAxes)

ax2.text(.5,.9,'Model A ',
        horizontalalignment='center',
        transform=ax2.transAxes)
plt.savefig("img/ResultsCombined.pdf")
plt.show()

print value_sine - value_square


