import matplotlib.pyplot as plt
import numpy as np
import sys
from scipy.special import ellipj, ellipk
import TNtools as TN

TN.PlotDefaults()

k2 = 0.5
K = ellipk(k2)
tau = np.linspace(0, 4*K, 100)

(sn, cn, dn, ph) = ellipj(tau, k2)
ax = plt.subplot(111)
ax.plot(tau, dn)

minimum = np.sqrt(1-k2)

ax.set_xticks([0, K, 2*K, 3*K, 4*K])
ax.set_xticklabels(["$0$", "$K$", "$2K$", "$3K$", "$4K$"])
ax.set_xlim(0, 4*K)
ax.set_xlabel(r"$\tau$")

ax.set_yticks([0, minimum, 1])
ax.set_yticklabels(['$0$', "$\sqrt{1-k^{2}}$", "1"])
ax.set_ylim(0, 1)
ax.set_ylabel(r"$\mathrm{dn}(\tau, k)$", labelpad=-50)

ax.grid()

plt.tight_layout()
plt.savefig("img/dn_tau_k2.pdf")
plt.show()

