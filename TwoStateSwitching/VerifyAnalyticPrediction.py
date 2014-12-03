import numpy as np
import matplotlib.pyplot as plt
from LyneTools import GenerateModel, MagnitudeOfPeaks

R = 0.5
upsilon = 0.3
nu_dot_1 = 1e-12

time, tres, T = GenerateModel(R=R, nu_dot_1=nu_dot_1, upsilon=upsilon)

print "From data : DeltaPhiMax={}".format(np.max(tres))


DeltaPhiMaxAnalytic= np.pi/32.0 * upsilon * T**2 * nu_dot_1 
DeltaTMaxAnalytic = DeltaPhiMaxAnalytic / (2 * np.pi)
print "Analytic {}".format(DeltaTMaxAnalytic)

plt.plot(time, tres)
plt.show()



