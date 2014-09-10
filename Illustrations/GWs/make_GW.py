import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import axes3d
import sys

N = 200

x = np.linspace(-15, 15, N)
y = np.linspace(-15, 15, N)
X, Y = np.meshgrid(x, y)
R = np.sqrt(X ** 2 + Y **2)
theta = np.arctan2(Y,X)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

R[R<5.6] = 0.0

if "1" in sys.argv:
    R[Y<-1.4] = 0
elif "2" in sys.argv:
    R[Y>=-1.4] = 0


Z = np.sin(1.2*R) / (0.2*R) +  0.5*np.cos(3 * theta + 0.4)


ax.plot_surface(X, Y, Z, rstride=4, cstride=4, 
                  color="cyan", alpha=0.3)
#ax.plot_wireframe(X, Y, Z, rstride=40, cstride=40, 
                  #color="k")

ax.set_zlim(-10, 10)

plt.axis('off')

if "1" in sys.argv:
    plt.savefig("GW1.pdf",bbox_inches='tight', transparent=True, dpi=400)
elif "2" in sys.argv:
    plt.savefig("GW2.pdf",bbox_inches='tight', transparent=True, dpi=400)
plt.show()
