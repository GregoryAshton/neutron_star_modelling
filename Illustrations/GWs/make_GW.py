import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import axes3d
import sys

if "1" in sys.argv:
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    x = np.linspace(-15, 15, 500)
    y = np.linspace(-15, 15, 500)
    X, Y = np.meshgrid(x, y)

    R = np.sqrt(X ** 2 + Y **2)
    theta = np.arctan2(Y,X)

    R[R<5.5] = 0
    R[Y<0] = 0

    Z = np.sin(2*R) / (0.2*R) +  0.5*np.cos(3 * theta + 0.4)

    ax.plot_wireframe(X, Y, Z, rstride=10, cstride=10)
    ax.set_zlim(-10, 10)

    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    plt.axis('off')
    plt.savefig("GW1.pdf",bbox_inches='tight')

    plt.show()

if "2" in sys.argv:
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    x = np.linspace(-15, 15, 500)
    y = np.linspace(-15, 15, 500)
    X, Y = np.meshgrid(x, y)

    R = np.sqrt(X ** 2 + Y **2)
    theta = np.arctan2(Y,X)

    R[R<5.5] = 0
    R[Y>0] = 0

    Z = np.sin(2*R) / (0.2*R) +  0.5*np.cos(3 * theta + 0.4)

    ax.plot_wireframe(X, Y, Z, rstride=10, cstride=10)
    ax.set_zlim(-10, 10)

    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    plt.axis('off')
    plt.savefig("GW2.pdf",bbox_inches='tight', transparent=True)

    plt.show()
