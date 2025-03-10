import numpy as np
import scipy as sc
import scipy.sparse.linalg

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from typing import Optional

class SetupElliptic:
    def __init__(self, dx: float, dy: float, lx: float, ly: float, v0, f):
        self.lx = lx
        self.ly = ly
        self.v0 = v0
        self.f = f
        dy = dy or dx
        ly = ly or lx

        self.nx = round(self.lx / dx)
        self.ny = round(self.ly / dy) 

        self.x_num = round(self.lx / dx) + 1
        self.y_num = round(self.ly / dy) + 1
        self.X, self.dx = np.linspace(*(0, lx), self.x_num, retstep=True)
        self.Y, self.dy = np.linspace(*(0, ly), self.y_num, retstep=True)

def eigenvector(setup, i, j):
    return lambda x,y: (2/np.sqrt(setup.lx*setup.ly))*np.sin(i*np.pi*(x/setup.lx))*np.sin(j*np.pi*(y/setup.ly))

def eigenvalue(setup, i, j):
    return 4*((np.sin(i*np.pi/(2*setup.nx))/setup.dx)**2 + (np.sin(j*np.pi/(2*setup.ny))/setup.dy)**2)

def inner_product(setup, u, v):
    s = 0;
    for x in setup.X:
        for y in setup.Y:
            s += u(x,y)*v(x,y)
    return setup.dx*setup.dy*s

def calc_u(coefs, x, y):
    r = 0
    for i in range(1,setup.x_num+1):
        for j in range(1,setup.y_num+1):
            r += coefs[i-1,j-1]*eigenvector(setup, i, j)(x, y)
    return r

def scheme_elliptic(setup):
    size = setup.x_num * setup.y_num
    coefs = np.zeros((setup.x_num, setup.y_num))
    for i in range(1,setup.x_num+1):
        for j in range(1,setup.y_num+1):
            coefs[i-1,j-1] = -inner_product(setup, setup.f, eigenvector(setup, i, j))/eigenvalue(setup, i, j)
           #print(coefs[i,j])
    print("coefs=", coefs)

    u = np.zeros((setup.x_num, setup.y_num))
    for i in range(setup.x_num):
        for j in range(setup.y_num):
            xi = i/setup.lx
            yj = j/setup.ly
            u[i,j] = calc_u(coefs, xi, yj)
            print(xi, yj, u[i,j])
    return u

setup = SetupElliptic(0.01, 0.1, 1.0, 1.0, 0, lambda x, y: (-2)*(np.pi ** 2)*np.sin(np.pi*x)*np.sin(np.pi*y))
u = scheme_elliptic(setup)
print(u)

"""
def plot_surface(u, setup, i = 1):
    fig = plt.figure(i)
    #ax = fig.gca(projection='3d')
    ax = plt.axes(projection='3d')

    # Prepare grid.
    X, Y = np.meshgrid(setup.X, setup.Y)

    # Plot the surface.
    surf = ax.plot_surface(X, Y, u, cmap=cm.coolwarm,
                           linewidth=0, antialiased=False)

    # Customize the z axis.
    # ax.set_zlim(-.2, 1.2)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.5, aspect=5)

    #plt.show()


#%time
# numerical solution
setup_elliptic = SetupElliptic(0.01)
numerical_elliptic = scheme_elliptic(setup_elliptic)
#print("numerical=", numerical_elliptic)


# exact solution
X, Y = np.meshgrid(setup_elliptic.X, setup_elliptic.Y)
exact_elliptic = setup_elliptic.exact(X, Y)
#print("exact=", exact_elliptic)


plot_surface(numerical_elliptic, setup_elliptic, 1)
plot_surface(exact_elliptic, setup_elliptic, 2)
plot_surface(numerical_elliptic - exact_elliptic, setup_elliptic, 3)
plt.show()
"""

"""
steps = [.1, .05, .01, .005, .001]
times = []
errors = []
for h in steps:
    pass
    # TODO
"""