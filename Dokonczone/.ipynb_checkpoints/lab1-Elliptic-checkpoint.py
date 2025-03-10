import numpy as np
import scipy as sc
import scipy.sparse.linalg

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from typing import Optional

#import sys
#np.set_printoptions(threshold=sys.maxsize)

class SetupElliptic:
    def __init__(self, dx: float, dy: Optional[float] = None):
        self.x_range = (0.0,1.0)
        self.y_range = (0.0,1.0)
        self.boundary_condition = 0.0

        dy = dy or dx

        # There are 2 hard problems in computer science:
        # cache invalidation, naming things, and off-by-1 errors.
        self.x_num = round((self.x_range[1] - self.x_range[0]) / dx) + 1
        self.y_num = round((self.y_range[1] - self.y_range[0]) / dy) + 1
        self.X, self.dx = np.linspace(*self.x_range, self.x_num, retstep=True)
        self.Y, self.dy = np.linspace(*self.y_range, self.y_num, retstep=True)

    @staticmethod
    def f(x, y):
        return (-2)*(np.pi ** 2)*np.sin(np.pi*x)*np.sin(np.pi*y)

    @staticmethod
    def exact(x, y):
        return np.sin(np.pi*x)*np.sin(np.pi*y)

def scheme_elliptic(setup):
    size = setup.x_num * setup.y_num
    print("size=", size)
    hx = 1.0/setup.dx ** 2
    hy = 1.0/setup.dy ** 2

    u = None
    A = np.zeros((size,size))
    np.fill_diagonal(A, 2*(hx + hy))

    #print("Enter first loop")
    for i in range(0, setup.y_num):
      for j in range(0,setup.x_num-1):
        index = i*setup.x_num + j
        A[index, index+1] = -hx
        A[index+1, index] = -hx
    for i in range(setup.x_num, size):
      A[i-setup.x_num, i] = -hy
      A[i, i-setup.x_num] = -hy
    print(A)
    
    vect_f = np.zeros(size)
    for i in range(0,setup.x_num):
      xi = i*setup.dx
      for j in range(0,setup.y_num):
        yj = j*setup.dy
        index = i*setup.y_num + j
        vect_f[index] = (-1)*setup.f(xi,yj)
    u = np.linalg.solve(A, vect_f)
    u = u.reshape(setup.x_num, setup.y_num)
    return u

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
steps = [.1, .05, .01, .005, .001]
times = []
errors = []
for h in steps:
    pass
    # TODO
"""