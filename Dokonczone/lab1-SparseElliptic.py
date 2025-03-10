import numpy as np
import scipy.sparse.linalg

#from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
#import matplotlib.pylab as plab
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

class SetupElliptic:
    def __init__(self, dx: float, dy: float, f, x_range, y_range, boundary_right, boundary_up, boundary_left, boundary_down):
        self.x_range = x_range
        self.y_range = y_range
        self.boundary_left = boundary_left
        self.boundary_right = boundary_right
        self.boundary_up = boundary_up
        self.boundary_down = boundary_down        

        self.f = f
        self.x_num = round((self.x_range[1] - self.x_range[0]) / dx) + 1
        self.y_num = round((self.y_range[1] - self.y_range[0]) / dy) + 1
        self.X, self.dx = np.linspace(*self.x_range, self.x_num, retstep=True)
        self.Y, self.dy = np.linspace(*self.y_range, self.y_num, retstep=True)

def scheme_elliptic(setup):
    ex = np.ones(setup.x_num-2)
    A = scipy.sparse.dia_matrix((np.array([-1/setup.dx**2 * ex, 2*(1/setup.dx**2 + 1/setup.dy**2) * ex, -1/setup.dx**2 * ex]), [-1, 0, 1]), shape=(setup.x_num-2, setup.x_num-2))
    B = scipy.sparse.dia_matrix((-1/setup.dy**2 * ex, [0]), shape=(setup.x_num-2, setup.x_num-2))
    X_data = [ [None]*(setup.y_num-2) for i in range(setup.y_num-2)]
    X_data[0][0] = A
    for i in range(1,setup.y_num-2):
        X_data[i][i] = A
        X_data[i-1][i] = B
        X_data[i][i-1] = B
    X = scipy.sparse.bmat(X_data, format="csc")
    f_vec = np.zeros((setup.x_num-2) * (setup.y_num-2))
    n = 0
    for j in range(setup.y_num-2):
        for i in range(setup.x_num-2):
            f_vec[n] += setup.f(setup.X[i+1], setup.Y[j+1])
            if (i == 0):
                f_vec[n] += 1/setup.dx**2 * setup.boundary_left(setup.Y[j+1])
            if (i == setup.x_num-3):
                f_vec[n] += 1/setup.dx**2 * setup.boundary_right(setup.Y[j+1])
            if (j == 0):
                f_vec[n] += 1/setup.dy**2 * setup.boundary_down(setup.X[i+1])
            if (j == setup.y_num-3):
                f_vec[n] += 1/setup.dy**2 * setup.boundary_up(setup.X[i+1])
            n+=1

    u = scipy.sparse.linalg.spsolve(X, f_vec).reshape(setup.y_num-2, setup.x_num-2) # od 1 do Nx-1, 1 do Ny-1
    r = np.zeros((setup.y_num, setup.x_num))
    for j in range(setup.y_num):
        r[j,0] = setup.boundary_left(setup.Y[j])
        r[j,setup.x_num-1] = setup.boundary_right(setup.Y[j])
    for i in range(setup.x_num):
        r[0,i] = setup.boundary_down(setup.X[i])
        r[setup.y_num-1,i] = setup.boundary_up(setup.X[i])
    r[1:setup.y_num-1, 1:setup.x_num-1] = u
    return r

def plot_surface(u, setup, i = 1):
    fig = plt.figure(i)
    ax = plt.axes(projection='3d')
    X, Y = np.meshgrid(setup.X, setup.Y)
    surf = ax.plot_surface(X, Y, u, cmap=cm.coolwarm,
                           linewidth=0, antialiased=False)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    fig.colorbar(surf, shrink=0.5, aspect=5)


setup_elliptic = SetupElliptic(0.01, 0.01,
    lambda x,y: (2)*(np.pi ** 2)*np.sin(np.pi*x)*np.sin(np.pi*y),
    (0.0, 1.0), (0.0, 1.0),
    lambda y: 0.0, # boundary right
    lambda x: 0.0, # boundary up
    lambda y: 0.0, # boundary left
    lambda x: 0.0  # boundary down
)
numerical_elliptic = scheme_elliptic(setup_elliptic)


# exact solution
#X, Y = np.meshgrid(setup_elliptic.X, setup_elliptic.Y)
#exact_elliptic = setup_elliptic.exact(X, Y)
#print("exact=", exact_elliptic)


plot_surface(numerical_elliptic, setup_elliptic, 1)
#plot_surface(exact_elliptic, setup_elliptic, 2)
#plot_surface(numerical_elliptic - exact_elliptic, setup_elliptic, 3)
plt.show()



"""
steps = [.1, .05, .01, .005, .001]
times = []
errors = []
for h in steps:
    pass
    # TODO
"""