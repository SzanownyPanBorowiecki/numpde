# u_t = u_xx                 t > 0
# u(t,0) = u(t,1) = 0        t > 0
# u(0,x) = sin(pix)           x \in [0,1]
# 
import numpy as np
import scipy as sc
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from typing import Optional

class SetupParabolic:
    def __init__(self,
                 dx: float,
                 dt: Optional[float] = None,
                 output_freq: int = 100,
                 alpha: float = 0
                 ):
        self.x_range = (0.0,1.0)
        self.t_range = (0.0,1.0)
        self.boundary_condition = 0

        self.alpha = alpha  # just ignore it
        #dt = dt or dx ** 2
        dt = dt or 10*dx ** 2
        # There are 2 hard problems in computer science:
        # cache invalidation, naming things, and off-by-1 errors.
        self.x_num = round((self.x_range[1] - self.x_range[0]) / dx) + 1
        self.t_num = round((self.t_range[1] - self.t_range[0]) / dt) + 1

        self.X, self.dx = np.linspace(*self.x_range, self.x_num, retstep=True)
        self.T, self.dt = np.linspace(*self.t_range, self.t_num, retstep=True)

        self.output_freq = output_freq

    @staticmethod
    def initial(x):
        return np.sin(np.pi*x)

    @staticmethod
    def exact(t, x):
        return np.exp(-(np.pi ** 2)*t)*np.sin(np.pi*x)

def create_matrices(setup):
    u = np.empty(setup.x_num)

    # TODO initial condition
    for i in range(0,setup.x_num):
        u[i] = setup.initial(setup.x_range[0]+i*setup.dx)
        

    # TODO boundary condition

    # TODO Matrix A
    #s = setup.dt / (setup.dx ** 2)
    s = 1/setup.dx ** 2
    Lh = np.zeros((setup.x_num-2, setup.x_num-2))
    #np.fill_diagonal(Lh, 1-2*s)
    np.fill_diagonal(Lh, -2)
    for i in range(1, setup.x_num-2):
        Lh[i-1,i] = 1
        Lh[i,i-1] = 1
    Lh = s*Lh
    I = np.zeros((setup.x_num-2, setup.x_num-2))
    X = np.linalg.inv(1/setup.dt*I - setup.alpha*Lh)*(1/setup.dt*I + (1-setup.alpha)*Lh)
    print("A=", A)
    return A, u

def create_RH(u, setup):
    d = np.zeros(setup.x_num - 2)
    # TODO right hand side
    return d


def scheme_parabolic(setup):
    u_matrix = np.zeros((setup.x_num, setup.t_num // setup.output_freq))
    A, u = create_matrices(setup)
    for t in range(1, setup.t_num):  # suboptimal, but simple
        u[1:-1] = np.matmul(A, u[1:-1])
        if (t-1) % setup.output_freq == 0:
            u_matrix[:, (t-1) // setup.output_freq] = u[:]

    return u_matrix
    #return u

def plot_surface(u, setup, i = 1, title=None):
    fig = plt.figure(i)
    ax = plt.axes(projection='3d')
    if (title != None):
        ax.set_title(title)
    # Prepare grid.
    T = np.linspace(*setup_parabolic.t_range, setup_parabolic.t_num // setup_parabolic.output_freq)
    X = np.linspace(*setup_parabolic.x_range, setup_parabolic.x_num)
    T, X = np.meshgrid(T, X)

    # Plot the surface.
    surf = ax.plot_surface(T, X, u, cmap=cm.coolwarm,
                           linewidth=0, antialiased=False)

    # Customize the z axis.
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.5, aspect=5)


# %%time
# numerical solution
setup_parabolic = SetupParabolic(0.01,0.5*(0.01 ** 2),100)
numerical_parabolic = scheme_parabolic(setup_parabolic)


print(numerical_parabolic)

# exact solution
X = np.linspace(*setup_parabolic.x_range, setup_parabolic.x_num)
T = np.linspace(*setup_parabolic.t_range, setup_parabolic.t_num // setup_parabolic.output_freq)
X, T = np.meshgrid(X, T)
exact_parabolic = setup_parabolic.exact(T, X).T
#print("exact=", exact_parabolic)

plot_surface(numerical_parabolic, setup_parabolic,1, "numerical")
plot_surface(exact_parabolic, setup_parabolic,2, "exact")
plot_surface(numerical_parabolic - exact_parabolic, setup_parabolic,3, "numerical - exact")



"""----------------------
        Animacja
----------------------"""
#%matplotlib inline
#import numpy as np
#import matplotlib.pyplot as plt
import seaborn as sns

# create a figure and axes
fig = plt.figure(figsize=(6,5))
ax1 = plt.subplot(1,1,1)   
# ax2 = plt.subplot(1,2,2)

# set up the subplots as needed
ax1.set_xlim(setup_parabolic.x_range)            
ax1.set_ylim((0, 1))
ax1.set_xlabel('X')
ax1.set_ylabel('Magnitude')

# create objects that will change in the animation. These are
# initially empty, and will be given new values for each frame
# in the animation.
txt_title = ax1.set_title('')
line1, = ax1.plot([], [], 'orange', lw=2)     # ax.plot returns a list of 2D line objects

# animation function. This is called sequentially
def drawframe(n):
    x = setup_parabolic.X
    y = numerical_parabolic[:, n]
    line1.set_data(x, y)
    txt_title.set_text('Frame = {0:4d}'.format(n))
    return (line1,)

from matplotlib import animation

# blit=True re-draws only the parts that have changed.
anim = animation.FuncAnimation(fig, drawframe, frames=numerical_parabolic.shape[1], interval=20)

plt.show()
