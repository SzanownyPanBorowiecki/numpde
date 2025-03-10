# Konończyć dla schematu jawnego i niejawnego!
# Następne zajęcia: zaczniemy robić schemat wariacyjny

import numpy as np
import scipy as sc
import scipy.sparse.linalg
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from scipy.sparse import lil_matrix, csr_matrix
from typing import Optional
import time

class SetupHyperbolic:
	def __init__(self, 
				 dx: float, 
				 dt: Optional[float] = None,
				 output_freq: Optional[int] = 10,
				 alpha: float = 0.5
				 ):
		self.x_range = (0.0, 1.0)
		self.t_range = (0.0, 2.0)
		self.boundary_condition = (lambda t: 0.0, lambda t: 0.0)

		self.alpha = alpha 
		dt = dt or dx ** 2
		# There are 2 hard problems in computer science:
		# cache invalidation, naming things, and off-by-1 errors.
		self.x_num = round((self.x_range[1] - self.x_range[0]) / dx) + 1
		self.t_num = round((self.t_range[1] - self.t_range[0]) / dt) + 1

		self.X, self.dx = np.linspace(*self.x_range, self.x_num, retstep=True)
		self.T, self.dt = np.linspace(*self.t_range, self.t_num, retstep=True)
		
		self.output_freq = output_freq

	@staticmethod
	def initial(x):
		return 0.125 * np.sin(np.pi * x)

	@staticmethod
	def initial_der(x):
		return 0

	@staticmethod
	def exact(t, x):
		return 0.125 * np.sin(np.pi * x) * np.cos(np.pi * t)

def create_matrices(setup):
	u0 = np.empty(setup.x_num)
	u0t = np.empty(setup.x_num)
	for i in range(setup.x_num):
		u0[i] = setup.initial(setup.x_range[0]+i*setup.dx)
		u0t[i] = setup.initial_der(setup.x_range[0]+i*setup.dx)

	#A = None
	Lh = np.zeros((setup.x_num-2, setup.x_num-2))
	np.fill_diagonal(Lh, -2)

	#print(setup.x_num-2)
	for i in range(1, setup.x_num-2):
		Lh[i-1,i] = 1
		Lh[i,i-1] = 1

	Lh /= (setup.dx ** 2)
	print(Lh)
	I = np.identity(setup.x_num-2)
	X = (1/setup.dt ** 2) * I - setup.alpha * Lh
	A = 2*I + np.linalg.inv(X)*Lh
	print("A=")
	print(A)
	# TODO
	return A, u0, u0t


def create_RH(u, u_prev, setup):
	d = np.zeros(setup.x_num - 2)
	r = setup.dt ** 2 / setup.dx ** 2

	# TODO

	return d


def fictional(x, setup):
	pass


def scheme_hyperbolic(setup):
	u_matrix = np.zeros((setup.x_num, setup.t_num // setup.output_freq))
	A, u0, u0t = create_matrices(setup)
	#u_prev = np.zeros_like(u)
	#for i, xi in enumerate(setup.X):
	#	u_prev[i] = fictional(xi, setup)

	u = np.empty(setup.x_num)
	u1 = np.empty(setup.x_num)
	u1[0] = u0[0] + (setup.dt/(2*(setup.dx ** 2)))*(-2*u0[0]+u0[1])
	u1[setup.x_num-1] = u0[setup.x_num-1] + (setup.dt/(2*(setup.dx ** 2)))*(u0[setup.x_num-2]-2*u0[setup.x_num-1])
	for i in range(1,setup.x_num-1):
		u1[i] = u0[i] + (setup.dt/(2*(setup.dx ** 2)))*(u0[i-1]-2*u0[i]+u0[i+1])

	for t in range(setup.t_num - 1):
		#d = create_RH(u, u_prev, setup)

		u[1:-1] = A @ u1[1:-1] - u0[1:-1]
		if (t-1) % setup.output_freq == 0:
			u_matrix[:, (t-1) // setup.output_freq] = u[:]
		u0 = u1
		u1 = u
	print(u_matrix)
	return u_matrix

#%%time
# numerical solution
#setup_hyperbolic = SetupHyperbolic(0.0001, 0.001)
setup_hyperbolic = SetupHyperbolic(0.1, 0.00001, alpha=0.5)
numerical_hyperbolic = scheme_hyperbolic(setup_hyperbolic)

def plot_surface(u, setup, i = 1, title=None):
	fig = plt.figure(i)
	ax = plt.axes(projection='3d')
	if (title != None):
		ax.set_title(title)
	# Prepare grid.
	T = np.linspace(*setup.t_range, setup.t_num // setup.output_freq)
	X = np.linspace(*setup.x_range, setup.x_num)
	T, X = np.meshgrid(T, X)

	# Plot the surface.
	surf = ax.plot_surface(T, X, u, cmap=cm.coolwarm,
						   linewidth=0, antialiased=False)

	# Customize the z axis.
	ax.zaxis.set_major_locator(LinearLocator(10))
	ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
	# Add a color bar which maps values to colors.
	fig.colorbar(surf, shrink=0.5, aspect=5)


# exact solution
X = np.linspace(*setup_hyperbolic.x_range, setup_hyperbolic.x_num)
T = np.linspace(*setup_hyperbolic.t_range, setup_hyperbolic.t_num // setup_hyperbolic.output_freq)
X, T = np.meshgrid(X, T)
exact_hyperbolic = setup_hyperbolic.exact(T, X).T

plot_surface(numerical_hyperbolic, setup_hyperbolic, 1)
#plot_surface(exact_hyperbolic, setup_hyperbolic, 2)
#plot_surface(numerical_hyperbolic - exact_hyperbolic, setup_hyperbolic, 3)
plt.show()
