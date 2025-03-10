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

from pprint import pprint

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
		#return 0
		return 0.125 * np.sin(np.pi * x)

	@staticmethod
	def initial_der(x):
		return 0

	@staticmethod
	def exact(t, x):
		return 0.125 * np.sin(np.pi * x) * np.cos(np.pi * t)

def create_matrices(setup):
	u = np.empty(setup.x_num)
	A = None
	# TODO
	return A, u


def create_RH(u, u_prev, setup):
	d = np.zeros(setup.x_num - 2)
	r = setup.dt ** 2 / setup.dx ** 2

	# TODO

	return d


def fictional(x, setup):
	pass


def scheme_hyperbolic(setup):
	N = setup.x_num
	T = setup.t_num

	dx = setup.dx
	dt = setup.dt
	
	c = (dt/dx)**2
	a = setup.alpha

	u_matrix = np.zeros((N, T // setup.output_freq))
	p = np.zeros(N)
	q = np.zeros(N)
	for i, xi in enumerate(setup.X):
		p[i] = setup.initial(xi)
		q[i] = setup.initial_der(xi)

	if a == 0:
		v0 = np.zeros(N)
		v0[0] = setup.boundary_condition[0](0)
		v0[N-1] = setup.boundary_condition[1](0)
		for i in range(1,N-2):
			v0[i] = p[i] - dt*q[i] + (c/2)*(p[i-1] - 2*p[i] + p[i+1])

		v1 = p
		u_matrix[:, 0] = v1[:]
		for t in range(1,T-1):
			tn = setup.T[t]
			v2 = np.zeros(N)
			v2[0] = setup.boundary_condition[0](tn)
			v2[N-1] = setup.boundary_condition[1](tn)
			for i in range(1,N-2):
				v2[i] = 2*v1[i] - v0[i] + c*(v1[i-1] - 2*v1[i] + v1[i+1])
			if t % setup.output_freq == 0:
				u_matrix[:, t // setup.output_freq] = v2[:]
			v0 = v1
			v1 = v2

	if a != 0:
		"""
		v0 = np.zeros(N)
		v0[0] = setup.boundary_condition[0](0)
		v0[1] = \
			(2+1/(c*a))*v0[0] \
			+ (1/a - 2 - 1/(c*a))*p[0] \
			- (1/(2*a)-1)*p[1] \
			+ (1/(c*a)+2)*q[0] \
			- q[1]
		v0[N-1] = setup.boundary_condition[1](0)
		for i in range(1,N-2):
			v0[i+1] = (2+1/(c*a))*v0[i] \
			+ v0[i-1] \
			+ (1/a - 2 - 1/(c*a))*p[i] \
			+ (1 - 1/(2*a))*(p[i-1] + p[i+1]) \
			+ (2 + 1/(c*a))*q[i] \
			- (q[i-1] + q[i+1])

		v1 = p
		print("v0=", v0)
		print("v1=", v1)

		u_matrix[:, 0] = v1[:]
		for t in range(1, T-1):
			tn = setup.T[t]
			v2 = np.zeros(N)
			v2[0] = setup.boundary_condition[0](tn)
			v2[1] = 
			"""

			"""
			v2[1] = \
				2*v2[0] \
				+ (1/(c*a))*(v2[0]-2*v1[0]+v0[0]) \
				- (1/a - 2)*(-2*v1[0] + v1[1]) \
				- (-2*v0[0] + v0[1])
			"""

			"""
			v2[N-1] = setup.boundary_condition[1](tn)
			for i in range(2, N-3):
				v2[i+1] = \
				2*v2[i] - v2[i-1] \
				+ (1/(c*a))*(v2[i] - 2*v1[i] + v0[i]) \
				- (1/a - 2)*(v1[i-1] - 2*v1[i] + v1[i+1]) \
				- (v0[i-1] - 2*v0[i] + v0[i+1])
				print(v2[i+1])
			return
			if t % setup.output_freq == 0:
				u_matrix[:, t // setup.output_freq] = v2[:]
			v0 = v1
			v1 = v2
		"""
	return u_matrix

def plot_surface(u, setup, i = 1, title=None):
	fig = plt.figure(i)
	ax = plt.axes(projection='3d')
	if (title != None):
		ax.set_title(title)
	
	T = np.linspace(*setup.t_range, setup.t_num // setup.output_freq)
	X = np.linspace(*setup.x_range, setup.x_num)

	T, X = np.meshgrid(T,X)
	surf = ax.plot_surface(T, X, u, cmap=cm.coolwarm,
						   linewidth=0, antialiased=False)

	ax.zaxis.set_major_locator(LinearLocator(10))
	ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
	fig.colorbar(surf, shrink=0.5, aspect=5)



#%%time
# numerical solution
setup_hyperbolic = SetupHyperbolic(0.001, 0.001, alpha=0)
numerical_hyperbolic = scheme_hyperbolic(setup_hyperbolic)
plot_surface(numerical_hyperbolic, setup_hyperbolic, 1, "numerical")

# exact solution
X = np.linspace(*setup_hyperbolic.x_range, setup_hyperbolic.x_num)
T = np.linspace(*setup_hyperbolic.t_range, setup_hyperbolic.t_num // setup_hyperbolic.output_freq)
X, T = np.meshgrid(X, T)
exact_hyperbolic = setup_hyperbolic.exact(T, X).T
plot_surface(exact_hyperbolic, setup_hyperbolic, 2, "exact")


plt.show()
