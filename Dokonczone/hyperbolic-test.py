import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

x_range = (0.0, 1.0)
t_range = (0.0, 1.0)

x_num = round((x_range[1] - x_range[0]) / dx) + 1
t_num = round((t_range[1] - t_range[0]) / dt) + 1

X, dx = np.linspace(*x_range, x_num, retstep=True)
T, dt = np.linspace(*t_range, t_num, retstep=True)

