# -*- coding: utf-8 -*-
"""
N PARTICLES PLOT

Plots the trajectory of all the particles in the file "results_R.txt".

"""

import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt


data1 = np.loadtxt("result/result_9_9_R.txt", unpack=True)

mpl.rcParams['legend.fontsize'] = 8

fig = plt.figure()
ax1 = fig.gca(projection='3d')
ax1.set_xlabel('X')
ax1.set_ylabel('Y')
ax1.set_zlabel('Z')

for i in range(int(len(data1)/3)):
    ax1.plot(data1[3*i], data1[3*i+1], data1[3*i+2], label=str(i))

plt.legend(["Alpha Cantauri AB", "Proxima Cantauri", "Planet"])
plt.show()
