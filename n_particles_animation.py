# -*- coding: utf-8 -*-
"""
N PARTICLE ANIMATION (Requare ffmpeg!)

Creates N (line 46) frames for the animation from the file "results_R.txt". 
The number of colors (line 41) has to be the same as the number of bodies. 
The limits of the figure (line 36-38 and 63-65) should be changed for every different animation. 
"""

from os import system
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

system("mkdir tmp")

data = np.loadtxt("results_R.txt", unpack=True)

mpl.rcParams.update({'font.size': 6})
mpl.rcParams['legend.fontsize'] = 8
plt.ioff()  # does not open all the figures
fig = plt.figure(figsize=(8, 4.5))
ax1 = fig.gca(projection='3d')
ax1.set_xlabel('X axis')
ax1.set_ylabel('Y axis')
ax1.set_zlabel('Z axis')
ax1.cla()

# The number of colors has to be the same as the number of bodies
colors = ["red", "orange", "blue"]

N = 500  # number of frames

counter = 0
scale = int(len(data[0])/N)

while counter < N:
    if counter <= 50:
        initial = 0
    else:
        initial = (counter - 50)*scale
        data[3*i], data[3*i+1], data[3*i+2]
    for i in range(len(colors)):
        ax1.scatter(data[3*i, counter*scale], data[3*i+1, counter *
                    scale], data[3*i+2, counter*scale], color=colors[i])
    for i in range(len(colors)):
        ax1.plot(data[3*i, initial:counter*scale], data[3*i+1, initial:counter *
                 scale], data[3*i+2, initial:counter*scale], color=colors[i])
    fig.savefig("tmp/tmp%04d.png" % (counter),
                format="png", dpi=240.)  # 6.4 4.8
    ax1.cla()
    counter += 1

system("ffmpeg -framerate 30 -pattern_type glob -i 'tmp/*.png' -c:a copy -shortest -c:v libx264 -pix_fmt yuv420p out.mp4")
