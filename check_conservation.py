# -*- coding: utf-8 -*-
"""
ESTIMATION OF ERRORS

This program needs have the files "results_R.txt" and "results_V.txt" in
its folder. 

To estimate the error of a simulation, run: 
>>> ERROR()

If there is the possibility of two bodies colliding or getting very close, run:
>>> PLOT_ERROR()
to see if there is a significant change in the error of both energy and angular momentum. 

"""


import numpy as np
import matplotlib.pyplot as plt

def ERROR():
    """
    Determines the error using:
    Error X (%) = 100*abs(Xfinal - Xinitial)/Xinitial
    Returns error in %. 
    """
    m, r, v = get_data()
    Etotal_i = E(m, r[0], v[0])
    Etotal_f = E(m, r[-1], v[-1])
    LT_i, Lx_i, Ly_i, Lz_i = L(m, r[0], v[0])
    LT_f, Lx_f, Ly_f, Lz_f = L(m, r[-1], v[-1])
    
    error_E = abs(Etotal_i-Etotal_f)/abs(Etotal_i)*100
    error_LT = abs(LT_i-LT_f)/LT_i*100  
    error_L = np.array([abs(Lx_i-Lx_f)/Lx_i*100, abs(Ly_i-Ly_f)/Ly_i*100, abs(Lz_i-Lz_f)/Lz_i*100])
    
    print("\nThe error of the final state compared to the initial one is:")
    print("\nError X (%) = 100*abs(Xfinal - Xinitial)/Xinitial")
    print("\nError E (%) = {0:0.15f}".format(error_E))
    print("\nError L (%) = {0:0.15f}".format(error_LT))
    
    return error_E, error_LT, error_L
    
    
def PLOT_ERROR():
    """
    Plots the evolution of the error vs iteration. 
    Error X (%) = 100*abs(Xfinal - Xinitial)/Xinitial
    """
    m, r, v = get_data()
    Etotal = np.array([E(m, r[i], v[i]) for i in range(len(v))])
    error_E = abs(Etotal - Etotal[0])/Etotal[0]*100
    LT = []
    Lx = []
    Ly = []
    Lz = []
    for i in range(len(r)):
        LTi, Lxi, Lyi, Lzi = L(m, r[i], v[i])
        LT += [LTi]
    LT = np.array(LT)
    error_L = abs(LT - LT[0])/LT[0]*100

    t = range(len(Etotal))
    fig1 = plt.figure(1, figsize=[14,7])
    ax11 = fig1.add_subplot(121)
    ax11.plot(t, error_E, ".")
    ax11.set_xlabel("iteration [#]")
    ax11.set_ylabel("error in energy [%]")
    
    ax12 = fig1.add_subplot(122)
    ax12.plot(t, error_L, ".")
    ax12.set_xlabel("iteration [#]")
    ax12.set_ylabel("error in modulus of angular momentum [%]")

    fig1.tight_layout()

    return error_E, error_L


def E(m, r, v):
    """
    Calculates total energy of a given system in a given instant. 
    All data has to be normalized. 
    E = T + V
    v : vx, vy, vz X n bodies (vector of 3N scalars)
    m : m X n bodies (vector of N scalars)
    """

    Ec = sum([0.5*m[i]*(v[3*i  ]**2 + v[3*i+1]**2 + v[3*i+2]**2) for i in range(len(m))])
    Ep = 0
    for i in np.arange(0, len(m) - 1):
        for j in np.arange(i + 1, len(m)):
            Ep += m[i]*m[j]/np.sqrt((r[3*i  ] - r[3*j  ])**2 + \
                                    (r[3*i+1] - r[3*j+1])**2 + \
                                    (r[3*i+2] - r[3*j+2])**2 )
    ET = Ec - Ep
    return ET
    
    
def L(m, r, v):
    """
    Calculates total angular momentum of a given system in a given instant. 
    All data has to be normalized. 
    v : vx, vy, vz X n bodies (vector of 3N scalars)
    r : rx, ry, rz X n bodies (vector of 3N scalars)
    m : m X n bodies (vector of N scalars)
    """
    Lx = sum([m[i]*(r[3*i+1]*v[3*i+2] - r[3*i+2]*v[3*i+1]) for i in range(len(m))])
    Ly = sum([m[i]*(r[3*i+2]*v[3*i  ] - r[3*i  ]*v[3*i+2]) for i in range(len(m))])
    Lz = sum([m[i]*(r[3*i  ]*v[3*i+1] - r[3*i+1]*v[3*i  ]) for i in range(len(m))])
    LT = np.sqrt(Lx**2 + Ly**2 + Lz**2)
    
    return LT, Lx, Ly, Lz
    
    
def get_data():
    """
    Gets all data from "initial_data.csv", "results_R.txt" and "results_V.txt". 
    The output are 3 vectors: r, v, m
    """
    #LOADS INITIAL DATA
    RVM = np.loadtxt("initial_data.csv", unpack=True, delimiter=",")    
    m = np.array(RVM[-1])
    RVM = np.transpose(RVM)
    r0 = []
    for i in range(len(m)):
        r0 += RVM[i, 0:3].tolist()
    r0 = np.array(r0)
    v0 = []
    for i in range(len(m)):
        v0 += RVM[i, 3:6].tolist()
    v0 = np.array(v0)
    #LOADS R
    f = open("results_R.txt", "r")
    r = f.read()
    f.close()
    r = np.array([[float(j) for j in i.split(" ")] for i in r.split("\n")[:-1]])
    r = np.array([r0.tolist()] + r.tolist())
    #LOADS V
    f = open("results_V.txt", "r")
    v = f.read()
    f.close()
    v = np.array([[float(j) for j in i.split(" ")] for i in v.split("\n")[:-1]])
    v = np.array([v0.tolist()] + v.tolist())
    
    return m, r, v
