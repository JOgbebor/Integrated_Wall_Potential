# -*- coding: utf-8 -*-
"""
Created on Fri Apr 28 10:39:00 2023

@author: Jason Ogbebor
"""

"""

Fit lennard-jones interaction parameters for LAMMPS wall potentials such that 
the walls mimic a graphitic slit pore of prescribed width (W).

The goal is to simulate a confined fluid in LAMMPS, but to avoid including the 
wall atoms explicitly, so that the overall properties of the system are those 
of the confined fluid, not the fluid and walls.

"N" is the number of layers in the solid wall. Not all solid pores are layered.
For slit pores of amorphous solid, N = 1. But n (surface density) must change 
to reflect wall thickness as well.

"pad" determines how close to the walls to calculate the potential. It prevents
U from approaching infinity and making the plots difficult to read.

This script can be modified for different fluid-solid combinations by
changing the fluid-fluid, and solid-solid parameters. eps_sf and sig_sf assume
Lawrence-Berthelot mixing rules, but can be rewritten if the parameters are 
known explicitly.

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scipy.constants as c
from scipy.optimize import minimize

if True: # Plotting style guide

    style = {
        "figure.figsize": (10, 8),
        "font.size": 36,
        "axes.labelsize": "26",
        "axes.titlesize": "26",
        "xtick.labelsize": "26",
        "ytick.labelsize": "26",
        "legend.fontsize": "18",
        "xtick.direction": "in",
        "ytick.direction": "in",
        "lines.linewidth": 2,
        "lines.markersize": 12,
        "xtick.major.size": 18,
        "ytick.major.size": 18,
        "xtick.top": True,
        "ytick.right": True,
    }
    
    matplotlib.rcParams.update(style)
    markers = ["o", "X",  "^", "P", "d", "*", "s", ".", "x", ">"] * 5
    matplotlib.rcParams["font.family"] = ["serif"]

#============================================================================#

N = 6 # number of layers of graphene in each direction
# (ensure the fluid sees graphene until the cutoff distance)

W = 15 # A
d = W/2 # A
pad = 2.8

Points = 1000
Z = np.linspace((-d+pad), (d-pad), Points) # vertical positions relative to the center of the pore, which is z=0

n = 0.3765 # surface density of atoms in each layer (atoms / Ang^2)

eps_ss = 0.0556414 # kcal/mol
eps_ff = 0.2376681 # kcal/mol
eps_sf = np.sqrt(eps_ss*eps_ff) # kcal/mol
sigma_ss = 3.4 # Ang
sigma_ff = 3.4 # Ang
r_0 = (sigma_ss + sigma_ff) / 2 # sigma_sf, Ang

#============================================================================#

U_walls = []

U = [0] * Points

for nw in np.linspace(0,N-1,N):
    
    T1 = [(((r_0 / (d + (nw*3.35) + z))**10)) for z in Z]
    T2 = [(((r_0 / (d + (nw*3.35) - z))**10)) for z in Z]
    T3 = [(((r_0 / (d + (nw*3.35) + z))**4)) for z in Z]
    T4 = [(((r_0 / (d + (nw*3.35) - z))**4)) for z in Z]

    U_wall = [(4*c.pi*n*eps_sf*(r_0**2))*((1/5)*(t1+t2)-(1/2)*(t3+t4)) for t1,t2,t3,t4 in zip(T1,T2,T3,T4)]

    U_walls.append(U_wall)

    U = [u+u_w for u,u_w in zip(U,U_wall)]

#============================================================================#

fig = plt.figure()
plt.plot(Z, U, color="blue")
plt.ylabel("Integrated Potential")
plt.xlabel("z-distance")
#plt.xlim(-15,15)
#plt.ylim(-3.5,5)
#plt.savefig("file_name.pdf", bbox_inches="tight")
plt.show(fig)
plt.close(fig)

#============================================================================#

R_sq = [] # keep track of the fitting process

Z_fit = np.linspace(pad, (d+pad), int(Points/2)) # half of the pore size (only looking at one wall)

lj_type = "10-4-3" # type of wall potential to fit

def func(x):

    eps, sig = x

    if lj_type == "9-3":
        U_left = [eps*((2/15)*(sig/(d+z))**9 - (sig/(d+z))**3) for z in Z]
        U_right = [eps*((2/15)*(sig/(d-z))**9 - (sig/(d-z))**3) for z in Z]
        U_fit = [left+right for left,right in zip(U_left,U_right)]

    elif lj_type == "10-4-3":
        U_left = [2*c.pi*eps*((2/5)*(sig/(d+z))**10 - (sig/(d+z))**4 - ((np.sqrt(2)*sig**3) / (3*((d+z) + (0.61/np.sqrt(2))*sig)**3))) for z in Z]
        U_right = [2*c.pi*eps*((2/5)*(sig/(d-z))**10 - (sig/(d-z))**4 - ((np.sqrt(2)*sig**3) / (3*((d-z) + (0.61/np.sqrt(2))*sig)**3))) for z in Z]
        U_fit = [left+right for left,right in zip(U_left,U_right)]

    elif lj_type == "12-6":
        U_left = [4*eps*((sig/(d+z))**12 - (sig/(d+z))**6) for z in Z]
        U_right = [4*eps*((sig/(d-z))**12 - (sig/(d-z))**6) for z in Z]
        U_fit = [left+right for left,right in zip(U_left,U_right)]

    else:
        raise Exception("Unrecognized lj_type")

    r_sq = np.corrcoef(U_fit, U)[0, 1]
    R_sq.append(r_sq)

    return (1-r_sq)

x0 = [0.448, 3] # Initial guesses (play around with these values to see how it affects the fit)

# minimize() doesn't change epsilon x0[0] that much, so it must be manipulated manually

params = minimize(func, x0, tol=1e-6, bounds=((0.1,10),(1,5)))

eps = params.x[0]
sig = params.x[1]

if lj_type == "9-3":
    U_left = [eps*((2/15)*(sig/(d+z))**9 - (sig/(d+z))**3) for z in Z]
    U_right = [eps*((2/15)*(sig/(d-z))**9 - (sig/(d-z))**3) for z in Z]
    U_fit = [left+right for left,right in zip(U_left,U_right)]

if lj_type == "10-4-3":
    U_left = [2*c.pi*eps*((2/5)*(sig/(d+z))**10 - (sig/(d+z))**4 - ((np.sqrt(2)*sig**3) / (3*((d+z) + (0.61/np.sqrt(2))*sig)**3))) for z in Z]
    U_right = [2*c.pi*eps*((2/5)*(sig/(d-z))**10 - (sig/(d-z))**4 - ((np.sqrt(2)*sig**3) / (3*((d-z) + (0.61/np.sqrt(2))*sig)**3))) for z in Z]
    U_fit = [left+right for left,right in zip(U_left,U_right)]

if lj_type == "12-6":
    U_left = [4*eps*((sig/(d+z))**12 - (sig/(d+z))**6) for z in Z]
    U_right = [4*eps*((sig/(d-z))**12 - (sig/(d-z))**6) for z in Z]
    U_fit = [left+right for left,right in zip(U_left,U_right)]

#============================================================================#

epsilon = round(eps,4)
sigma = round(sig,3)
rsq = round(R_sq[-1],4)

if W%10 == 0:
    H = int(W/10)
else:
    H = W/10

fig = plt.figure()
plt.plot(Z, U, color="blue", label="Integrated Potential")
plt.plot(Z, U_fit, color="red", linestyle="--", label="Fit")
plt.ylabel("Potential Energy")
plt.xlabel("z-distance")
plt.title(rf"LJ Type: {lj_type}, W = {W} $\AA$", x=0.305, y=1.08)
plt.title(rf"$\epsilon = ${epsilon} kcal/mol, $\sigma =${sigma} $\AA$, R$^2 =${rsq}", loc="left")
#plt.xlim(-23,-15)
#plt.ylim(-3.5,5)
plt.legend(loc="best")
#plt.savefig(f"C:/Argon_Carbon/{H}nm/Integrated_Potential.pdf", bbox_inches="tight")
plt.show(fig)
plt.close(fig)
