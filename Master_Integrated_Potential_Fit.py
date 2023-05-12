# -*- coding: utf-8 -*-
"""
Created on Fri Apr 28 09:58:21 2023

@author: jeogb
"""

import numpy as np

# For interaction parameters: any units will work, but keep track
# They determine the units of the output

#============================================================================#

N = 6 # Number of wall layers in each direction

Z = 1000 # Resolution of calculation (how many points across the pore)

eps_ff = 0.2941 # fluid-fluid interaction potential (energy units)
sig_ff = 3.73 # fluid-fluid interaction potential (distance units)

eps_ss = 0.0556 # solid-solid interaction potential (energy units)
sig_ss = 3.40 # solid-solid interaction potential (distance units)

x = (40.878 + 1.84204) # x-length of simulation box (same units as sig)
y = (48.611 + 0.709) # y-length of simulation box (same units as sig)
A = 793 # number of atoms per individual wall (approximate if not constant)

W = 30 # Pore width (same units as sig)
d = W/2 # Half pore width
pad = 0.8 # Padding from wall to prevent potential from approaching (same units as sig)

#============================================================================#

file = open('C:/Users/jeogb/.spyder-py3/Integrated_Potential_Fit.py', 'w')

Date_Time_Auth_Lines = [
    '# -*- coding: utf-8 -*-\n'
    '"""\n'
    'Created on Fri Apr 28 10:39:00 2023\n'
    '\n'
    '@author: Jason Ogbebor\n'
    '"""\n'
    '\n'
    ]
file.writelines(Date_Time_Auth_Lines)

Description = [
    '"""\n'
    '\n'
    'Fit sigma and epsilon so that a single wall can mimic six layers of graphene\n'
    'interacting with\n'
    '\n'
    '"""\n'
    '\n'
    ]
file.writelines(Description)

Preamble_Lines = [
    'import numpy as np\n'
    'import matplotlib.pyplot as plt\n'
    'import matplotlib\n'
    'import scipy.constants as c\n'
    'from scipy.optimize import minimize\n'
    '\n'
    'if True: # Plotting style guide\n'
    '\n'
    '    style = {\n'
    '        "figure.figsize": (10, 8),\n'
    '        "font.size": 36,\n'
    '        "axes.labelsize": "26",\n'
    '        "axes.titlesize": "26",\n'
    '        "xtick.labelsize": "26",\n'
    '        "ytick.labelsize": "26",\n'
    '        "legend.fontsize": "18",\n'
    '        "xtick.direction": "in",\n'
    '        "ytick.direction": "in",\n'
    '        "lines.linewidth": 2,\n'
    '        "lines.markersize": 12,\n'
    '        "xtick.major.size": 18,\n'
    '        "ytick.major.size": 18,\n'
    '        "xtick.top": True,\n'
    '        "ytick.right": True,\n'
    '    }\n'
    '    \n'
    '    matplotlib.rcParams.update(style)\n'
    '    markers = ["o", "X",  "^", "P", "d", "*", "s", ".", "x", ">"] * 5\n'
    '    matplotlib.rcParams["font.family"] = ["serif"]\n'
    '\n']
file.writelines(Preamble_Lines)

Division = [
    '#============================================================================#\n'
    '\n']

##############################################################################
file.writelines(Division)
##############################################################################

Initial_parameters = [
    f'N = {N}\n'
    '\n'
    f'W = {W} # A\n'
    'd = W/2 # A\n'
    f'pad = {pad}\n'
    '\n'
    f'Points = {Z}\n'
    f'Z = np.linspace((-d+pad), (d-pad), Points) # vertical positions relative to the center of the pore, which is z=0\n'
    '\n'
    f'x = {x}\n'
    f'y = {y}\n'
    f'n = {A} / (x*y)\n'
    '\n'
    f'eps_ss = {eps_ss} # kJ/mol\n'
    f'eps_ff = {eps_ff} # kJ/mol\n'
    'eps_fs = np.sqrt(eps_ss*eps_ff)\n'
    f'sigma_ff = {sig_ff} # A\n'
    f'sigma_ss = {sig_ss} # A\n'
    'r_0 = (sigma_ff + sigma_ss) / 2 # sigma_fs\n'
    '\n'
    ]
file.writelines(Initial_parameters)

##############################################################################
file.writelines(Division)
##############################################################################

Wall_Potential = [
    'U_walls = []\n'
    '\n'
    f'U = [0] * Points\n'
    '\n'
    f'for nw in np.linspace(0,N-1,N):\n'
    '    \n'
    '    T1 = [(((r_0 / (d + (nw*3.35) + z))**10)) for z in Z]\n'
    '    T2 = [(((r_0 / (d + (nw*3.35) - z))**10)) for z in Z]\n'
    '    T3 = [(((r_0 / (d + (nw*3.35) + z))**4)) for z in Z]\n'
    '    T4 = [(((r_0 / (d + (nw*3.35) - z))**4)) for z in Z]\n'
    '\n'
    '    U_wall = [(4*c.pi*n*eps_fs*(r_0**2))*((1/5)*(t1+t2)-(1/2)*(t3+t4)) for t1,t2,t3,t4 in zip(T1,T2,T3,T4)]\n'
    '\n'
    '    U_walls.append(U_wall)\n'
    '\n'
    '    U = [u+u_w for u,u_w in zip(U,U_wall)]\n'
    '\n'
    ]
file.writelines(Wall_Potential)

##############################################################################
file.writelines(Division)
##############################################################################

Plots = [
    'fig = plt.figure()\n'
    'plt.plot(Z, U, color="blue")\n'
    'plt.ylabel("Integrated Potential")\n'
    'plt.xlabel("z-distance")\n'
    'plt.xlim(-15,15)\n'
    'plt.ylim(-3.5,5)\n'
    '#plt.savefig("file_name.pdf", bbox_inches="tight")\n'
    'plt.show(fig)\n'
    'plt.close(fig)\n'
    '\n'
    ]
file.writelines(Plots)

##############################################################################
file.writelines(Division)
##############################################################################

Fit = [
    'R_sq = [] # keep track of the fitting process\n'
    '\n'
    'Z_fit = np.linspace(pad, (d+pad), int(Points/2)) # half of the pore size (only looking at one wall)\n'
    '\n'
    'lj_type = "lj93" # type of wall potential to fit\n'
    '\n'
    'def func(x):\n'
    '\n'
    '    eps, sig = x\n'
    '\n'
    '    if lj_type == "lj93":\n'
    '        U_fit = [eps*((2/15)*(sig/z)**9 - (sig/z)**3) for z in Z_fit]\n'
    '\n'
    '    elif lj_type == "lj1043":\n'
    '        U_fit = [2*c.pi*eps*((2/5)*(sig/z)**10 - (sig/z)**4 - ((np.sqrt(2)*sig**3) / (3*(z + (0.61/np.sqrt(2))*sig)**3))) for z in Z_fit]\n'
    '\n'
    '    elif lj_type == "lj126":\n'
    '        U_fit = [4*eps*((sig/z)**12 - (sig/z)**6) for z in Z_fit]\n'
    '\n'
    '    else:\n'
    '        raise Exception("Unrecognized lj_type")\n'
    '\n'
    '    r_sq = np.corrcoef(U_fit, U[0:500])[0, 1]\n'
    '    R_sq.append(r_sq)\n'
    '\n'
    '    return (1-r_sq)\n'
    '\n'
    'x0 = [2.535, 4.4] # Initial guesses (play around with these values to see how it affects the fit)\n'
    '\n'
    'params = minimize(func, x0, tol=1e-6, bounds=((0,3),(2.5,5)))\n'
    '\n'
    'eps = params.x[0]\n'
    'sig = params.x[1]\n'
    '\n'
    'if lj_type == "lj93":\n'
    '    U_fit = [eps*((2/15)*(sig/z)**9 - (sig/z)**3) for z in Z_fit]\n'
    '\n'
    'if lj_type == "lj1043":\n'
    '    U_fit = [2*c.pi*eps*((2/5)*(sig/z)**10 - (sig/z)**4 - ((np.sqrt(2)*sig**3) / (3*(z + (0.61/np.sqrt(2))*sig)**3))) for z in Z_fit]\n'
    '\n'
    'if lj_type == "lj126":\n'
    '    U_fit = [4*eps*((sig/z)**12 - (sig/z)**6) for z in Z_fit]\n'
    ]
file.writelines(Fit)

##############################################################################
file.writelines(Division)
##############################################################################

Plots = [
    'fig = plt.figure()\n'
    'plt.plot(Z[0:500], U[0:500], color="blue", label="Integrated Potential")\n'
    'plt.plot(Z[0:500], U_fit, color="red", linestyle="--", label="Fit")\n'
    'plt.ylabel("Potential Energy")\n'
    'plt.xlabel("z-distance")\n'
    'plt.xlim(-15,0)\n'
    'plt.ylim(-3.5,5)\n'
    'plt.legend(loc="best")\n'
    '#plt.savefig("file_name.pdf", bbox_inches="tight")\n'
    'plt.show(fig)\n'
    'plt.close(fig)\n'
    '\n'
    ]
file.writelines(Plots)

file.close()

exec(open("Integrated_Potential_Fit.py").read())


