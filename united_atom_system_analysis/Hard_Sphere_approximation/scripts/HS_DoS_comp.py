# Below is the code for computing thermodynamic parameters of a 2PT system

import sys
import numpy as np
from scipy import integrate
from scipy import signal
import scipy.constants
import cmath
import math
import matplotlib.pyplot as plt
from scipy.fftpack import fft, fftfreq
from scipy import optimize

infile = sys.argv[1]  # system argument to specify input file
timefile = sys.argv[2]  # time file
outfile = sys.argv[3]  # system argument to specify output file
enerfile = sys.argv[4]  # system argument to specify energy of md simulations
temp = sys.argv[5]  # system temperature
volume = sys.argv[6]  # volume original in nm^3
sf = sys.argv[7]  # volume scaling to original
minm = sys.argv[8]  # minimum system size
mass = sys.argv[9]  # mass of one particle in g/mol
rn = sys.argv[10]  # system size
nr = sys.argv[11]  # system run


# Snippet to get run number and number of particles in the system

numstr = 0  # string to store number of components in the systems tested
numrun = 0  # string to store run number for the system being tested
for j in range(len(infile)):
    if (infile[j] == "K") and (infile[j + 1] == "-"):
        u1 = j + 2
    elif (infile[j] == "a") and (infile[j + 1] == "v"):
        u2 = j - 1
        u3 = j + 3
    elif (infile[j] == "p") and (infile[j + 1] == "r"):
        u4 = j - 1
    else:
        continue
# numstr = infile[u1:u2]
numstr = rn
print(numstr)
# numrun = infile[u3:u4]
numrun = nr
print(numrun)
# Snippet to get run number and number of particles in the system

# --------------------------------------------------------------------------------------------------------------------

# reading file and performing operations

# data = np.loadtxt(fname = infile)
data = np.loadtxt(fname=infile, skiprows=17, usecols=1)
enermdall = np.loadtxt(fname=enerfile)
t = np.loadtxt(fname=timefile)
sf = float(sf)
# all parameters and constants

# t = [] # list to store time
out_list = []  # list to store auto-correlation values
kb = scipy.constants.k  # Boltzmann Constant k
Na = scipy.constants.Avogadro  # Avogadro's number
T = float(temp)  # system temperature
ff = []  # storing frequency values
rc = (
    []
)  # integral for partition function corresponding to classical harmonic oscillator
rq = []  # integral for partition function corresponding to quantum harmonic oscillator
y = []  # store DOS overall
ys = []  # store DOS solid contribution
yg = []  # store DOS gas contribution
m = float(mass) * (1e-3) / Na  # mass of one methane molecule
vol = float(volume) * (1e-27) * float(sf)  # volume of the system
beta = 1 / (kb * T)
h = scipy.constants.Planck  # Planck's constant
W_Acho = (
    []
)  # weighting function for helmholtz free energy corresponding to solid assuming classical harmonic oscillator
W_Acho.append(0)
W_Aqho = (
    []
)  # weighting function for helmholtz free energy corresponding to solid assuming quantum harmonic oscillator
W_Aqho.append(0)
W_Scho = (
    []
)  # weighting function for entropy corresponding to solid assuming classical harmonic oscillator
W_Scho.append(0)
W_Sqho = (
    []
)  # weighting function for entropy corresponding to solid assuming quantum harmonic oscillator
W_Sqho.append(0)
W_Eqho = (
    []
)  # weighting function for energy corresponding to solid assuming quantum harmonic oscillator
W_Eqho.append(0)
ysm1 = (
    []
)  # modified solid DOS parameter used to calculate differnt properties assuming classical harmonic oscillator
ysm1.append(0)
ysm2 = (
    []
)  # modified solid DOS parameter used to calculate differnt properties assuming quantum harmonic oscillator
ysm2.append(0)
ygm = []  # modified gas DOS parameter used to calculate different properties
ygm.append(0)

# --------------------------------------------------------------------------------------------------------------------

# storing time and correlation values to corresponding lists

N = int(numstr)  # number of components in the system
den = N / (vol)  # particle number density of the system
# length = len(data) # calculating length of the data lists
out_list = data[0:] * m * N * (1e6)  # storing mass weighted auto correlation values
dt = (t[1] - t[0]) * (1e-12)  # time step length

# --------------------------------------------------------------------------------------------------------------------

# Normalizing data such that C(0) = 3*N*kb*T

out_list = (out_list / out_list[0]) * 3 * N * kb * T
length = len(out_list)  # calculating length of the data lists
out_list = np.hstack(
    (out_list, out_list[::-1])
)  # making autocorrelation values symmetric
length *= 2

# --------------------------------------------------------------------------------------------------------------------
"""
# performing discrete fourier cosine transform on auto correlation values (1)(takes a lot of time to run and is similar to dct in python)

angles = 2*np.pi*np.arange(length, dtype=np.float_)/length
y = np.empty(length, dtype = np.float_) # store DOS overal
for i in range(length):
    y[i] = np.sum(out_list*np.cos(i*angles)) 
#print(y)

#--------------------------------------------------------------------------------------------------------------------
"""

# performing discrete fourier transform on auto correlation values using inbult python functions(2)

yn = np.empty(length, dtype=np.float_)  # store DOS overal
yn = fft(out_list)
yn = np.real(yn)
y = signal.savgol_filter(yn, 21, 5)

# --------------------------------------------------------------------------------------------------------------------

# calculating total degrees of freedom and self diffusion coefficient of the system

yn = yn * (beta) * 2 * dt  # multiplying by the factor (2/kT)*dt
y = y * (beta) * 2 * dt  # multiplying by the factor (2/kT)*dt
ff = fftfreq(len(yn), d=dt)  # v (frequency)
dif_sp = y[0] * (T * kb) / (12 * m * N)  # Self diffusion coefficient (p)
y = y[: len(y) // 2].reshape(-1, 1)
ff = ff[: len(ff) // 2].reshape(-1, 1)
# print(tdos)

# --------------------------------------------------------------------------------------------------------------------

# Calculate solid and gas fraction of dos (here y) along with several diffusion coefficients
# Solid and gas dos vs frequency also calculated here

delta = (
    math.sqrt((np.pi * kb * T) / m)
    * ((6 / np.pi) ** (2 / 3))
    * (den ** (1 / 3))
    * (2 * y[0])
    / (9 * N)
)  # normalized diffusivity constant

# Snippet used to calculate packing fraction and fraction of gas

"""
def packing_fraction(pf):
    return (2*(delta**3)*(pf**5)-(delta**2)*(pf**(7/3))-6*(pf**(10/3))*(delta**2)+2*delta*(pf**(2/3))+6*delta*(pf**(5/3))-2)
pf_value = optimize.root_scalar(packing_fraction,x0 = 0,x1 = 1)
pf_assval = pf_value.root
"""

# Snippet used to calculate packing fraction and fraction of gas

q = 0


def fluidicity(f):
    return (
        2 * (delta ** (0 - 9 / 2)) * (f ** (4.5 * q + 15 / 2))
        - 6 * (delta ** (0 - 3)) * (f ** (3 * q + 5))
        - (delta ** (0 - 3 / 2)) * (f ** (2.5 * q + 7 / 2))
        + 6 * (delta ** (0 - 3 / 2)) * (f ** (1.5 * q + 5 / 2))
        + 2 * (f ** (q + 1))
        - 2
    )


f_value = optimize.root_scalar(fluidicity, bracket=[0, 1], method="brentq")
frac = f_value.root

# frac = delta*(pf_assval**(2/3)) # fraction of gas
pf_assval = (frac ** (3 / 2)) * (delta ** (-3 / 2))  # packing fraction
dia = (pf_assval * 6 / (np.pi * den)) ** (1 / 3)  # hard sphere diameter
dif_hsfp = (kb * T * y[0]) / (
    12 * m * N * frac
)  # diffusivity of the gas component assuming gas as a hard sphere (fp)
dif_cep = (
    (3 / 8) * math.sqrt((kb * T) / (np.pi * m)) / (den * dia**2)
)  # hard sphere Chapman Enskog diffusion (p)
dif_cefp = dif_cep / frac  # hard sphere Chapman Enskog diffusion (fp)

"""
# some verification

c1 = frac*pf_assval
z = (1 + c1 + c1**2 - c1**3)/(1 - c1)**3 # Carnahan Starling EOS
coeff = (4*c1)/(z-1)
k = coeff * dif_cefp/dif_hsfp
print(k)
"""

# Calculating separate total solid and gas dos

p = (np.pi * y[0] * ff) / (2 * frac * 3 * N)
yg = y[0] / (1 + (p**2))  # calculating gas DOS function vs frequency
ys = y - yg  # calculating solid DOS function vs frequency


out1 = open("total-vdos-N-" + numstr + ".txt", "w")
out2 = open("solid-vdos-N-" + numstr + ".txt", "w")
out3 = open("gas-vdos-N-" + numstr + ".txt", "w")
out4 = open("frequency.txt", "w")
out1.truncate()
out2.truncate()
out3.truncate()
out4.truncate()
for i in range(len(y)):
    print(float(y[i]) * 3e10, file=out1)
    print(float(ys[i]) * 3e10, file=out2)
    print(float(yg[i]) * 3e10, file=out3)
    print(float(ff[i]) / 3e10, file=out4)
out1.close()
out2.close()
out3.close()
out4.close()


# 2017 paper implementation
"""
flag = 1
for i in range((int)((len(y))/2))[::-1]:
    if ((y[i]-yg[i])>0 and (y[i-1]-yg[i-1]>0)):
        ka = i
        break
    if (y[i]>0 and y[i-1]>0 and flag == 1):
        kc = i    
        flag = 0 
#print(ff[ka]/3e10)
#print(ff[kc]/3e10)
"""
tdosg = integrate.simps(yg.T, ff[: len(yg)].T)  # total number of gas DOS
tdoss = integrate.simps(ys.T, ff[: len(ys)].T)  # total number of solid DOS
tdos = tdosg + tdoss

# finding partition function for the system and different thermodynamical properties

# calculating reference energy

enermd = enermdall[int((int(N / int(minm)) - 1) / 2)]  # converting to appropriate units
print(enermd)
refener = enermd - ((3 * N) * (1 - (0.5 * frac)) / (beta))  # reference energy

# finding entropy contribution for hard sphere fluid

# Step 1 : find excess entropy contribution per molecule assuming monatomic ideal gas at same temperature and density (united atom model), per molecule basis

r = (((2 * np.pi * m * kb * T) / (h**2)) ** (1.5)) * (np.exp(2.5)) / (den * frac)
si = kb * np.log(r)  # entropy value per molecule basis

# Step 2 : use Carnahan Starling equation to calculate the hard sphere excess entropy, per molecule basis

c1 = frac * pf_assval
z = (1 + c1 + c1**2 - c1**3) / (1 - c1) ** 3  # Carnahan Starling EOS
r1 = np.log(z)
r2 = (c1 * (3 * c1 - 4)) / ((1 - c1) ** 2)
r3 = r1 + r2
# sh = si + kb*r3 # excess entropy value of hard sphere fluid per molecule basis
sh = si + kb * r2  # excess entropy value of hard sphere fluid per molecule basis

# Calculating Entropy of the system

W_Shs = sh / (3 * kb)  # weighting function for entropy corresponding to hard sphere gas
W_Scho = (1 - np.log(beta * h * ff[1 : len(ff)])).reshape(-1, 1)
ysm1 = W_Scho * ys[1 : len(ff)]
a1 = ((beta * h * ff[1 : len(ff)]) / (np.exp(beta * h * ff[1 : len(ff)]) - 1)).reshape(
    -1, 1
)
a2 = (np.log(1 - np.exp(0 - beta * h * ff[1 : len(ff)]))).reshape(-1, 1)
W_Sqho = a1 - a2
ysm2 = W_Sqho * ys[1 : len(ff)]
S_cho1 = (
    integrate.simps(ysm1.T, ff[1 : len(ff)].T)
) * kb  # Entropy (assuming classical harmonic oscillator)
S_qho1 = (
    integrate.simps(ysm2.T, ff[1 : len(ff)].T)
) * kb  # Entropy (assuming quantum harmonic oscillator)
S_hs = (
    integrate.simps(yg[1 : len(yg)].T * W_Shs, ff[1 : len(yg)].T)
) * kb  # Entropy (assuming quantum harmonic oscillator)
S_cho = S_cho1 + S_hs
S_qho = S_qho1 + S_hs

W_Ahs = 0.5 - (
    sh / (3 * kb)
)  # weighting function for helmholtz free energy corresponding to hard sphere gas
W_Acho = (np.log(beta * h * ff[1 : len(ff)])).reshape(-1, 1)
ysm1 = W_Acho * ys[1 : len(ff)]
q = (1 - np.exp(0 - beta * h * ff[1 : len(ff)])) / np.exp(
    0 - beta * h * ff[1 : len(ff)]
).reshape(-1, 1)
W_Aqho = np.log(q)
ysm2 = W_Aqho * ys[1 : len(ff)]

A_cho1 = refener + (integrate.simps(ysm1.T, ff[1 : len(ff)].T)) * beta ** (
    -1
)  # Helmholtz free energy (assuming classical harmonic oscillator)
A_qho1 = refener + (integrate.simps(ysm2.T, ff[1 : len(ff)].T)) * beta ** (
    -1
)  # Helmholtz free energy (assuming quantum harmonic oscillator)
A_hs = (integrate.simps(yg[1 : len(yg)].T * W_Ahs, ff[1 : len(yg)].T)) * beta ** (-1)
A_cho = A_cho1 + A_hs
A_qho = A_qho1 + A_hs
Q = 0 - (A_cho) / (
    kb * T
)  # partition function of the system (assuming classical harmonic oscillator)
Q_qho = 0 - (A_qho) / (
    kb * T
)  # partition function of the system (assuming quantum harmonic oscillator)

# --------------------------------------------------------------------------------------------------------------------

# Printing results in a file
out = open(outfile, "w")
out.truncate(0)
print("Self diffusion coefficient = ", np.float(dif_sp), file=out)
print("Fluidicity = ", np.float(frac), file=out)
print("Delta = ", np.float(delta), file=out)
print("Density = ", np.float(den), file=out)
print("HS Diameter (in nm) = ", np.float(dia * 1e9), file=out)
print("Packing fraction of gas component = ", np.float(pf_assval * frac), file=out)
print("Total Density of states in the system = ", np.float(tdos), file=out)
print("Total number of gas DOS in the system = ", np.float(tdosg), file=out)
print("Total number of solid DOS in the system = ", np.float(tdoss), file=out)
print(
    "Entropy of the system (assuming classical harmonic oscillator for solid part of the system) per molecule basis = ",
    np.float(S_cho / (N * kb)),
    file=out,
)
print(
    "Entropy of the system (assuming quantum harmonic oscillator for solid part of the system) per molecule basis = ",
    np.float(S_qho / (N * kb)),
    file=out,
)
print(
    "Solid Entropy contribution (assuming classical harmonic oscillator) per molecule basis = ",
    np.float(S_cho1 / (N * kb)),
    file=out,
)
print(
    "Solid Entropy contribution (assuming quantum harmonic oscillator) per molecule basis = ",
    np.float(S_qho1 / (N * kb)),
    file=out,
)
print(
    "Gas Entropy contribution per molecule basis = ",
    np.float(S_hs / (N * kb)),
    file=out,
)
print(
    "Entropy of the system (assuming classical harmonic oscillator for solid part of the system) scales with N = ",
    np.float(S_cho / (sf * kb)),
    file=out,
)
print(
    "Entropy of the system (assuming quantum harmonic oscillator for solid part of the system) scales with N = ",
    np.float(S_qho / (sf * kb)),
    file=out,
)
print(
    "Solid Entropy contribution (assuming classical harmonic oscillator) scales with N = ",
    np.float(S_cho1 / (sf * kb)),
    file=out,
)
print(
    "Solid Entropy contribution (assuming quantum harmonic oscillator) scales with N = ",
    np.float(S_qho1 / (sf * kb)),
    file=out,
)
print("Entropy hs wf = ", np.float(W_Shs), file=out)
print("Gas Entropy contribution scales with N = ", np.float(S_hs / (sf * kb)), file=out)
print(
    "Helmholtz free energy of the system (assuming classical harmonic oscillator for solid part of the system) per molecule basis = ",
    np.float(A_cho / (N * T * kb)),
    file=out,
)
print(
    "Helmholtz free energy of the system (assuming quantum harmonic oscillator for solid part of the system) per molecule basis = ",
    np.float(A_qho / (N * T * kb)),
    file=out,
)
print(
    "Solid Helmholtz free energy contribution (assuming classical harmonic oscillator) per molecule basis = ",
    np.float(A_cho1 / (N * T * kb)),
    file=out,
)
print(
    "Solid Helmholtz free energy contribution (assuming quantum harmonic oscillator) per molecule basis = ",
    np.float(A_qho1 / (N * T * kb)),
    file=out,
)
print(
    "Gas Helmholtz free energy contribution per molecule basis = ",
    np.float(A_hs / (N * T * kb)),
    file=out,
)
print(
    "Helmholtz free energy of the system (assuming classical harmonic oscillator for solid part of the system) scales with N = ",
    np.float(A_cho / (sf * T * kb)),
    file=out,
)
print(
    "Helmholtz free energy of the system (assuming quantum harmonic oscillator for solid part of the system) scales with N = ",
    np.float(A_qho / (sf * T * kb)),
    file=out,
)
print(
    "Solid Helmholtz free energy contribution (assuming classical harmonic oscillator) scales with N = ",
    np.float(A_cho1 / (sf * T * kb)),
    file=out,
)
print(
    "Solid Helmholtz free energy contribution (assuming quantum harmonic oscillator) scales with N = ",
    np.float(A_qho1 / (sf * T * kb)),
    file=out,
)
print("Helmholtz hs wf = ", np.float(W_Ahs), file=out)
print(
    "Gas Helmholtz free energy contribution scales with N = ",
    np.float(A_hs / (sf * T * kb)),
    file=out,
)
print(
    "Total Energy of the system (assuming classical harmonic oscillator for solid part of the system) per molecule basis = ",
    np.float(A_cho / (N * T * kb) + S_cho / (N * kb)),
    file=out,
)
print(
    "Total Energy of the system (assuming quantum harmonic oscillator for solid part of the system) per molecule basis = ",
    np.float(A_qho / (N * T * kb) + S_qho / (N * kb)),
    file=out,
)
print(
    "Solid contribution to Total Energy (assuming classical harmonic oscillator) per molecule basis = ",
    np.float(A_cho1 / (N * T * kb) + S_cho1 / (N * kb)),
    file=out,
)
print(
    "Solid contribution to Total Energy (assuming quantum harmonic oscillator) per molecule basis = ",
    np.float(A_qho1 / (N * T * kb) + S_qho1 / (N * kb)),
    file=out,
)
print(
    "Gas contribution to Total Energy per molecule basis = ",
    np.float(A_hs / (N * T * kb) + S_hs / (N * kb)),
    file=out,
)
print(
    "Total Energy of the system (assuming classical harmonic oscillator for solid part of the system) scales with N = ",
    np.float(A_cho / (sf * T * kb) + S_cho / (sf * kb)),
    file=out,
)
print(
    "Total Energy of the system (assuming quantum harmonic oscillator for solid part of the system) scales with N = ",
    np.float(A_qho / (sf * T * kb) + S_qho / (sf * kb)),
    file=out,
)
print(
    "Solid contribution to Total Energy (assuming classical harmonic oscillator) scales with N = ",
    np.float(A_cho1 / (sf * T * kb) + S_cho1 / (sf * kb)),
    file=out,
)
print(
    "Solid contribution to Total Energy (assuming quantum harmonic oscillator) scales with N = ",
    np.float(A_qho1 / (sf * T * kb) + S_qho1 / (sf * kb)),
    file=out,
)
print(
    "Gas contribution to Total Energy scales with N = ",
    np.float(A_hs / (sf * T * kb) + S_hs / (sf * kb)),
    file=out,
)
print("ln of partition function cho hs =", np.float(Q), file=out)
print("ln of partition function qho hs =", np.float(Q_qho), file=out)

out.close()

# --------------------------------------------------------------------------------------------------------------------

print("done!")  # indication that the particular run is completed
