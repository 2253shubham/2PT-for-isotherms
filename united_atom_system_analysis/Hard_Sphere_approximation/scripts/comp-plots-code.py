# %run comp-plots-code.py averaged_out/average_data/vac-CH4-308K-10-average-prod1.xvg averaged_out/ref_prop/ref_dos_avg.txt averaged_out/ref_prop/ref_prop_avg.txt averaged_out/average_data/time.txt  

# Below is the code for computing thermodynamic parameters of a 2PT system

import sys
import numpy as np 
from scipy import integrate
from scipy import signal
import scipy.constants
import cmath
import math
import matplotlib.pyplot as plt
from scipy.fftpack import fft
from scipy import optimize

Na = scipy.constants.Avogadro # Avogadro's number
kb = scipy.constants.k # Boltzmann Constant k
data1 = sys.argv[1]
data2 = sys.argv[2]
#dataout = sys.argv[3]
#temp = sys.argv[4]

data3 = sys.argv[3]
data4 = sys.argv[4]
data5 = sys.argv[5]
data6 = sys.argv[6]
dataout = sys.argv[7]
data1a = np.loadtxt(fname = data1)
data2a = np.loadtxt(fname = data2)
data3a = np.loadtxt(fname = data3)
data4a = np.loadtxt(fname = data4)
data5a = np.loadtxt(fname = data5)
data6a = np.loadtxt(fname = data6)

dataav = (data1a + data2a + data3a + data4a + data5a + data6a)/6

#--------------------------------------------------------------------------------------------------------------------

N_org = np.linspace(18,738,21)/18
#N_org = np.insert(N_org, 0, 0)
N_org = np.append(N_org,[756/18, 792/18, 828/18, 864/18])
N_org = N_org[0:]
N = np.linspace(0,864,865)/18
N = N[0:]

datat = np.vstack((data1a, data2a, data3a, data4a, data5a, data6a))
dataav = np.mean(datat, axis=0)
datastd = np.std(datat, axis=0)
datacf = (datastd/6**0.5)*2.447

# all parameters and constants

out = open(dataout+"-avg.txt","w")
out.truncate()
for i in range(len(dataav)):
    print(N_org[i],dataav[i],datacf[i],datastd[i], file = out)
    #print("\n", file = out)
plt.figure(1)
plt.plot()

plt.figure(1)
plt.errorbar(N_org, dataav, datastd, linestyle="none", marker='*', capsize=2)
plt.ylabel(dataout)
plt.xlabel("Loading [molecules/unit cell]")
plt.show()
plt.savefig(dataout+"_vs_N_plot.png")
