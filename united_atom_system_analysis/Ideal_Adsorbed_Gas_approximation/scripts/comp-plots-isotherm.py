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

Na = scipy.constants.Avogadro  # Avogadro's number
kb = scipy.constants.k  # Boltzmann Constant k
data1 = sys.argv[1]
data2 = sys.argv[2]
# dataout = sys.argv[3]
# temp = sys.argv[4]

data3 = sys.argv[3]
data4 = sys.argv[4]
data5 = sys.argv[5]
data6 = sys.argv[6]
data7 = sys.argv[7]
data8 = sys.argv[8]
data9 = sys.argv[9]
data10 = sys.argv[10]
dataout = sys.argv[11]
data1a = np.loadtxt(fname=data1, usecols=1)  # 0 for CN
data2a = np.loadtxt(fname=data2, usecols=1)  # 0 for CN
data3a = np.loadtxt(fname=data3, usecols=1)  # 0 for CN
data4a = np.loadtxt(fname=data4, usecols=1)  # 0 for CN
data5a = np.loadtxt(fname=data5, usecols=1)  # 0 for CN
data6a = np.loadtxt(fname=data6, usecols=1)  # 0 for CN
data7a = np.loadtxt(fname=data7, usecols=1)  # 0 for CN
data8a = np.loadtxt(fname=data8, usecols=1)  # 0 for CN
data9a = np.loadtxt(fname=data9, usecols=1)  # 0 for CN
data10a = np.loadtxt(fname=data10, usecols=1)  # 0 for CN
datag = np.loadtxt(fname=data1, usecols=0)  # 1 for CN

# dataav = (data1a + data2a + data3a + data4a + data5a + data6a + data7a + data8a + data9a + data10a)/10

# --------------------------------------------------------------------------------------------------------------------

N_org = np.arange(20, 361, 20) / 24
N_org = N_org[0:]
N = np.linspace(0, 360, 361)
N = N[0:]

datat = np.vstack(
    (data1a, data2a, data3a, data4a, data5a, data6a, data7a, data8a, data9a, data10a)
)
dataav = np.mean(datat, axis=0)
datastd = np.std(datat, axis=0)
datacf = (datastd / 10**0.5) * 2.228  # change this
# datacf = (datastd/4**0.5)*2.776

# all parameters and constants

out = open(dataout + "-avg.txt", "w")
out.truncate()
for i in range(0, len(dataav), 65):  # 23 for CN
    print(datag[i], dataav[i], datacf[i], datastd[i], file=out)
    # print("\n", file = out)
plt.figure(1)
plt.plot()

plt.figure(1)
plt.errorbar(datag, dataav, yerr=datastd, linestyle="none", marker="*", capsize=2)
# plt.errorbar(dataav, data7, xerr=datastd, linestyle="none", marker='*', capsize=2) # for CN
plt.xscale("log")
plt.ylabel(dataout)
plt.xlabel("Fugacity [kPa]")
plt.show()
plt.savefig(dataout + "_vs_fugacity_plot.png")
