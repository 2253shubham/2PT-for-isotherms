import numpy as np
import matplotlib.pyplot as plt
import scipy.constants
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.optimize import curve_fit
import sys

# exp_data = sys.argv[1]
GCN_cho = []
GCN_qho = []

kb = scipy.constants.k  # Boltzmann Constant k
Na = scipy.constants.Avogadro  # Avogadro's number
T = 300  # system temperature
beta = 1 / (kb * T)
h = scipy.constants.Planck  # Planck's constant
m = 18.06 * (1e-3) / Na  # mass of one water molecule
N_org = np.arange(25, 375, 25)
N_org = N_org[0:]
Npre_cho = []
chepot_cho = []
pressure_cho = []
Npre_qho = []
chepot_qho = []
pressure_qho = []
u1 = []
u2 = []
p1 = []
p2 = []
p = []
a1 = []
a2 = []
min1 = []
min2 = []
slope = []
# slope.append(0)
sl1 = []
sl1.append(0)
N = np.linspace(0, 350, 351)
N = N[0:]
P_cho = []
P_qho = []
Gc = []
Gq = []
Ndc = []
Ndq = []
q1 = []

# chemical potential
r = (((2 * np.pi * m * kb * T) / (h**2)) ** (1.5)) * kb * T
u0 = (
    0
    - kb * T * np.log(r)
    - kb
    * T
    * np.log(
        (((np.pi) ** (0.5)) * 0.5)
        * ((T**3) / (13.706211107457026 * 20.998505222742857 * 39.467682045442515))
        ** (0.5)
    )
)
u11 = u0 + kb * T * np.log(1e5)
u21 = u0 + kb * T * np.log(3e9)
u = np.linspace(u11, u21, 250)
out1 = open("chem_pot._values_1e1_3e4_Pa.txt", "w")
out1.truncate(0)
for i in range(len(u)):
    print(u[i], file=out1)
out1.close()
# ------------------------------------------------------------

datao1 = np.loadtxt("lnQ.txt")
# datao2 = np.loadtxt("Average_ln_Partition_Function_qho.txt")
data3 = np.loadtxt("chem_pot._values_1e1_3e4_Pa.txt")
length3 = len(data3)

spl1 = InterpolatedUnivariateSpline(N_org, datao1)
# spl2 = InterpolatedUnivariateSpline(N_org,datao2)
data1 = spl1(N)
# data2 = spl2(N)
length1 = len(data1)

out1 = open("ln-part-fn-cho-interpolated.txt", "w")
out1.truncate(0)
for i in range(len(data1)):
    print(data1[i], file=out1)
out1.close()


def func(x, a, b, c, d, e):
    return a + b * x + c * (x**2) + d * (x**3) + e * (x**4)


popt, pconv = curve_fit(func, N, data1)
a = popt[0]
b = popt[1]
c = popt[2]
d = popt[3]
e = popt[4]

plt.figure(3)
plt.plot(N_org, datao1, "b--", label="original data")
plt.plot(N, data1, "r", label="interpolated data")
plt.plot(N, func(N, *popt), "k--", label="curve_fit_4th_order")
plt.xlabel("N")
plt.ylabel("lnQ")
plt.legend()

for i in data3:
    r1 = 0
    r2 = 0
    print()
    for j in range(length1):
        r1 += np.exp(data1[j] + beta * i * N[j])
        # r2 += np.exp(data2[j] + beta*i*N[j])
    Gc.append(r1)
    # Gq.append(r2)
out1 = open("GCN-part-fn-cho.txt", "w")
out2 = open("ln-GCN-part-fn-cho.txt", "w")
out1.truncate(0)
out2.truncate(0)
for i in range(length3):
    print(Gc[i], file=out1)
    print(np.log(Gc[i]), file=out2)
out1.close()
out2.close()

for i in range(length1):
    slope.append(b + 2 * c * N[i] + 3 * d * (N[i] ** 2) + 4 * e * (N[i] ** 3))

out1 = open("Slope-lnQ-vs-N-cho_new.txt", "w")
# out1 = open("Slope-lnQ-vs-N-cho.txt","w")
# out1 = open("Slope_lnQ_vs_N_qho.txt","w")
out1.truncate(0)
for i in range(length1):
    print(slope[i], file=out1)
out1.close()

plt.figure(2)
plt.plot(N, slope)
plt.xlabel("N")
plt.ylabel("lnQvsN_slope")

data = np.loadtxt("ln-GCN-part-fn-cho.txt")
length = len(data)
slope = []
data3mod = data3[0:] * 1e20
"""
def func(x, a, b, c, d, e, f):
	return a + b*x + c*(x**2) + d*(x**3) + e*(x**4) + f*(x**5) 

popt , pconv = curve_fit(func,data3mod,data)
a = popt[0]
b = popt[1]
c = popt[2]
d = popt[3]
e = popt[4]
f = popt[5]


for i in range(len(data)):
	slope.append((b + 2*c*data3mod[i] + 3*d*(data3mod[i]**2) + 4*e*(data3mod[i]**3) + 5*f*(data3mod[i]**4))*1e20)
	# + 6*g*(data3mod[i]**5) g = popt[6] + g*(x**6) , g

plt.figure(4)
plt.plot(data3,data,"r--",label = "data")
plt.plot(data3,func(data3mod,*popt),"k--",label = "curve_fit_5th_order")
plt.xlabel("chem.pot.")
plt.ylabel("lnGCN")
plt.legend()

out1 = open("Slope-lnGCNpf-vs-chempot-cho_newp.txt","w")
out1.truncate(0)
for i in range(length):
	print(slope[i],file = out1)
out1.close()

data4 = np.loadtxt("Slope-lnGCNpf-vs-chempot-cho_newp.txt")
#data5 = np.loadtxt("Slope_lnGCNpf_vs_chempot_qho.txt")
length4 = len(data4)
for i in range(length4):
	a = data4[i]/beta
	Npre_cho.append(a)
	chepot_cho.append(data3[i])
	pressure_cho.append(np.exp((data3[i]-u0)/(kb*T)))
"""
plt.figure(1)
# plt.plot(pressure_cho,Npre_cho,"b--",label = "cho_GCN_5th-order_poly_fitted")
ax = plt.gca()
ax.set_xscale("log")
# plt.plot(pressure_cho,Npre_qho,label = "qho_GCN")
plt.xlabel("Pressure in SI units")
plt.ylabel("Loading")

datac1 = np.loadtxt(fname="Slope-lnQ-vs-N-cho_new.txt")
# data2 = np.loadtxt(fname = "Slope_lnQ_vs_N_qho.txt")
lengthc1 = len(datac1)
for i in range(lengthc1):
    u1.append(-kb * T * datac1[i])
    # u2.append(-kb*T*data2[i])
    p1.append(np.exp((u1[i] - u0) / (kb * T)))
    # p2.append(np.exp((u2[i]-u0)/(kb*T)))
plt.figure(1)
plt.plot(p1, N, label="cho_CN_4th-order_poly_fitted")
ax.set_xscale("log")
plt.show()
