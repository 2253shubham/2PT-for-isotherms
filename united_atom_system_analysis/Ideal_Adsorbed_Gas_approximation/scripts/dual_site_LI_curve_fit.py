import numpy as np
import matplotlib.pyplot as plt
import scipy.constants
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.optimize import curve_fit
import sys
from sympy import *

exp_data = sys.argv[1]
kb = scipy.constants.k  # Boltzmann Constant k
Na = scipy.constants.Avogadro  # Avogadro's number
T = 308  # system temperature
beta = 1 / (kb * T)
beta = beta * 1e-20
h = scipy.constants.Planck  # Planck's constant
m = 39.948 * (1e-3) / Na  # mass of one methane molecule
N_org = np.linspace(8, 96, 12) / 24
N_org = np.insert(N_org, 0, 0)
N_org = N_org[0:]
N = np.linspace(0, 96, 97) / 24
N = N[0:]
N[0] = 1e-7
Npre_cho = []
chepot_cho = []
pressure_cho = []
r = (((2 * np.pi * m * kb * T) / (h**2)) ** (1.5)) * kb * T
u0 = 0 - kb * T * np.log(r)
u1 = []
p1 = []

r = (((2 * np.pi * m * kb * T) / (h**2)) ** (1.5)) * kb * T
u0 = 0 - kb * T * np.log(r)
u11 = u0 + kb * T * np.log(5e5)
u21 = u0 + kb * T * np.log(5e9)
u = np.linspace(u11, u21, 1000)
out1 = open("chem_pot._values_Pa.txt", "w")
out2 = open("considered_pressure_values.txt", "w")
out1.truncate(0)
out2.truncate(0)
for i in range(len(u)):
    print(u[i], file=out1)
    print(np.exp((u[i] - u0) / (kb * T)), file=out2)
out1.close()
out2.close()

data1 = np.loadtxt("ln-part-fn-cho.txt") / 24
datalnQ = data1
data2 = np.loadtxt("chem_pot._values_Pa.txt") * 1e20
# data3 = np.loadtxt("N_predicted_from_GCN.txt")

"""
def func(x, b, c):
	return (b*np.exp(beta*(x+c)))/(1+np.exp(beta*(x+c))) 

popt, pconv = curve_fit(func, data2, data3[:,1])

b1 = popt[0]
c1 = popt[1]

def func(x, a, b, c, d):
	return (a*np.exp(beta*(x+c)))/(1+np.exp(beta*(x+c))) + (b*np.exp(beta*(x+d)))/(1+np.exp(beta*(x+d)))

popt, pconv = curve_fit(func, data2, data3[:,1], p0=[1, 1, 1, 1])
a1 = popt[0]
b1 = popt[1]
c1 = popt[2]
d1 = popt[3]


plt.figure(1)
plt.plot(data2,data3[:,1],"b--",label = "cho_GCN_5th-order_poly_fitted")
plt.plot(data2,func(data2,*popt),"k--",label = "curve_fit_single_site_LI_model")
#ax = plt.gca()
#ax.set_xscale("log")
#plt.plot(pressure_cho,Npre_qho,label = "qho_GCN")
plt.xlabel("Chemical potential (* 1e-20)")
plt.ylabel("Loading")
plt.show()
"""


def func(x, a, b):
    return a * np.log(a) - x * np.log(x) - (a - x) * np.log(a - x) + x * np.log(b)


popt1, pconv1 = curve_fit(
    func, N_org[1:], data1[1:], p0=[800, 10]
)  # p0=[100,100,-8,5,-100]
a1 = popt1[0]
b1 = popt1[1]

out1 = open("single-site-LI-parameters.txt", "w")
out1.truncate(0)
print(-a1, file=out1)
print(-b1, file=out1)
out1.close()

N_org[0] = 1e-7

q2 = a1 * np.log(a1) - N * np.log(N) - (a1 - N) * np.log(a1 - N) + N * np.log(b1)

"""
def func(x, a, b, c, d, e):
    t1 = (x-a)*np.exp(beta*c)+(x-b)*np.exp(beta*d)
    t2 = (a+b-x)*np.exp(beta*(c+d))
    K = (t1+((t1**2)+(4*t2*x))**(0.5))/(2*t2)
    return -x*np.log(K) + a*np.log(1+K*np.exp(beta*c)) + b*np.log(1+K*np.exp(beta*d)) + e

popt2, pconv2 = curve_fit(func, N_org, data1, p0=[100,100,-8,5,-100]) # p0=[100,100,-8,5,-100]
a2 = popt2[0]
b2 = popt2[1]
c2 = popt2[2]
d2 = popt2[3]
e2 = popt2[4]

print(a2,b2,c2,d2,e2)

out1 = open("dual-site-LI-parameters.txt","w")
out1.truncate(0)
print(-a2,file = out1)
print(-b2,file = out1)
print(-c2,file = out1)
print(-d2,file = out1)
print(-e2,file = out1)
out1.close()

t1m = (N-a2)*np.exp(beta*c2)+(N-b2)*np.exp(beta*d2)
t2m = (a2+b2-N)*np.exp(beta*(c2+d2))
Km = (t1m+((t1m**2)+(4*t2m*N))**(0.5))/(2*t2m)
q2 = -N*np.log(Km) + a2*np.log(1+Km*np.exp(beta*c2)) + b2*np.log(1+Km*np.exp(beta*d2)) + e2
"""

data1 = q2
length1 = len(data1)


def func2(x, a, b, c, d, e):
    return a + b * x + c * (x**2) + d * (x**3) + e * (x**4)


popt2, pconv2 = curve_fit(func2, N_org, datalnQ)
a2 = popt2[0]
b2 = popt2[1]
c2 = popt2[2]
d2 = popt2[3]
e2 = popt2[4]

print(a2, b2, c2, d2, e2)

out1 = open("fourth-order-polynomial-parameters.txt", "w")
out1.truncate(0)
print(-a2, file=out1)
print(-b2, file=out1)
print(-c2, file=out1)
print(-d2, file=out1)
print(-e2, file=out1)

out1.close()

q4th2 = a2 + b2 * N + c2 * (N**2) + d2 * (N**3) + e2 * (N**4)
data4th1 = q4th2
length1 = len(data4th1)

plt.figure(1)
plt.plot(N_org, datalnQ, "*", label="original data")
plt.plot(
    N_org[1:], func(N_org[1:], *popt1), "r*", label="curve_fit_single_site_LI_model"
)
plt.plot(N, q2, "r--", label="extrapolated data from single_site_model")
plt.plot(N_org[1:], func2(N_org[1:], *popt2), "b*", label="curve_fit_4th order fit")
plt.plot(N, q4th2, "b--", label="extrapolated data from 4th order fit")
plt.xlabel("N")
plt.ylabel("lnQ")
plt.legend()
plt.show()

out1 = open("ln-part-fn-cho-interpolated.txt", "w")
out1.truncate(0)
for i in range(len(data1)):
    print(data1[i], file=out1)
out1.close()

out1 = open("Total-Helmholtz-interpolated.txt", "w")
out1.truncate(0)
for i in range(len(data1) - 1):
    print(-data1[i + 1] / (i + 1), file=out1)
out1.close()

data3 = np.loadtxt("chem_pot._values_Pa.txt")
length3 = len(data3)
beta = beta * 1e20

lGc = []
for i in data3:
    r1 = 0
    r2 = 0
    print()
    # for j in range(length1):
    # r1 += np.exp(data1[j] + beta*i*N[j])
    # r2 += np.exp(data2[j] + beta*i*N[j])
    # Gc.append(r1)
    lGc.append(np.log(np.sum(np.exp(data1 + beta * i * N))))
    # Gq.append(r2)
# out1 = open("GCN-part-fn-cho.txt","w")
out2 = open("ln-GCN-part-fn-cho.txt", "w")
# out1.truncate(0)
out2.truncate(0)
for i in range(length3):
    # print(Gc[i],file = out1)
    # print(np.log(Gc[i]),file = out2)
    print(lGc[i], file=out2)
# out1.close()
out2.close()

"""
r1 = []
r2 = []
for i in [-3.25*1e-20,-3.15*1e-20,-3.05*1e-20,-2.95*1e-20,-2.75*1e-20]:
    r1.append(beta*i*N)  
    r2.append(data1+beta*i*N)

plt.figure(8)
#plt.plot(N,data1,"*",label = "lnQ")
plt.plot(N,r2[0],"r*",label = "-3.25*1e-20")
plt.plot(N,r2[1],"g*",label = "-3.15*1e-20")
plt.plot(N,r2[2],"b*",label = "-3.05*1e-20")
plt.plot(N,r2[3],"k*",label = "-2.95*1e-20")
plt.plot(N,r2[4],"m*",label = "-2.75*1e-20")
#plt.plot(N_org,func(N_org,*popt2),"k*",label = "actual_curve_fit_4th_order")
plt.xlabel("N")
plt.ylabel("lnQ+beta*mu*N")
plt.legend()


plt.figure(9)
plt.plot(data1,r1[0],"r*",label = "-3.25*1e-20")
plt.plot(data1,r1[1],"g*",label = "-3.15*1e-20")
plt.plot(data1,r1[2],"b*",label = "-3.05*1e-20")
plt.plot(data1,r1[3],"k*",label = "-2.95*1e-20")
plt.plot(data1,r1[4],"m*",label = "-2.75*1e-20")
#plt.plot(N_org,func(N_org,*popt2),"k--",label = "actual_curve_fit_4th_order")
plt.xlabel("lnQ")
plt.ylabel("beta*mu*N")
plt.legend()
"""

slope = np.log(b1 * (a1 - N) / N)
slope4th = b2 + 2 * c2 * N + 3 * d2 * (N**2) + 4 * e2 * (N**3)

"""
for i in range(length1):
    slope.append(b2 + 2*c2*N[i] + 3*d2*(N[i]**2) + 4*e2*(N[i]**3))
"""

# slp = spl1.derivative()
# slope = slp(N)
out1 = open("Slope-lnQ-vs-N-cho_new.txt", "w")
# out1 = open("Slope-lnQ-vs-N-cho.txt","w")
# out1 = open("Slope_lnQ_vs_N_qho.txt","w")
out1.truncate(0)
for i in range(length1):
    # print(f_prime(N[i]),file = out1)
    print(slope[i], file=out1)
out1.close()

plt.figure(2)
plt.plot(N, slope)
plt.xlabel("N")
plt.ylabel("lnQvsN_slope")

data = np.loadtxt("ln-GCN-part-fn-cho.txt")
length = len(data)

spl2 = InterpolatedUnivariateSpline(data3, data)
data2 = spl2(data3)
length = len(data2)
slp2 = spl2.derivative()
slope = slp2(data3)


# slope = []
data3mod = data3[0:] * 1e20


def func(x, a, b, c, d, e, f):
    return a + b * x + c * (x**2) + d * (x**3) + e * (x**4) + f * (x**5)


popt, pconv = curve_fit(func, data3mod, data)
a = popt[0]
b = popt[1]
c = popt[2]
d = popt[3]
e = popt[4]
f = popt[5]


# for i in range(len(data)):
#    slope.append((b + 2*c*data3mod[i] + 3*d*(data3mod[i]**2) + 4*e*(data3mod[i]**3) + 5*f*(data3mod[i]**4))*1e20)
# + 6*g*(data3mod[i]**5) g = popt[6] + g*(x**6) , g


plt.figure(4)
plt.plot(data3, data, "*", label="data")
# plt.plot(data3,func(data3mod,*popt),"r--",label = "curve_fit_5th_order")
plt.plot(data3, data2, "k--", label="spl fit")
plt.xlabel("chem.pot.")
plt.ylabel("lnGCN")
plt.legend()

datap = np.loadtxt("considered_pressure_values.txt")
plt.figure(6)
plt.plot(datap, data, "*", label="data")
# plt.plot(data3,func(data3mod,*popt),"k--",label = "curve_fit_5th_order")
plt.plot(datap, data2, "k--", label="spl fit")
plt.xlabel("fugacity in Pa")
plt.ylabel("lnGCN")
plt.legend()

out1 = open("Slope-lnGCNpf-vs-chempot-cho_newp.txt", "w")
out1.truncate(0)
for i in range(length):
    print(slope[i], file=out1)
out1.close()

plt.figure(5)
plt.plot(data3, slope, "*", label="data")
# plt.plot(data3,func(data3mod,*popt),"k--",label = "curve_fit_5th_order")
# plt.plot(data3,data2,"k--",label = "spl fit")
plt.xlabel("chem.pot.")
plt.ylabel("slope lnGCN vs chem. pot.")
plt.legend()

data4 = np.loadtxt("Slope-lnGCNpf-vs-chempot-cho_newp.txt")
length4 = len(data4)
for i in range(length4):
    a = data4[i] / beta
    Npre_cho.append(a)
    chepot_cho.append(data3[i])
    pressure_cho.append(np.exp((data3[i] - u0) / (kb * T)))


out = open("fugacities_from_chem_pot_in_Pa.txt", "w")
out.truncate(0)
for i in range(len(pressure_cho)):
    print(np.float(pressure_cho[i]), file=out)
out.close()


"""
#for methane
Tc = 190.4
Pc = 46*10**5
omega = 0.011
R = 8.314
kappa = 0.37464+1.54226*omega-0.26992*omega**2
Tr = T/Tc
alpha = (1+kappa*(1-Tr**0.5))**2
a = 0.45724*alpha*(R**2*Tc**2)/(Pc)
b = 0.07780*R*Tc/Pc
fugacity_coeff = []
pressure_actual = []
zep = 0
for i in pressure_cho:
    c1 = i
    c2 = b*i - R*T
    c3 = a - 2*b*R*T - 3*i*b**2
    c4 = i*b**3 - a*b + R*T*b**2
    coeffs = [c1, c2, c3, c4]
    d = np.roots(coeffs)
    e = float(d[np.isreal(d)].real)
    lnf = i*e/(R*T) - 1 + np.log((R*T/i)/(e-b)) - (a/(2*(2**(0.5))*R*T*b))*np.log((e+b+b*2**0.5)/(e+b-b*2**0.5))
    fugacity_coeff.append(np.exp(lnf))
    pressure_actual.append(i/fugacity_coeff[zep])
    zep = zep + 1

out = open("pressures_from_chem_pot_in_kPa.txt", "w")
out.truncate(0)
for i in range(len(pressure_actual)):
    print(np.float(pressure_actual[i]/1000), file=out)
out.close()
"""

out = open("N_predicted_from_GCN.txt", "w")
out.truncate(0)
for i in range(len(Npre_cho)):
    print(np.float(pressure_cho[i] / 1e3), np.float(Npre_cho[i]), file=out)
out.close()

plt.figure(1)
plt.plot(pressure_cho, np.array(Npre_cho), "g--", label="cho_GCN_spl_fit")
ax = plt.gca()
ax.set_xscale("log")
# plt.plot(pressure_cho,Npre_qho,label = "qho_GCN")
plt.xlabel("Pressure in SI units")
plt.ylabel("Loading")

datac1 = np.loadtxt(fname="Slope-lnQ-vs-N-cho_new.txt")
lengthc1 = len(datac1)
for i in range(lengthc1):
    u1.append(-kb * T * datac1[i])
    p1.append(np.exp((u1[i] - u0) / (kb * T)))

u4th = -kb * T * slope4th
p4th = np.exp((u4th - u0) / (kb * T))

"""
R = 8.314
fugacity_coeff = []
pressure_actual = []
zep = 0
for i in p1:
    c1 = i
    c2 = b*i - R*T
    c3 = a - 2*b*R*T - 3*i*b**2
    c4 = i*b**3 - a*b + R*T*b**2
    coeffs = [c1, c2, c3, c4]
    d = np.roots(coeffs)
    e = float(d[np.isreal(d)].real)
    lnf = i*e/(R*T) - 1 + np.log((R*T/i)/(e-b)) - (a/(2*(2**(0.5))*R*T*b))*np.log((e+b+b*2**0.5)/(e+b-b*2**0.5))
    fugacity_coeff.append(np.exp(lnf))
    pressure_actual.append(i/fugacity_coeff[zep])
    zep = zep + 1
"""
out = open("P_predicted_from_CN.txt", "w")
out.truncate(0)
for i in range(len(p1)):
    print(np.float(p1[i] / 1e3), np.float(N[i]), file=out)
out.close()

plt.figure(1)
plt.plot(p1[8:97:8], N[8:97:8], "r*", label="cho_CN_single_site_model")
plt.plot(p4th[8:97:8], N[8:97:8], "b*", label="cho_CN_4th_order_fit")
ax = plt.gca()
ax.set_xscale("log")

exp_pr = np.loadtxt(fname=exp_data, skiprows=1, usecols=0)
exp_N = np.loadtxt(fname=exp_data, skiprows=1, usecols=1)
exp_N = exp_N[0:] * 12
plt.figure(1)
ax = plt.gca()
ax.set_xscale("log")

LIl = np.linspace(0, int(a1), int(a1) + 1)
bpc = b1 * np.exp(u0 / (kb * T))

plt.plot(
    LIl / ((a1 - LIl) * bpc),
    LIl,
    label="Langmuir Isotherm fitted by calculated parameters",
)
plt.plot(exp_pr, exp_N / 12, label="mfi_data")
plt.xlim(5e3, 5e10)
plt.legend()
plt.show()


data1 = np.loadtxt(fname="ln-part-fn-cho-interpolated.txt")
data2 = np.loadtxt("chem_pot._values_Pa.txt")
l2 = len(data2)
u1 = data2[0]
u2 = data2[::-1][0]
j = np.linspace(u1, u2, 5)
FL = np.zeros((len(data2), 1))
for i in range(len(data2)):
    tot = np.exp(data1 + beta * data2[i] * N)
    FL[i] = N[int(np.float(np.where(tot == tot.max())[0]))]
out = open("N_predicted_from_FL.txt", "w")
out.truncate(0)
for i in range(len(FL)):
    print(np.float(np.exp((data3[i] - u0) / (kb * T)) / 1e3), np.float(FL[i]), file=out)
out.close()


data1 = np.loadtxt(fname="ln-part-fn-cho-interpolated.txt")
data2 = np.loadtxt("chem_pot._values_Pa.txt")
l2 = len(data2)
u1 = data2[0]
u2 = data2[::-1][0]
j = np.linspace(u1, u2, 5)
tot = np.zeros((97, 5))
for i in range(len(j)):
    tot[:, i] = (data1 + beta * j[i] * N) * (-1 / beta)
out = open("plot7_manuscript_data.txt", "w")
out.truncate(0)
for i in range(len(N)):
    print(
        np.float(N[i]),
        np.float(tot[i, 0]),
        np.float(tot[i, 1]),
        np.float(tot[i, 2]),
        np.float(tot[i, 3]),
        np.float(tot[i, 4]),
        file=out,
    )
out.close()
out = open("plot7_manuscript_data_scaled_by_kbT.txt", "w")
out.truncate(0)
for i in range(len(N)):
    print(
        np.float(N[i]),
        np.float(tot[i, 0] * beta),
        np.float(tot[i, 1] * beta),
        np.float(tot[i, 2] * beta),
        np.float(tot[i, 3] * beta),
        np.float(tot[i, 4] * beta),
        file=out,
    )
out.close()
