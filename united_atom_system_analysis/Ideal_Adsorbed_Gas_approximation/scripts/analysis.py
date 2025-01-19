import numpy as np
import matplotlib.pyplot as plt
import scipy.constants
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.optimize import curve_fit
import sys

exp_data = sys.argv[1]
GCN_cho = []
GCN_qho = []

kb = scipy.constants.k  # Boltzmann Constant k
Na = scipy.constants.Avogadro  # Avogadro's number
T = 308  # system temperature
beta = 1 / (kb * T)
h = scipy.constants.Planck  # Planck's constant
m = 16.04 * (1e-3) / Na  # mass of one methane molecule
# N_org = np.arange(0,190,10)
# N_org = np.arange(0,190,20)
N_org = np.arange(0, 370, 20)
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
N = np.linspace(0, 360, 361)
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
u0 = 0 - kb * T * np.log(r)
u11 = u0 + kb * T * np.log(1e4)
u21 = u0 + kb * T * np.log(2.5e6)
u = np.linspace(u11, u21, 1000)
out1 = open("chem_pot._values_1e4_3e6_Pa.txt", "w")
out2 = open("considered_pressure_values_1e4_3e6.txt", "w")
out1.truncate(0)
out2.truncate(0)
for i in range(len(u)):
    print(u[i], file=out1)
    print(np.exp((u[i] - u0) / (kb * T)), file=out2)
out1.close()
out2.close()
# ------------------------------------------------------------

datao1 = np.loadtxt("ln-part-fn-cho.txt")
# datao2 = np.loadtxt("Average_ln_Partition_Function_qho.txt")
data3 = np.loadtxt("chem_pot._values_1e4_3e6_Pa.txt")
length3 = len(data3)

spl1 = InterpolatedUnivariateSpline(N_org, datao1)
# spl2 = InterpolatedUnivariateSpline(N_org,datao2)
data1 = spl1(N)
# data2 = spl2(N)
length1 = len(data1)


def func(x, a, b, c, d, e):
    return a + b * x + c * (x**2) + d * (x**3) + e * (x**4)


popt, pconv = curve_fit(func, N, data1)
a = popt[0]
b = popt[1]
c = popt[2]
d = popt[3]
e = popt[4]

print(a, b, c, d, e)

popt2, pconv2 = curve_fit(func, N_org, datao1)
a2 = popt2[0]
b2 = popt2[1]
c2 = popt2[2]
d2 = popt2[3]
e2 = popt2[4]

print(a2, b2, c2, d2, e2)

q1 = a + b * N + c * (N**2) + d * (N**3) + e * (N**4)
q2 = a2 + b2 * N + c2 * (N**2) + d2 * (N**3) + e2 * (N**4)
data1 = q2
for i in range(len(N)):
    print(N[i], q1[i], q2[i], q1[i] - q2[i], (q1[i] - q2[i]) / q1[i])

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

plt.figure(3)
plt.plot(N_org, datao1, "b--", label="original data")
plt.plot(N, data1, "r", label="interpolated data")
plt.plot(N, func(N, *popt), "k--", label="spline_plus_curve_fit_4th_order")
plt.plot(N_org, func(N_org, *popt2), "k--", label="actual_curve_fit_4th_order")
plt.xlabel("N")
plt.ylabel("lnQ")
plt.legend()

plt.figure(6)
plt.plot(N, (q1 - q2) / q1)
plt.show()

"""
for i in range(len(datao1)):
    q1.append(datao1[i]/N_org[i])
splq1 = InterpolatedUnivariateSpline(N_org,q1)
q2 = splq1(N)

plt.figure(7)
plt.plot(N_org,q1,label = "lnQ/N vs N org")
plt.plot(N,q2,label = "lnQ/N vs N interpolated")
plt.legend()
"""
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
"""
out1 = open("GCN_partition_function_qho.txt","w")
out2 = open("ln_GCN_partition_function_qho.txt","w")
out1.truncate(0)
out2.truncate(0)
for i in range(length3):
    print(Gq[i],file = out1)
    print(np.log(Gq[i]),file = out2)
out1.close()
out2.close()
"""
"""
datau1 = np.loadtxt("GCN-part-fn-cho.txt")
#datau2 = np.loadtxt("GCN_partition_function_qho.txt")
for i in range(length3):
    r1 = 0
    r = []
    #r2 = 0
    for j in range(length1):
        r1 = r1 + (np.exp(data1[j]+beta*data3[i]*N[j]))*N[j]/datau1[i]
        #r2 = r2 + (np.exp(data2[j]+beta*data3[i]*N[j]))*N[j]/datau2[i]
        r.append((np.exp(data1[j]+beta*data3[i]*N[j]))*N[j]/datau1[i])
    Ndc.append(r1)
    #Ndq.append(r2)
p = np.exp((data3[0:]-u0)/(kb*T))

plt.figure(1)
ax = plt.gca()
plt.plot(p,Ndc,label = "dis_GCN_cho")
ax.set_xscale("log")
#plt.plot(p,Ndq,label = "dis_qho")
"""
"""
dataAo = np.loadtxt("Average-Helmholtz-cho.txt") 
dataAsc = dataAo[0:]*Na/1e3
out1 = open("Scaled-Helmholtz-cho.txt","w")
out1.truncate(0)
for i in range(len(dataAsc)):
    print(dataAsc[i],file = out1)
out1.close()
splA = InterpolatedUnivariateSpline(N_org,dataAo)
dataA = splA(N)
"""
"""
p = [1e4,1e5,1e6,1e7]
u = u0 + kb*T*np.log(p[0:])
for i in range(len(p)):
    r = []
    for j in range(len(dataA)):
        r.append(u[i]*N[j])
"""
"""
    plt.figure(2)
    plt.plot(N,r,label = "uN for P = %d"%p[i] + " Pa")
    plt.xlabel("N")
    plt.ylabel("A or uN")
plt.figure(2)
plt.plot(N,dataA,label = "A")
plt.legend()
"""
"""
for i in range(1,length1-1):
    s1 = (data1[i] - data1[i-1])/(N[i] - N[i-1])
    s2 = (data1[i+1] - data1[i])/(N[i+1] - N[i])
    s3 = (s1 + s2)/2
    slope.append(s3)
slope[0] = (data1[1] - data1[0])/(N[1] - N[0])
slope.append((data1[length1-1] - data1[length1-2])/(N[length1-1] - N[length1-2]))
"""
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

# fig2 = plt.figure(2)
# fig2.savefig("Slope_lnQ_vs_N.png")

"""
for i in range(1,len(datao1)-1):
    s1 = (datao1[i] - datao1[i-1])/(N_org[i] - N_org[i-1])
    s2 = (datao1[i+1] - datao1[i])/(N_org[i+1] - N_org[i])
    s3 = (s1 + s2)/2
    sl1.append(s3)
sl1[0] = (datao1[1] - datao1[0])/(N_org[1] - N_org[0])
sl1.append((datao1[len(datao1)-1] - datao1[len(datao1)-2])/(N_org[len(datao1)-1] - N_org[len(datao1)-2]))
"""
"""
plt.figure(4)
plt.plot(N,slope,label = "lnQ/N from interpolated lnQ vs N data")
plt.plot(N_org,sl1,label = "lnQ/N from original lnQ vs N data")
plt.xlabel("N")
plt.ylabel("lnQ/N")
plt.legend()
"""

data = np.loadtxt("ln-GCN-part-fn-cho.txt")
length = len(data)
slope = []
data3mod = data3[0:] * 1e20
"""
def func(x, a, b, c):
    return a + b*np.exp(c*x) 

popt , pconv = curve_fit(func,data3mod,data)
a = popt[0]
b = popt[1]
c = popt[2]

for i in range(len(data)):
    slope.append(b*c*1e20*np.exp(c*data3mod[i]))
"""


def func(x, a, b, c, d, e, f):
    return a + b * x + c * (x**2) + d * (x**3) + e * (x**4) + f * (x**5)


popt, pconv = curve_fit(func, data3mod, data)
a = popt[0]
b = popt[1]
c = popt[2]
d = popt[3]
e = popt[4]
f = popt[5]


for i in range(len(data)):
    slope.append(
        (
            b
            + 2 * c * data3mod[i]
            + 3 * d * (data3mod[i] ** 2)
            + 4 * e * (data3mod[i] ** 3)
            + 5 * f * (data3mod[i] ** 4)
        )
        * 1e20
    )
    # + 6*g*(data3mod[i]**5) g = popt[6] + g*(x**6) , g

plt.figure(4)
plt.plot(data3, data, "r--", label="data")
plt.plot(data3, func(data3mod, *popt), "k--", label="curve_fit_5th_order")
plt.xlabel("chem.pot.")
plt.ylabel("lnGCN")
plt.legend()

"""
for i in range(1,length-1):
    s1 = (data[i] - data[i-1])/(data3[i] - data3[i-1])
    s2 = (data[i+1] - data[i])/(data3[i+1] - data3[i])
    s3 = (s1 + s2)/2
    slope.append(s3)
slope[0] = (data[1] - data[0])/(data3[1] - data3[0])
slope.append((data[length-1] - data[length-2])/(data3[length-1] - data3[length-2]))
out1 = open("Slope-lnGCNpf-vs-chempot-cho.txt","w")
out1.truncate(0)
for i in range(length):
    print(slope[i],file = out1)
out1.close()
"""
out1 = open("Slope-lnGCNpf-vs-chempot-cho_newp.txt", "w")
out1.truncate(0)
for i in range(length):
    print(slope[i], file=out1)
out1.close()
"""
plt.figure(3)
plt.plot(N,slope)
plt.xlabel("N")
plt.ylabel("lnGCN_vs_N_slope")
fig2 = plt.figure(2)
fig2.savefig("Slope_lnGCN_vs_N.png")
"""
data4 = np.loadtxt("Slope-lnGCNpf-vs-chempot-cho_newp.txt")
# data5 = np.loadtxt("Slope_lnGCNpf_vs_chempot_qho.txt")
length4 = len(data4)
for i in range(length4):
    a = data4[i] / beta
    Npre_cho.append(a)
    chepot_cho.append(data3[i])
    pressure_cho.append(np.exp((data3[i] - u0) / (kb * T)))
"""
for i in range(length4):
    a = data5[i]/beta
    Npre_qho.append(a)
    chepot_qho.append(data3[i])
    pressure_qho.append(np.exp((data3[i]-u0)/(kb*T)))
"""

out = open("fugacities_from_chem_pot_in_Pa.txt", "w")
out.truncate(0)
for i in range(len(pressure_cho)):
    print(np.float(pressure_cho[i]), file=out)
out.close()

# for methane
Tc = 190.4
Pc = 46 * 10**5
omega = 0.011
R = 8.314
kappa = 0.37464 + 1.54226 * omega - 0.26992 * omega**2
Tr = T / Tc
alpha = (1 + kappa * (1 - Tr**0.5)) ** 2
a = 0.45724 * alpha * (R**2 * Tc**2) / (Pc)
b = 0.07780 * R * Tc / Pc
fugacity_coeff = []
pressure_actual = []
zep = 0
for i in pressure_cho:
    c1 = i
    c2 = b * i - R * T
    c3 = a - 2 * b * R * T - 3 * i * b**2
    c4 = i * b**3 - a * b + R * T * b**2
    coeffs = [c1, c2, c3, c4]
    d = np.roots(coeffs)
    e = float(d[np.isreal(d)].real)
    lnf = (
        i * e / (R * T)
        - 1
        + np.log((R * T / i) / (e - b))
        - (a / (2 * (2 ** (0.5)) * R * T * b))
        * np.log((e + b + b * 2**0.5) / (e + b - b * 2**0.5))
    )
    fugacity_coeff.append(np.exp(lnf))
    pressure_actual.append(i / fugacity_coeff[zep])
    zep = zep + 1

out = open("pressures_from_chem_pot_in_kPa.txt", "w")
out.truncate(0)
for i in range(len(pressure_cho)):
    print(np.float(pressure_actual[i] / 1000), file=out)
out.close()

out = open("N_predicted_from_GCN.txt", "w")
out.truncate(0)
for i in range(len(Npre_cho)):
    print(np.float(pressure_cho[i] / 1e3), np.float(Npre_cho[i] / 24), file=out)
out.close()

plt.figure(1)
plt.plot(pressure_cho, Npre_cho, "b--", label="cho_GCN_5th-order_poly_fitted")
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

fugacity_coeff = []
pressure_actual = []
zep = 0
for i in p1:
    c1 = i
    c2 = b * i - R * T
    c3 = a - 2 * b * R * T - 3 * i * b**2
    c4 = i * b**3 - a * b + R * T * b**2
    coeffs = [c1, c2, c3, c4]
    d = np.roots(coeffs)
    e = float(d[np.isreal(d)].real)
    lnf = (
        i * e / (R * T)
        - 1
        + np.log((R * T / i) / (e - b))
        - (a / (2 * (2 ** (0.5)) * R * T * b))
        * np.log((e + b + b * 2**0.5) / (e + b - b * 2**0.5))
    )
    fugacity_coeff.append(np.exp(lnf))
    pressure_actual.append(i / fugacity_coeff[zep])
    zep = zep + 1

out = open("P_predicted_from_CN.txt", "w")
out.truncate(0)
for i in range(len(p1)):
    print(np.float(p1[i] / 1e3), np.float(N[i] / 24), file=out)
out.close()

plt.figure(1)
plt.plot(p1, N, label="cho_CN_4th-order_poly_fitted")
ax.set_xscale("log")
# plt.plot(p2,N, label = "qho_CN")
"""
datar3 = np.loadtxt(fname = "Average-Entropy-cho.txt")
datar3 = datar3[0:]*kb
out1 = open("Scaled-Entropy-cho.txt","w")
out1.truncate(0)
for i in range(len(datar3)):
    print(datar3[i]*N_org[i]*Na*12/1e3,file = out1)
out1.close()
"""
"""
c = 0
data1 = np.loadtxt(fname = "chem_pot._values_1e4_3e6_Pa.txt")
datar2 = np.loadtxt(fname = "Average_MD_Energy.txt")
datar3 = np.loadtxt(fname = "Average_Entropy_cho.txt")
#datar4 = np.loadtxt(fname = "Average_Entropy_qho.txt")

for i in range(len(datar2)):
    datar3[i] = datar3[i]*N_org[i]
    #datar4[i] = datar4[i]*N_org[i]

spl2 = InterpolatedUnivariateSpline(N_org,datar2)
spl3 = InterpolatedUnivariateSpline(N_org,datar3)
#spl4 = InterpolatedUnivariateSpline(N_org,datar4)
#data4 = spl4(N)
data2 = spl2(N)
data3 = spl3(N)
length1 = len(data1)
length3 = len(data2)
for j in range(length1):
    for i in range(length3):
        a1.append(data2[i] - T*data3[i] - data1[j]*N[i])
        #a2.append(data2[i] - T*data4[i] - data1[j]*N[i])
    if(j==((int)((length1)/10)*c)):
        
        plt.figure(2)
        plt.plot(N,a1,label = np.exp((data1[j]-u0)/(kb*T)))
        plt.xlabel("N")
        plt.ylabel("A-uN")
        plt.legend(loc = 'upper left')
        plt.figure(3)
        plt.plot(N,a2,label = np.exp((data1[j]-u0)/(kb*T)))
        plt.xlabel("N")
        plt.ylabel("A-uN")
        plt.legend(loc = 'upper left')  
        
        c = c + 1
    min1.append(N[a1.index(min(a1))])
    #min2.append(N[a2.index(min(a2))])
    a1 = []
    #a2 = []
plt.legend(loc = 'upper left')
p1 = np.exp((data1[0:]-u0)/(kb*T))
"""
exp_pr = np.loadtxt(fname=exp_data, skiprows=1, usecols=0)
exp_N = np.loadtxt(fname=exp_data, skiprows=1, usecols=1)
exp_N = exp_N[0:] * 12
plt.figure(1)
ax = plt.gca()
ax.set_xscale("log")
# plt.plot(p1,min1, label = "cho_A-uN_up")
# plt.plot(p1,min2, label = "qho_A-uN")
plt.plot(exp_pr, exp_N, label="mfi_data")
plt.xlim(5e3, 3e6)
plt.legend()
plt.show()
# fig = plt.figure(1)
# fig.savefig("loading_vs_P.png")

"""
gas_helm = np.loadtxt(fname = "Average-Gas-Helmholtz-cho.txt")
sol_helm = np.loadtxt(fname = "Average-Solid-Helmholtz-cho.txt")
"""
"""
gas_int = InterpolatedUnivariateSpline(N_org,gas_helm)
sol_int = InterpolatedUnivariateSpline(N_org,sol_helm)
gas_fin = gas_int(N)
sol_fin = sol_int(N)
"""
"""
def func_gas(x, a, b, c, d, e):
    return a + b*x + c*(x**2) + d*(x**3) + e*(x**4) 

popt , pconv = curve_fit(func_gas,N,gas_fin)
a = popt[0]
b = popt[1]
c = popt[2]
d = popt[3]
e = popt[4]

for i in range(len(gas_helm)):
    slope.append((b + 2*c*N_org[i] + 3*d*(N_org[i]**2) + 4*e*(N_org[i]**3)))

def func_sol(x, p, q, r, s, t):
    return p + q*x + r*(x**2) + s*(x**3) + t*(x**4) 

popt , pconv = curve_fit(func_gas,N,gas_fin)
p = popt[0]
q = popt[1]
r = popt[2]
s = popt[3]
t = popt[4]

for i in range(len(gas_helm)):
    slope.append((q + 2*r*N_org[i] + 3*s*(N_org[i]**2) + 4*t*(N_org[i]**3)))

u_ex = -8.66171*(1e3)/Na # excess gibbs free energy
for i in range
r = (((2*np.pi*m*kb*T)/(h**2))**(1.5))*(vol/1)
u_ig = -kb*T*np.log(r) # ideal gibbs free energy
u = u_ig + u_ex # total gibbs free energy
s = -(u - enermd)/T
print(u_ig)
"""

data1 = np.loadtxt(fname="ln-part-fn-cho-interpolated.txt")
data2 = np.loadtxt("chem_pot._values_1e4_3e6_Pa.txt")
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
    print(
        np.float(np.exp((data3[i] - u0) / (kb * T)) / 1e3),
        np.float(FL[i] / 24),
        file=out,
    )
out.close()


data1 = np.loadtxt(fname="ln-part-fn-cho-interpolated.txt")
data2 = np.loadtxt("chem_pot._values_1e4_3e6_Pa.txt")
l2 = len(data2)
u1 = data2[0]
u2 = data2[::-1][0]
j = np.linspace(u1, u2, 5)
tot = np.zeros((361, 5))
for i in range(len(j)):
    tot[:, i] = (data1 + beta * j[i] * N) * (-1 / beta)
out = open("plot7_manuscript_data.txt", "w")
out.truncate(0)
for i in range(len(N)):
    print(
        np.float(N[i] / 24),
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
        np.float(N[i] / 24),
        np.float(tot[i, 0] * beta),
        np.float(tot[i, 1] * beta),
        np.float(tot[i, 2] * beta),
        np.float(tot[i, 3] * beta),
        np.float(tot[i, 4] * beta),
        file=out,
    )
out.close()
