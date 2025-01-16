import numpy as np
import matplotlib.pyplot as plt
import scipy.constants
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.optimize import curve_fit
import sys

exp_data = sys.argv[1]
GCN_cho = []
GCN_qho = []

kb = scipy.constants.k # Boltzmann Constant k
Na = scipy.constants.Avogadro # Avogadro's number
T = 87 # system temperature
beta = 1/(kb*T)
h = scipy.constants.Planck # Planck's constant
m = 39.948*(1e-3)/Na # mass of one methane molecule
N_org = np.arange(0,1100,50)
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
#slope.append(0)
sl1 = []
sl1.append(0)
N = np.linspace(0,1050,1051)
N = N[0:]
P_cho = []
P_qho = []
Gc = []
lGc = []
Gq = []
Ndc = []
Ndq = []
q1 = []

# chemical potential
r = (((2*np.pi*m*kb*T)/(h**2))**(1.5))*kb*T
u0 = 0 - kb*T*np.log(r)
u11 = (u0 + kb*T*np.log(1e-2))
u21 = (u0 + kb*T*np.log(4e2))
u = np.linspace(u11,u21,1000)
out1 = open("chem_pot._values_Pa.txt","w")
out2 = open("considered_pressure_values.txt","w")
out1.truncate(0)
out2.truncate(0)
for i in range(len(u)):
    print(u[i], file=out1)
    print(np.exp((u[i]-u0)/(kb*T)), file=out2)
out1.close()
out2.close()
#------------------------------------------------------------

datao1 = np.loadtxt("ln-part-fn-cho.txt")
data3 = np.loadtxt("chem_pot._values_Pa.txt")
length3 = len(data3)

spl1 = InterpolatedUnivariateSpline(N_org,datao1)
data1 = spl1(N)
length1 = len(data1)

"""
def func(x, a, b, c, d, e, f, g):
    return a + b*x + c*(x**2) +d*(x**3) + e*(x**4) + f*(x**5) + g*(x**6)

popt2, pconv2 = curve_fit(func, N_org, datao1)
a2 = popt2[0]
b2 = popt2[1]
c2 = popt2[2]
d2 = popt2[3]
e2 = popt2[4]
f2 = popt2[5]
g2 = popt2[6]

print(a2,b2,c2,d2,e2,f2,g2)

out1 = open("fourth-order-polynomial-parameters.txt","w")
out1.truncate(0)
print(-a2,file = out1)
print(-b2,file = out1)
print(-c2,file = out1)
print(-d2,file = out1)
print(-e2,file = out1)
print(-f2,file = out1)
print(-g2,file = out1)
out1.close()

q2 = a2 + b2*N + c2*(N**2) +d2*(N**3) + e2*(N**4) + f2*(N**5) + g2*(N**6) 
data1 = q2
length1 = len(data1)
"""

out1 = open("ln-part-fn-cho-interpolated.txt","w")
out1.truncate(0)
for i in range(len(data1)):
    print(data1[i],file = out1)
out1.close()

out1 = open("Total-Helmholtz-interpolated.txt","w")
out1.truncate(0)
for i in range(len(data1)-1):
    print(-data1[i+1]/(i+1),file = out1)
out1.close()

plt.figure(3)
plt.plot(N_org,datao1,"b--",label = "original data")
plt.plot(N,data1,"r",label = "interpolated data")
#plt.plot(N_org,func(N_org,*popt2),"k--",label = "actual_curve_fit_4th_order")
plt.xlabel("N")
plt.ylabel("lnQ")
plt.legend()

for i in data3:
    r1 = 0  
    r2 = 0
    print()
    #for j in range(length1):
        #r1 += np.exp(data1[j] + beta*i*N[j])
        #r2 += np.exp(data2[j] + beta*i*N[j])
    #Gc.append(r1)
    lGc.append(np.log(np.sum(np.exp(data1 + beta*i*N))))
    #Gq.append(r2)
#out1 = open("GCN-part-fn-cho.txt","w")
out2 = open("ln-GCN-part-fn-cho.txt","w")
#out1.truncate(0)
out2.truncate(0)
for i in range(length3):
    #print(Gc[i],file = out1)
    print(lGc[i],file = out2)
#out1.close()
out2.close()


for i in range(length1):
    slope.append(b2 + 2*c2*N[i] + 3*d2*(N[i]**2) + 4*e2*(N[i]**3) + 5*f2*(N[i]**4) + 6*g2*(N[i]**5))

#slp = spl1.derivative()
#slope = slp(N)
out1 = open("Slope-lnQ-vs-N-cho_new.txt","w")
#out1 = open("Slope-lnQ-vs-N-cho.txt","w")
#out1 = open("Slope_lnQ_vs_N_qho.txt","w")
out1.truncate(0)
for i in range(length1):
    print(slope[i],file = out1)
out1.close()

plt.figure(2)
plt.plot(N,slope)
plt.xlabel("N")
plt.ylabel("lnQvsN_slope")

data = np.loadtxt("ln-GCN-part-fn-cho.txt")
length = len(data) 
slope = []
data3mod = data3[0:]*1e20

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
length4 = len(data4)
for i in range(length4):
    a = data4[i]/beta
    Npre_cho.append(a)
    chepot_cho.append(data3[i])
    pressure_cho.append(np.exp((data3[i]-u0)/(kb*T)))


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
    print(np.float(pressure_cho[i]/1e3),np.float(Npre_cho[i]/8), file=out)
out.close()

plt.figure(1)
plt.plot(pressure_cho,Npre_cho,"b--",label = "cho_GCN_5th-order_poly_fitted")
ax = plt.gca()
ax.set_xscale("log")
#plt.plot(pressure_cho,Npre_qho,label = "qho_GCN")
plt.xlabel("Pressure in SI units")
plt.ylabel("Loading")

datac1 = np.loadtxt(fname = "Slope-lnQ-vs-N-cho_new.txt")
lengthc1 = len(datac1)
for i in range(lengthc1):
    u1.append(-kb*T*datac1[i])
    p1.append(np.exp((u1[i]-u0)/(kb*T)))

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
    print(np.float(p1[i]/1e3),np.float(N[i]/8), file=out)
out.close()

plt.figure(1)
plt.plot(p1,N, label = "cho_CN_4th-order_poly_fitted")
ax = plt.gca()
ax.set_xscale("log")

exp_pr = np.loadtxt(fname = exp_data,skiprows = 1,usecols = 0)
exp_N = np.loadtxt(fname = exp_data,skiprows = 1,usecols = 1)
exp_N = exp_N[0:]*12
plt.figure(1)
ax = plt.gca()
ax.set_xscale("log")

plt.plot(exp_pr,exp_N,label = "mfi_data")
plt.xlim(5e3,3e6)
plt.legend()
plt.show()


data1 = np.loadtxt(fname="ln-part-fn-cho-interpolated.txt")
data2 = np.loadtxt("chem_pot._values_Pa.txt")
l2 = len(data2)
u1 = data2[0]
u2 = data2[::-1][0]
j = np.linspace(u1,u2,5)
FL = np.zeros((len(data2),1))
for i in range(len(data2)):
    tot = np.exp(data1 + beta*data2[i]*N)
    FL[i] = N[int(np.float(np.where(tot==tot.max())[0]))]
out = open("N_predicted_from_FL.txt", "w")
out.truncate(0)
for i in range(len(FL)):
    print(np.float(np.exp((data3[i]-u0)/(kb*T))/1e3),np.float(FL[i]/8), file=out)
out.close()


data1 = np.loadtxt(fname="ln-part-fn-cho-interpolated.txt") 
data2 = np.loadtxt("chem_pot._values_Pa.txt") 
l2 = len(data2) 
u1 = data2[0] 
u2 = data2[::-1][0] 
j = np.linspace(u1,u2,5) 
tot = np.zeros((1051,5)) 
for i in range(len(j)): 
   tot[:,i] = (data1 + beta*j[i]*N)*(-1/beta)
out = open("plot7_manuscript_data.txt", "w") 
out.truncate(0) 
for i in range(len(N)): 
   print(np.float(N[i]/8),np.float(tot[i,0]),np.float(tot[i,1]),np.float(tot[i,2]),np.float(tot[i,3]),np.float(tot[i,4]), file=out) 
out.close()
out = open("plot7_manuscript_data_scaled_by_kbT.txt", "w") 
out.truncate(0) 
for i in range(len(N)): 
   print(np.float(N[i]/8),np.float(tot[i,0]*beta),np.float(tot[i,1]*beta),np.float(tot[i,2]*beta),np.float(tot[i,3]*beta),np.float(tot[i,4]*beta), file=out) 
out.close()

