import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
import sys
import math
import scipy.constants
from scipy.fftpack import fft, ifft, fftfreq
#import MDAnalysis as mda
import sys
from scipy import integrate
from scipy import optimize
from time import perf_counter 
from scipy.optimize import curve_fit
#from numba import jit
import psutil

# Defining Consatants
kb = scipy.constants.k # Boltzmann Constant k
h = scipy.constants.Planck # Planck's constant
Na = scipy.constants.Avogadro # Avogadro's number
T = 300 # Temperature, will make sure to make this as user input
beta = 1/(kb*T) # beta
tot_mass = 18.01528024673462*(1e-3)/Na

# rotational fluidicity
def rot_hs_fluidicity(m, den, y0, T):
    delta = math.sqrt((np.pi*kb*T)/m)*((6/np.pi)**(2/3))*(den**(1/3))*(2*y0)/9 # normalized diffusivity constant
    def fluidicity(f):
        return (2*(delta**(0-9/2))*(f**(15/2))-(delta**(0-3/2))*(f**(7/2))-6*(f**(5))*(delta**(0-3))+2*f+6*(delta**(0-3/2))*(f**(5/2))-2)
    f_value = optimize.root_scalar(fluidicity,bracket = [0,1],method = 'brentq')
    frac = f_value.root
    return frac

# rotational dos
def ref_rotational_dos(fl, tdos, y, ff):
    p = ((np.pi*y[0]*ff)/(2*fl*tdos)).reshape(-1,1)
    yg = (y[0]/(1+(p**2))).reshape(-1,1) # calculating gas DOS function vs frequency
    ys = y-yg # calculating solid DOS function vs frequency
    return yg, ys

# hard sphere rotational entropy
def hs_rot_entropy_wf(tot_mass, den, rot_fl, k):
    r = (((((np.pi)*(np.e)**3)*(T**3))/(13.706211107457026*20.998505222742857*39.467682045442515))**(0.5))*0.5
    si = kb*np.log(r)  # entropy value per molecule basis
    return (si/(3*kb))

# total rotational entropy
def rot_entropy(W_Shs, ff, y, ys, yg):
    ygm = W_Shs*yg
    W_Scho = (1 - np.log(beta*h*ff[1:len(ff)]))
    outd1 = open("cho_ent_wf.txt", "w")
    outd1.truncate(0)
    for i in range(len(W_Scho)): 
        print(np.float(ff[i+1]*1e-12),np.float(W_Scho[i]), file = outd1) 
    outd1.close() 
    ysm1 = W_Scho*ys[1:len(ys)]
    outd1 = open("cho_ent_wf_times_sol_dos.txt", "w")
    outd1.truncate(0)
    for i in range(len(ysm1)): 
        print(np.float(ff[i+1]*1e-12),np.float(ysm1[i]), file = outd1) 
    outd1.close() 
    a1 = (beta*h*ff[1:len(ff)])/(np.exp(beta*h*ff[1:len(ff)])-1)
    a2 = np.log(1 - np.exp(0 - beta*h*ff[1:len(ff)]))
    W_Sqho = a1-a2
    outd1 = open("qho_ent_wf.txt", "w")
    outd1.truncate(0)
    for i in range(len(W_Sqho)): 
        print(np.float(ff[i+1]*1e-12),np.float(W_Sqho[i]), file = outd1) 
    outd1.close() 
    ysm2 = W_Sqho*ys[1:len(ff)]
    outd1 = open("qho_ent_wf_times_sol_dos.txt", "w")
    outd1.truncate(0)
    for i in range(len(ysm2)): 
        print(np.float(ff[i+1]*1e-12),np.float(ysm2[i]), file = outd1) 
    outd1.close()     
    S_gas = (integrate.simps(ygm[1:len(ygm)].T,ff[1:len(ygm)].T))*kb # Gas Entropy
    S_cho1 = (integrate.simps(ysm1.T,ff[1:len(ff)].T))*kb # Solid Entropy (assuming classical harmonic oscillator)
    S_cho =  S_cho1 + S_gas # Entropy (assuming classical harmonic oscillator)
    S_qho1 = (integrate.simps(ysm2.T,ff[1:len(ff)].T))*kb # Solid Entropy (assuming quantum harmonic oscillator)
    S_qho =  S_qho1 + S_gas # Entropy (assuming quantum harmonic oscillator)
    print(np.float(S_gas))
    print(np.float(S_cho1))
    print(np.float(S_qho1))
    return np.float(S_cho), np.float(S_qho), np.float(S_gas), np.float(S_cho1), np.float(S_qho1)

# rotational internal energy (assumption --> reference energy = 0)
def rot_energy(W_Ehs, ff, y, ys, yg):
    b1 = beta*h*ff[1:len(ff)]/2
    b2 = (beta*h*ff[1:len(ff)])/(np.exp(beta*h*ff[1:len(ff)])-1)
    W_Eqho = b1+b2
    outd1 = open("qho_ener_wf.txt", "w")
    outd1.truncate(0)
    for i in range(len(W_Eqho)): 
        print(np.float(ff[i+1]*1e-12),np.float(W_Eqho[i]), file = outd1) 
    outd1.close()
    ysm2 = W_Eqho*ys[1:len(ff)]
    outd1 = open("qho_energy_wf_times_sol_dos.txt", "w")
    outd1.truncate(0)
    for i in range(len(ysm2)): 
        print(np.float(ff[i+1]*1e-12),np.float(ysm2[i]), file = outd1) 
    outd1.close()    
    E_cho = (integrate.simps(yg[1:len(yg)].T*W_Ehs,ff[1:len(yg)].T) + integrate.simps(ys[1:len(ys)].T,ff[1:len(ys)].T))*(beta**(-1)) # Energy (assuming classical harmonic oscillator)
    E_qho = (integrate.simps(yg[1:len(yg)].T*W_Ehs,ff[1:len(yg)].T) + integrate.simps(ysm2[1:len(ysm2)].T,ff[1:len(ysm2)].T))*(beta**(-1)) # Energy (assuming classical harmonic oscillator) 
    E_gas = np.float(integrate.simps(yg[1:len(yg)].T*W_Ehs,ff[1:len(yg)].T)*(beta**(-1)))
    E_cho1 = np.float(integrate.simps(ys[1:len(ys)].T,ff[1:len(ys)].T)*(beta**(-1)))
    E_qho1 = np.float(integrate.simps(ysm2[1:len(ysm2)].T,ff[1:len(ysm2)].T)*(beta**(-1)))
    print(E_gas)
    print(E_cho1)
    print(E_qho1)
    return np.float(E_cho), np.float(E_qho), np.float(E_gas), np.float(E_cho1), np.float(E_qho1)

# rotational helmholtz energy (assumption --> reference energy = 0)
def rot_helmholtz(W_Ahs, ff, y, ys, yg):
    W_Acho = np.log(beta*h*ff[1:len(ff)])
    outd1 = open("cho_helm_wf.txt", "w")
    outd1.truncate(0)
    for i in range(len(W_Acho)): 
        print(np.float(ff[i+1]*1e-12),np.float(W_Acho[i]), file = outd1) 
    outd1.close()  
    ysm1 = W_Acho*ys[1:len(ff)]
    outd1 = open("cho_helm_wf_times_sol_dos.txt", "w")
    outd1.truncate(0)
    for i in range(len(ysm1)): 
        print(np.float(ff[i+1]*1e-12),np.float(ysm1[i]), file = outd1) 
    outd1.close()    
    W_Aqho = np.log((1 - np.exp(0 - beta*h*ff[1:len(ff)]))/np.exp(0 - beta*h*ff[1:len(ff)]/2))
    outd1 = open("qho_helm_wf.txt", "w")
    outd1.truncate(0)
    for i in range(len(W_Aqho)): 
        print(np.float(ff[i+1]*1e-12),np.float(W_Aqho[i]), file = outd1) 
    outd1.close()     
    ysm2 = W_Aqho*ys[1:len(ff)]
    outd1 = open("qho_helm_wf_times_sol_dos.txt", "w")
    outd1.truncate(0)
    for i in range(len(ysm2)): 
        print(np.float(ff[i+1]*1e-12),np.float(ysm2[i]), file = outd1) 
    outd1.close()   
    ygm = W_Ahs*yg
    A_cho = (integrate.simps(ygm[1:len(ygm)].T,ff[1:len(ygm)].T) + (integrate.simps(ysm1.T,ff[1:len(ff)].T)))*(beta**(-1)) # Helmholtz free energy (assuming classical harmonic oscillator)
    A_qho = (integrate.simps(ygm[1:len(ygm)].T,ff[1:len(ygm)].T) + (integrate.simps(ysm2.T,ff[1:len(ff)].T)))*(beta**(-1)) # Helmholtz free energy (assuming quantum harmonic oscillator)
    A_gas = np.float(integrate.simps(ygm[1:len(ygm)].T,ff[1:len(ygm)].T)*(beta**(-1)))
    A_cho1 = np.float(integrate.simps(ysm1[1:len(ysm1)].T,ff[1:len(ysm1)].T)*(beta**(-1)))
    A_qho1 = np.float(integrate.simps(ysm2[1:len(ysm2)].T,ff[1:len(ysm2)].T)*(beta**(-1)))
    print(A_gas)
    print(A_cho1)
    print(A_qho1)
    return np.float(A_cho), np.float(A_qho), np.float(A_gas), np.float(A_cho1), np.float(A_qho1)


def trans_iag_helmholtz(vol, m, T):
    A_ig = []
    Nr = np.arange(0,200,1)
    for i in Nr:
        x = -kb*T*i*((np.log(vol*((2*np.pi*m*kb*T)/(h**2))**(1.5)))) #+ np.log(*(((np.pi*(T**3))/(13.706211107457026*20.998505222742857*39.467682045442515))**(0.5))*0.5))
        for j in range(i):
            x = x + kb*T*np.log(j+1)
        A_ig.append(x)
        x = 0

    def func(x, a, b, c, d, e):
        return a + b*x + c*(x**2) + d*(x**3) + e*(x**4) 

    popt , pconv = curve_fit(func,Nr,A_ig)
    a = popt[0]
    b = popt[1]
    c = popt[2]
    d = popt[3]
    e = popt[4]

    return a, b, c, d, e

#--------------------------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    time_start = perf_counter() # starting time counter

    outfile1 = sys.argv[1]
    outfile2 = sys.argv[2]
    enerfile = sys.argv[3]
   #runnum = sys.argv[4]

    print("files read")

    """    
    tot_dos = np.array(np.loadtxt(fname = "tot_dos_for_5fs_dt_for_"+runnum+".txt", usecols=1)) # total density of states
    tot_dos = tot_dos.reshape(-1,1)
    tr_dos = np.array(np.loadtxt(fname = "tr_dos_for_5fs_dt_for_"+runnum+".txt", usecols=1)) # total density of states
    tr_dos = tr_dos.reshape(-1,1)
    rot_dos = np.array(np.loadtxt(fname = "rot_dos_for_5fs_dt_for_"+runnum+".txt", usecols=1)) # total density of states
    rot_dos = rot_dos.reshape(-1,1)
    ff = np.array(np.loadtxt(fname = "rot_dos_for_5fs_dt_for_"+runnum+".txt", usecols=0)) # total density of states
    ff = ff.reshape(-1,1)   
    """    
    tot_dos = np.array(np.loadtxt(fname = "avg_tot_dos_for_5fs_dt.txt", usecols=1)) # total density of states
    tot_dos = tot_dos.reshape(-1,1)
    tr_dos = np.array(np.loadtxt(fname = "avg_tr_dos_for_5fs_dt.txt", usecols=1)) # total density of states
    tr_dos = tr_dos.reshape(-1,1)
    rot_dos = np.array(np.loadtxt(fname = "avg_rot_dos_for_5fs_dt.txt", usecols=1)) # total density of states
    rot_dos = rot_dos.reshape(-1,1)
    ff = np.array(np.loadtxt(fname = "avg_rot_dos_for_5fs_dt.txt", usecols=0)) # total density of states
    ff = ff.reshape(-1,1)
    tot_dos_val = integrate.simps(tot_dos.T,ff[:len(tot_dos)].T)
    tr_dos_val = integrate.simps(tr_dos.T,ff[:len(tr_dos)].T)
    rot_dos_val = integrate.simps(rot_dos.T,ff[:len(rot_dos)].T)
    print(np.float(tot_dos_val))
    print(np.float(tr_dos_val))
    print(np.float(rot_dos_val)) 

    tr_dos = tr_dos*(1e-12)/(tr_dos_val/3)
    rot_dos = rot_dos*(1e-12)/(rot_dos_val/3)
    tot_dos = tot_dos*(1e-12)/(tot_dos_val/6)
    ff = ff*(1e12)
    outd1 = open("norm_rot_dos.txt", "w")
    outd1.truncate(0)
    for i in range(len(rot_dos)): 
        print(np.float(ff[i]*1e-12),np.float(rot_dos[i]*1e12), file = outd1) 
    outd1.close()  
    tot_dos_val = integrate.simps(tot_dos.T,ff[:len(tot_dos)].T)
    tr_dos_val = integrate.simps(tr_dos.T,ff[:len(tr_dos)].T)
    rot_dos_val = integrate.simps(rot_dos.T,ff[:len(rot_dos)].T)
    print(np.float(tot_dos_val))
    print(np.float(tr_dos_val))
    print(np.float(rot_dos_val)) 

    # copying translational part dos to a file
    out = open(outfile1,"w")
    out.truncate(0)
    for i in range(len(tr_dos)):
        print(np.float(tr_dos[i]), file = out)
    out.close()

    vol = 4.00440*3.97980*4.01490*(1e-27) # volume of the system
    den = 1/vol # calculating density

   # Total system data
    md_energy = (-244.025/14.3771)*1e3/Na + 3*kb*T # total md energy
    #md_energy = (np.float(enerfile)*1e3/Na + 750269*1e3/Na)/100 # total md energy*1e3/Na + 750269*1e3/Na
    u_ex = -6.64901*(1e3)/Na # excess gibbs free energy
    r = ((((2*np.pi*tot_mass*kb*T)/(h**2))**(1.5))*(vol/1))*(((np.pi*(T**3))/(13.706211107457026*20.998505222742857*39.467682045442515))**(0.5))*0.5 
    u_ig = -kb*T*np.log(r) # ideal gibbs free energy
    u = u_ig + u_ex # total gibbs free energy
    s = -(u - md_energy)/T # total entropy
    r_rot = (((np.pi*(T**3))/(13.706211107457026*20.998505222742857*39.467682045442515))**(0.5))*0.5 
    u_ig_rot = -kb*T*np.log(r_rot)
    
    # Rotational system data
    rot_dif_sp = np.float(rot_dos[0]*(T*kb)/(12*tot_mass*1)) # Rotational Self diffusion coefficient (p)
    rot_fl = rot_hs_fluidicity(tot_mass, den, rot_dos[0], T) # Rotational fluidicity
    print(rot_fl)
    gas_rot_dos, solid_rot_dos = ref_rotational_dos(rot_fl, rot_dos_val, rot_dos, ff) # gas and solid rotational dos
    outd1 = open("rot_sol_dos.txt", "w")
    outd1.truncate(0)
    for i in range(len(solid_rot_dos)): 
        print(np.float(ff[i]*1e-12),np.float(solid_rot_dos[i]*1e12), file = outd1) 
    outd1.close()  
    outd1 = open("rot_gas_dos.txt", "w")
    outd1.truncate(0)
    for i in range(len(gas_rot_dos)): 
        print(np.float(ff[i]*1e-12),np.float(gas_rot_dos[i]*1e12), file = outd1) 
    outd1.close()  
    W_Shs = hs_rot_entropy_wf(tot_mass, den, rot_fl, rot_dos[0]) # hs approx. rotational entropy weighing function 
    W_Ahs = 0.5 - W_Shs # hs approx. rotational helmholtz energy weighing function
    W_Ehs = 0.5 # hs approx. rotational internal energy weighing function
    rot_s_cho, rot_s_qho, rot_s_gas, rot_s_cho1, rot_s_qho1 = rot_entropy(W_Shs, ff, rot_dos, solid_rot_dos, gas_rot_dos) # rotational part entropy 
    rot_e_cho, rot_e_qho, rot_e_gas, rot_e_cho1, rot_e_qho1 = rot_energy(W_Ehs, ff, rot_dos, solid_rot_dos, gas_rot_dos) # rotational part internal energy
    rot_h_cho, rot_h_qho, rot_h_gas, rot_h_cho1, rot_h_qho1 = rot_helmholtz(W_Ahs, ff, rot_dos, solid_rot_dos, gas_rot_dos) # rotational part helmholtz energy
    print(rot_s_cho)
    
    # Translational system data
    tr_e = md_energy - 3*(1-0.5*rot_fl)*(beta**(-1))
    tr_dif_sp = np.float(tr_dos[0]*(T*kb)/(12*tot_mass*1)) # Translational Self diffusion coefficient (p)
    tr_s_cho = s - rot_s_cho # translational part entropy
    tr_s_qho = s - rot_s_qho # translational part entropy
    tr_h_cho = u - rot_h_cho # translational part helmholtz energy
    tr_h_qho = u - rot_h_qho # translational part helmholtz energy
    a, b, c, d, e = trans_iag_helmholtz(vol, tot_mass, T)
 
    time_stop = perf_counter() # stopping time counter
    run_time = time_stop - time_start # run_time

    out = open(outfile2,"w")
    out.truncate(0)
    print(tr_dif_sp, file = out)
    print(tr_e, file = out)
    print(tr_s_cho, file = out)
    print(u_ex + u_ig_rot - rot_h_cho, file = out)
    print(a, file = out)
    print(b, file = out)
    print(c, file = out)
    print(d, file = out)
    print(e, file = out)
    print(tot_mass, file = out)
    print(tr_s_qho, file = out)
    print(u_ex + u_ig_rot - rot_h_qho, file = out)
    out.close()

    print(run_time)