# python ../../dos_H2O.py ../../2PT_H2O_files/350/RUN2/analysis.tpr ../../2PT_H2O_files/350/RUN2/prod1.trr ../../2PT_H2O_files/350/RUN2/prod1.xtc RUN2-output.txt "-760565.0" "RUN2" ../ref_data/avg_ref_tr_dos.txt ../ref_data/avg_ref_data.txt 

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

# Defining Constants
kb = scipy.constants.k # Boltzmann Constant k
h = scipy.constants.Planck # Planck's constant
Na = scipy.constants.Avogadro # Avogadro's number
T = 300 # Temperature, will make sure to make this as user input
beta = 1/(kb*T) # beta
tot_mass = 18.01528024673462*(1e-3)/Na

def rot_hs_fluidicity(m, den, y0, T, k):
    delta = math.sqrt((np.pi*kb*T)/m)*((6/np.pi)**(2/3))*(den**(1/3))*(2*y0)/(9*k) # normalized diffusivity constant
    def fluidicity(f):
        return (2*(delta**(0-9/2))*(f**(15/2))-(delta**(0-3/2))*(f**(7/2))-6*(f**(5))*(delta**(0-3))+2*f+6*(delta**(0-3/2))*(f**(5/2))-2)
    f_value = optimize.root_scalar(fluidicity,bracket = [0,1],method = 'brentq')
    frac = f_value.root
    return frac

# rotational dos
def rotational_dos(fl, tdos, y, ff):
    p = ((np.pi*y[0]*ff)/(2*fl*tdos)).reshape(-1,1)
    yg = (y[0]/(1+(p**2))).reshape(-1,1) # calculating gas DOS function vs frequency
    ys = y-yg # calculating solid DOS function vs frequency
    return yg, ys

# hard sphere rotational entropy
def hs_rot_entropy_wf():
    r = (((((np.pi)*(np.e)**3)*(T**3))/(13.706211107457026*20.998505222742857*39.467682045442515))**(0.5))*0.5
    si = kb*np.log(r)  # entropy value per molecule basis
    return (si/(3*kb))

# entropy
def entropy(W_Shs, ff, y, ys, yg):
    ygm = W_Shs*yg
    W_Scho = (1 - np.log(beta*h*ff[1:len(ff)]))
    outd1 = open("cho_ent_wf.txt", "w")
    outd1.truncate(0)
    for i in range(len(W_Scho)): 
        print(np.float(ff[i+1]*1e-12),np.float(W_Scho[i]), file = outd1) 
    outd1.close() 
    ysm1 = W_Scho*ys[1:len(ys)]
    outd1 = open("cho_ent_wf_times_sol_dos_for_"+str(W_Shs)+"_W_Shs.txt", "w")
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
    outd1 = open("qho_ent_wf_times_sol_dos_for_"+str(W_Shs)+"_W_Shs.txt", "w")
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

# internal energy 
def energy(W_Ehs, ff, y, ys, yg):
    b1 = beta*h*ff[1:len(ff)]/2
    b2 = (beta*h*ff[1:len(ff)])/(np.exp(beta*h*ff[1:len(ff)])-1)
    W_Eqho = b1+b2
    outd1 = open("qho_ener_wf.txt", "w")
    outd1.truncate(0)
    for i in range(len(W_Eqho)): 
        print(np.float(ff[i+1]*1e-12),np.float(W_Eqho[i]), file = outd1) 
    outd1.close()
    ysm2 = W_Eqho*ys[1:len(ff)]
    outd1 = open("qho_energy_wf_times_sol_dos_for_"+str(W_Ehs)+"_W_Ehs.txt", "w")
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


# helmholtz energy 
def helmholtz(W_Ahs, ff, y, ys, yg):
    W_Acho = np.log(beta*h*ff[1:len(ff)])
    outd1 = open("cho_helm_wf.txt", "w")
    outd1.truncate(0)
    for i in range(len(W_Acho)): 
        print(np.float(ff[i+1]*1e-12),np.float(W_Acho[i]), file = outd1) 
    outd1.close()  
    ysm1 = W_Acho*ys[1:len(ff)]
    outd1 = open("cho_helm_wf_times_sol_dos_for_"+str(W_Ahs)+"_W_Ahs.txt", "w")
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
    outd1 = open("qho_helm_wf_times_sol_dos_for_"+str(W_Ahs)+"_W_Ahs.txt", "w")
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

#--------------------------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    time_start = perf_counter() # starting time counter

    outfile = sys.argv[1]
    enerfile = sys.argv[2]
    #runnum = sys.argv[3]
    Num = sys.argv[3]
    reffile = sys.argv[4] # system argument to specify translational gas reference dos file
    refpropfile = sys.argv[5] # system argument to specify translational reference state properties file

    print("files read")
    """    
    tot_dos = np.array(np.loadtxt(fname = "tot_dos_for_5fs_dt_for_"+runnum+".txt", usecols=1)) # total density of states
    tot_dos = tot_dos.reshape(-1,1)
    tr_dos = np.array(np.loadtxt(fname = "tr_dos_for_5fs_dt_for_"+runnum+".txt", usecols=1)) # total density of states
    tr_dos = tr_dos.reshape(-1,1)
    rot_dos = np.array(np.loadtxt(fname = "rot_dos_for_5fs_dt_for_"+runnum+".txt", usecols=1)) # total density of states
    rot_dos = rot_dos.reshape(-1,1)
    """
    tot_dos = np.array(np.loadtxt(fname = "avg_tot_dos_for_5fs_dt.txt", usecols=1)) # total density of states
    tot_dos = tot_dos.reshape(-1,1)
    tr_dos = np.array(np.loadtxt(fname = "avg_tr_dos_for_5fs_dt.txt", usecols=1)) # total density of states
    tr_dos = tr_dos.reshape(-1,1)
    rot_dos = np.array(np.loadtxt(fname = "avg_rot_dos_for_5fs_dt.txt", usecols=1)) # total density of states
    rot_dos = rot_dos.reshape(-1,1)
    ff = np.array(np.loadtxt(fname = "avg_rot_dos_for_5fs_dt.txt", usecols=0)) # total density of states
    ff = ff.reshape(-1,1)
    tr_dos = tr_dos*(1e-12)
    rot_dos = rot_dos*(1e-12)
    tot_dos = tot_dos*(1e-12)
    ff = ff*(1e12)
    tot_dos_val = integrate.simps(tot_dos.T,ff[:len(tot_dos)].T)
    tr_dos_val = integrate.simps(tr_dos.T,ff[:len(tr_dos)].T)
    rot_dos_val = integrate.simps(rot_dos.T,ff[:len(rot_dos)].T)
    print(np.float(tot_dos_val))
    print(np.float(tr_dos_val))
    print(np.float(rot_dos_val)) 
      
    # system information 
    md_energy = np.float(enerfile)
    print(outfile)
    print(md_energy)
    vol = 4.00440*3.97980*4.01490*(1e-27) # volume of the system
    N = int(np.float(Num))
    den = N/vol # calculating density

    # Translational gas reference system data
    tr_dos_ref = np.loadtxt(fname = reffile)
    dataref = np.loadtxt(fname = refpropfile) 
    tr_dif_sp_ref = np.float(dataref[0]) # value of self-diffusion coefficient for reference state
    tr_e_gas_ref = dataref[1] # value of reference gas md energy
    tr_s_gas_ref_cho = dataref[2] # value of reference gas entropy cho
    tr_s_gas_ref_qho = dataref[9] # value of reference gas entropy cho
    a = dataref[4] # value of curve parameter for Ideal Helmholtz energy
    b = dataref[5] # value of curve parameter for Ideal Helmholtz energy
    c = dataref[6] # value of curve parameter for Ideal Helmholtz energy
    d = dataref[7] # value of curve parameter for Ideal Helmholtz energy
    e = dataref[8] # value of curve parameter for Ideal Helmholtz energy

# Total system data
    md_energy = np.float(enerfile)*1e3/Na + 750269*1e3/Na # total md energy + 750269*1e3/Na
    tot_u_ex_ref_cho = dataref[3] # value of reference excess gibbs free energy
    tot_u_ex_cho = N*tot_u_ex_ref_cho # excess gibbs free energy
    tot_u_ex_ref_qho = dataref[10] # value of reference excess gibbs free energy
    tot_u_ex_qho = N*tot_u_ex_ref_qho # excess gibbs free energy
    
    # Rotational system data
    rot_dif_sp = np.float(rot_dos[0]*(T*kb)/(12*tot_mass*N)) # Rotational Self diffusion coefficient (p)
    rot_fl = rot_hs_fluidicity(tot_mass, den, rot_dos[0], T, N) # Rotational fluidicity
    print(rot_fl)
    gas_rot_dos, solid_rot_dos = rotational_dos(rot_fl, rot_dos_val, rot_dos, ff) # gas and solid rotational dos
    """    
    outd1 = open("solid_rot_dos_for_5fs_dt_for_"+runnum+".txt", "w") 
    outd2 = open("gas_rot_dos_for_5fs_dt_for_"+runnum+".txt", "w")
    """   
    outd1 = open("avg_solid_rot_dos_for_5fs_dt.txt", "w") 
    outd2 = open("avg_gas_rot_dos_for_5fs_dt.txt", "w")
    outd1.truncate(0) 
    outd2.truncate(0)
    for i in range(len(solid_rot_dos)): 
        print(np.float(ff[i]*1e-12),np.float(solid_rot_dos[i]*1e12), file = outd1) 
    outd1.close() 
    for i in range(len(gas_rot_dos)): 
        print(np.float(ff[i]*1e-12),np.float(gas_rot_dos[i]*1e12), file = outd2) 
    outd2.close() 
    W_Shs = hs_rot_entropy_wf() # hs approx. rotational entropy weighing function 
    W_Ahs = 0.5 - W_Shs # hs approx. rotational helmholtz energy weighing function
    W_Ehs = 0.5 # hs approx. rotational internal energy weighing function
    rot_s_cho, rot_s_qho, rot_s_gas, rot_s_cho1, rot_s_qho1 = entropy(W_Shs, ff, rot_dos, solid_rot_dos, gas_rot_dos) # rotational part entropy 
    rot_e_cho, rot_e_qho, rot_e_gas, rot_e_cho1, rot_e_qho1  = energy(W_Ehs, ff, rot_dos, solid_rot_dos, gas_rot_dos) # rotational part internal energy
    rot_h_cho, rot_h_qho, rot_h_gas, rot_h_cho1, rot_h_qho1  = helmholtz(W_Ahs, ff, rot_dos, solid_rot_dos, gas_rot_dos) # rotational part helmholtz energy
    rot_ln_q_cho = -rot_h_cho*beta # ln of rotational partition function
    rot_ln_q_qho = -rot_h_qho*beta # ln of rotational partition function
    print(rot_s_cho)
    
    # Translational system data
    tr_dif_sp = np.float(tr_dos[0]*(T*kb)/(12*tot_mass*N)) # Translational Self diffusion coefficient (p)
    tr_fl = tr_dif_sp/tr_dif_sp_ref # fluidicity
    gas_tr_dos = tr_fl*tr_dos_ref.reshape(-1,1)*N # calculating gas DOS function vs frequency
    solid_tr_dos = tr_dos - gas_tr_dos # calculating solid DOS function vs frequency
    """    
    outd1 = open("solid_tr_dos_for_5fs_dt_for_"+runnum+".txt", "w") 
    outd2 = open("gas_tr_dos_for_5fs_dt_for_"+runnum+".txt", "w")
    """
    outd1 = open("avg_solid_tr_dos_for_5fs_dt.txt", "w") 
    outd2 = open("avg_gas_tr_dos_for_5fs_dt.txt", "w")   
    outd1.truncate(0) 
    outd2.truncate(0)
    for i in range(len(solid_tr_dos)): 
        print(np.float(ff[i]*1e-12),np.float(solid_tr_dos[i]*1e12), file = outd1) 
    outd1.close() 
    for i in range(len(gas_tr_dos)): 
        print(np.float(ff[i]*1e-12),np.float(gas_tr_dos[i]*1e12), file = outd2) 
    outd2.close() 
    tr_u_ex_cho = tot_u_ex_cho # - rot_e_cho # excess translational gibbs energy cho
    tr_u_ex_qho = tot_u_ex_qho # - rot_e_cho # excess translational gibbs energy qho
    tr_h_gas_cho = tr_fl*tr_u_ex_cho + (a + b*(N*tr_fl) + c*((N*tr_fl)**2) + d*((N*tr_fl)**3) + e*((N*tr_fl)**4)) 
    tr_h_gas_qho = tr_fl*tr_u_ex_qho + (a + b*(N*tr_fl) + c*((N*tr_fl)**2) + d*((N*tr_fl)**3) + e*((N*tr_fl)**4)) 
    tr_h_solid_cho, tr_h_solid_qho, tr_h_gas, tr_h_solid_cho1, tr_h_solid_qho1= helmholtz(0, ff, tr_dos, solid_tr_dos, gas_tr_dos) # translational solid part entropy 
    tr_h_cho = tr_h_gas_cho + tr_h_solid_cho
    tr_h_qho = tr_h_gas_qho + tr_h_solid_qho
    tr_e_solid_cho, tr_e_solid_qho, tr_e_gas, tr_e_solid_cho1, tr_e_solid_qho1 = energy(0, ff, tr_dos, solid_tr_dos, gas_tr_dos) # translational solid part entropy 
    tr_e_cho = N*tr_fl*tr_e_gas_ref + tr_e_solid_cho # translational part internal energy
    tr_e_qho = N*tr_fl*tr_e_gas_ref + tr_e_solid_qho
    tr_s_solid_cho, tr_s_solid_qho, tr_s_gas, tr_s_solid_cho1, tr_s_solid_qho1 = entropy(0, ff, tr_dos, solid_tr_dos, gas_tr_dos) # translational solid part energy
    tr_s_cho = (N*tr_fl*tr_e_gas_ref - tr_h_gas_cho)/T + tr_s_solid_cho # translational part entropy
    tr_s_qho = (N*tr_fl*tr_e_gas_ref - tr_h_gas_qho)/T + tr_s_solid_qho
    tr_ln_q_cho = -tr_h_cho*beta # ln of translational partition function
    tr_ln_q_qho = -tr_h_qho*beta # ln of translational partition function
 
    # Total system data
    reference = md_energy - ((3*N)*(2 - 0.5*rot_fl - tr_fl))/(beta) - N*tr_fl*tr_e_gas_ref # reference energy
    tot_e_cho = reference + tr_e_cho + rot_e_cho # total internal energy
    tot_s_cho = tr_s_cho + rot_s_cho # total entropy
    tot_h_cho = (md_energy - tot_s_cho*T) # total helmholtz free energy
    tot_ln_q_cho = -tot_h_cho*beta # ln of total partition function
    tot_e_qho = reference + tr_e_qho + rot_e_qho # total internal energy
    tot_s_qho = tr_s_qho + rot_s_qho # total entropy
    tot_h_qho = (md_energy - tot_s_qho*T) # total helmholtz free energy
    tot_ln_q_qho = -tot_h_qho*beta # ln of total partition function
    
    time_stop = perf_counter() # stopping time counter
    run_time = time_stop - time_start # run_time
    print(run_time)

    out = open(outfile,"w")
    out.truncate(0)
    print("Data collected for N = " + str(N) + " molecules system", file = out)
    print("All values are in their SI units \n", file = out)
    print("\n Total system data :\n", file = out)
    print("Total Density of states in the system = ", np.float(tot_dos_val), file = out)
    print("MD energy of the system scales with N (normalized by beta) = ", md_energy/(kb*T), file = out)
    print("Total Helmholtz free energy of the system (assuming classical harmonic oscillator for solid part of the system) scales with N (normalized by beta) = ", tot_h_cho/(kb*T), file = out)
    print("Total Helmholtz free energy of the system (assuming quantum harmonic oscillator for solid part of the system) scales with N (normalized by beta) = ", tot_h_qho/(kb*T), file = out)
    print("Total Internal energy of the system (assuming classical harmonic oscillator for solid part of the system) scales with N (normalized by beta) = ", tot_e_cho/(kb*T), file = out)
    print("Total Internal energy of the system (assuming quantum harmonic oscillator for solid part of the system) scales with N (normalized by beta) = ", tot_e_qho/(kb*T), file = out)
    print("Total Entropy (assuming classical harmonic oscillator for solid part of the system) (normalizd by kb) scales with N = ", tot_s_cho/kb, file = out)
    print("Total Entropy (assuming quantum harmonic oscillator for solid part of the system) (normalizd by kb) scales with N = ", tot_s_qho/kb, file = out)
    print("Reference potential energy of the system scales with N (normalized by beta) = ", reference/(kb*T), file = out)
    print("ln of Total Partition function of the system (assuming classical harmonic oscillator for solid part of the system) scales with N = ", tot_ln_q_cho, file = out)
    print("ln of Total Partition function of the system (assuming quantum harmonic oscillator for solid part of the system) scales with N = ", tot_ln_q_qho, file = out)
    print("\n\n Rotational System Data: \n", file = out)
    print("Total Rotational Density of states in the system = ", np.float(rot_dos_val), file = out)
    print("Rotational Self Diffusion Coefficient of the system = ", rot_dif_sp, file = out)
    print("Rotational Fluidicity / Fraction of gas = ", rot_fl, file = out)
    print("Total number of gas Rotational DOS in the system = ", rot_fl*np.float(rot_dos_val), file = out)
    print("Total number of solid Rotational DOS in the system = ", (1-rot_fl)*np.float(rot_dos_val), file = out)
    print("Total Rotational Helmholtz free energy of the system (assuming classical harmonic oscillator for solid part of the system) scales with N (normalized by beta) = ", rot_h_cho/(kb*T), file = out)
    print("Total Rotational Helmholtz free energy of the system (assuming quantum harmonic oscillator for solid part of the system) scales with N (normalized by beta) = ", rot_h_qho/(kb*T), file = out)
    print("Gas part of Total Rotational Helmholtz free energy of the system (assuming hard sphere approximation for gas part) scales with N (normalized by beta) = ", rot_h_gas/(kb*T), file = out)    
    print("Solid part of Total Rotational Helmholtz free energy of the system (assuming classical harmonic oscillator for solid part of the system) scales with N (normalized by beta) = ", rot_h_cho1/(kb*T), file = out)
    print("Solid part of Total Rotational Helmholtz free energy of the system (assuming quantum harmonic oscillator for solid part of the system) scales with N (normalized by beta) = ", rot_h_qho1/(kb*T), file = out)    
    print("Total Rotational energy of the system (assuming classical harmonic oscillator for solid part of the system) scales with N (normalized by beta) = ", rot_e_cho/(kb*T), file = out)
    print("Total Rotational energy of the system (assuming quantum harmonic oscillator for solid part of the system) scales with N (normalized by beta) = ", rot_e_qho/(kb*T), file = out)
    print("Gas part of Total Rotational energy of the system (assuming hard sphere approximation for gas part) scales with N (normalized by beta) = ", rot_e_gas/(kb*T), file = out)    
    print("Solid part of Total Rotational energy of the system (assuming classical harmonic oscillator for solid part of the system) scales with N (normalized by beta) = ", rot_e_cho1/(kb*T), file = out)
    print("Solid part of Total Rotational energy of the system (assuming quantum harmonic oscillator for solid part of the system) scales with N (normalized by beta) = ", rot_e_qho1/(kb*T), file = out)    
    print("Total Rotational Entropy (assuming classical harmonic oscillator for solid part of the system) (normalized by kb) scales with N = ", rot_s_cho/kb, file = out)
    print("Total Rotational Entropy (assuming quantum harmonic oscillator for solid part of the system) (normalized by kb) scales with N = ", rot_s_qho/kb, file = out)
    print("Gas part of Total Rotational Entropy of the system (assuming hard sphere approximation for gas part) scales with N (normalized by beta) = ", rot_s_gas/(kb), file = out)    
    print("Solid part of Total Rotational Entropy of the system (assuming classical harmonic oscillator for solid part of the system) scales with N (normalized by beta) = ", rot_s_cho1/(kb), file = out)
    print("Solid part of Total Rotational Entropy of the system (assuming quantum harmonic oscillator for solid part of the system) scales with N (normalized by beta) = ", rot_s_qho1/(kb), file = out)
    print("ln of Rotational Partition function of the system (assuming classical harmonic oscillator for solid part of the system) scales with N = ", rot_ln_q_cho, file = out)
    print("ln of Rotational Partition function of the system (assuming quantum harmonic oscillator for solid part of the system) scales with N = ", rot_ln_q_qho, file = out)
    print("\n\n Translational System Data: \n", file = out)
    print("Total Translational Density of states in the system = ", np.float(tr_dos_val), file = out)
    print("Translational Self Diffusion Coefficient of the system = ", tr_dif_sp, file = out)
    print("Translational Fluidicity / Fraction of gas = ", tr_fl, file = out)
    print("Total number of gas Translational DOS in the system = ", tr_fl*np.float(tr_dos_val), file = out)
    print("Total number of solid Translational DOS in the system = ", (1-tr_fl)*np.float(tr_dos_val), file = out)
    print("Total Translational Helmholtz free energy of the system (assuming classical harmonic oscillator for solid part of the system) scales with N (normalized by beta) = ", tr_h_cho/(kb*T), file = out)
    print("Total Translational Helmholtz free energy of the system (assuming quantum harmonic oscillator for solid part of the system) scales with N (normalized by beta) = ", tr_h_qho/(kb*T), file = out)
    print("Gas part of Total Translational Helmholtz energy of the system (assuming IAG approximation for gas part and cho for reference) scales with N (normalized by beta) = ", tr_h_gas_cho/(kb*T), file = out) 
    print("Gas part of Total Translational Helmholtz energy of the system (assuming IAG approximation for gas part and qho for reference) scales with N (normalized by beta) = ", tr_h_gas_qho/(kb*T), file = out)   
    print("Solid part of Total Translational Helmholtz energy of the system (assuming classical harmonic oscillator for solid part of the system) scales with N (normalized by beta) = ", tr_h_solid_cho/(kb*T), file = out)
    print("Solid part of Total Translational Helmholtz energy of the system (assuming quantum harmonic oscillator for solid part of the system) scales with N (normalized by beta) = ", tr_h_solid_qho/(kb*T), file = out)
    print("Total Translational energy of the system (assuming classical harmonic oscillator for solid part of the system) scales with N (normalized by beta) = ", tr_e_cho/(kb*T), file = out)
    print("Total Translational energy of the system (assuming quantum harmonic oscillator for solid part of the system) scales with N (normalized by beta) = ", tr_e_qho/(kb*T), file = out)
    print("Gas part of Total Translational energy of the system (assuming IAG approximation for gas part) scales with N (normalized by beta) = ", (N*tr_fl*tr_e_gas_ref)/(kb*T), file = out)   
    print("Solid part of Total Translational energy of the system (assuming classical harmonic oscillator for solid part of the system) scales with N (normalized by beta) = ", tr_e_solid_cho/(kb*T), file = out)
    print("Solid part of Total Translational energy of the system (assuming quantum harmonic oscillator for solid part of the system) scales with N (normalized by beta) = ", tr_e_solid_qho/(kb*T), file = out)
    print("Total Translational Entropy (assuming classical harmonic oscillator for solid part of the system) (normalized by kb) scales wih N = ", tr_s_cho/kb, file = out)
    print("Total Translational Entropy (assuming quantum harmonic oscillator for solid part of the system) (normalized by kb) scales wih N = ", tr_s_qho/kb, file = out)
    print("Gas part of Total Translational Entropy of the system (assuming IAG approximation for gas part and cho for the solid part) scales with N (normalized by beta) = ", (N*tr_fl*tr_e_gas_ref - tr_h_gas_cho)/(kb*T), file = out) 
    print("Gas part of Total Translational Entropy of the system (assuming IAG approximation for gas part and qho for the solid part) scales with N (normalized by beta) = ", (N*tr_fl*tr_e_gas_ref - tr_h_gas_qho)/(kb*T), file = out)  
    print("Solid part of Total Translational Entropy of the system (assuming classical harmonic oscillator for solid part of the system) scales with N (normalized by beta) = ", tr_s_solid_cho/(kb), file = out)
    print("Solid part of Total Translational Entropy of the system (assuming quantum harmonic oscillator for solid part of the system) scales with N (normalized by beta) = ", tr_s_solid_qho/(kb), file = out)
    print("ln of Translational Partition function of the system (assuming classical harmonic oscillator for solid part of the system) scales with N = ", tr_ln_q_cho, file = out)
    print("ln of Translational Partition function of the system (assuming quantum harmonic oscillator for solid part of the system) scales with N = ", tr_ln_q_qho, file = out)
    out.close()

"""
    plt.figure(1)
    plt.plot(ff[:len(ff)-1], tot_dos, "k--", label="total DoS")
    plt.plot(ff[:len(ff)-1], tr_dos, "r--", label="translational DoS")
    plt.plot(ff[:len(ff)-1], rot_dos, "b--", label="rotational DoS")
    plt.xlabel("freq in s-1")
    plt.ylabel("S(v) in s")
    plt.legend()
    plt.show()

    plt.figure(2)
    plt.plot(ff[:len(ff)-1], rot_dos, "k--", label="total rotational DoS") 
    plt.plot(ff[:len(ff)-1], solid_rot_dos, "r--", label="rotational solid DoS") 
    plt.plot(ff[:len(ff)-1], gas_rot_dos, "b--", label="rotational gas DoS") 
    plt.xlabel("freq in s-1") 
    plt.ylabel("S(v) in s") 
    plt.legend() 
    plt.show()  

    plt.figure(3)
    plt.plot(ff[:len(ff)-1], tr_dos, "k--", label="total translational DoS") 
    plt.plot(ff[:len(ff)-1], solid_tr_dos, "r--", label="translational solid DoS") 
    plt.plot(ff[:len(ff)-1], gas_tr_dos, "b--", label="translational gas DoS") 
    plt.xlabel("freq in s-1") 
    plt.ylabel("S(v) in s") 
    plt.legend() 
    plt.show()  
"""



