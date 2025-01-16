import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
import sys
import math
import scipy.constants
from scipy.fftpack import fft, ifft, fftfreq
import MDAnalysis as mda
import sys
from scipy import integrate
from scipy import optimize
from time import perf_counter 
from scipy.optimize import curve_fit
import psutil

# Defining Consatants
kb = scipy.constants.k # Boltzmann Constant k
h = scipy.constants.Planck # Planck's constant
Na = scipy.constants.Avogadro # Avogadro's number
T = 298.15 # Temperature
beta = 1/(kb*T) # beta

# extracting atom weight and storing them in SI units (kg/particle)
def extract_atom_weight(b):
    mass_dict = dict({i+1: b[i].mass*(1e-3)/Na for i in range(len(b))})
    mass = np.array([x for x in mass_dict.values()])
    return mass.reshape(-1,1) # reshaping them in 4x1 array

# extacting velocity and time and storing them in SI units (ms-1, s)
def extract_velocity_and_time(univ, x):
    print(psutil.virtual_memory())
    vel = np.zeros((len(univ),len(x),3))
    time = np.zeros((len(vel),1))
    print("after preallocation of vel and time arrays")
    print(psutil.virtual_memory())
    for i in range(len(univ)):
        univ[i]
        vel[i] = x.velocities
        time[i] = univ[i].time
    print("after extracting vel and time data")
    print(psutil.virtual_memory())
    return vel*1e2, (time.reshape(-1))*(1e-12)
    return c[:n]

# extracting positions (m), c.o.m (m), MI (kg*m^2), MI inverse, and principal moments of inertia 
def extract_pos_com_MI_prinaxes(univ, x, y, red):
    print(psutil.virtual_memory())
    l = np.zeros((len(univ),len(red)*3,3)) # list to store moment of inertia of all molecules at all time steps)
    pos = np.zeros((len(univ),len(red.atoms),3)) # array to store positions of all atoms 
    print("after preallocation of pos and MI arrays")
    print(psutil.virtual_memory())
    for i in range(len(univ)): 
        univ[i]
        q = np.zeros((len(red),3,3)) # array to store moment of inertia of all molecules at one time step
        for j in range(len(red)):
            q[j] = red[j].atoms.moment_of_inertia() # moment of inertia
        l[i] = q.reshape(-1,3)
        pos[i] = red.atoms.positions
    print("after extracting pos and MI data")
    print(psutil.virtual_memory())
    return pos*(1e-10), l*(1e-20)/(1e3*Na)

# calculate c.o.m. information from extracted positions and velocities and store them in SI units
def com_data(pos, vel, mass):
    tot_mass = np.sum(mass) # total mass of a molecule
    m = np.tile(mass.transpose(), int(len(vel[0])/len(mass))).transpose()
    m_v = np.zeros((len(vel), len(vel[0]),3))
    m_p = np.zeros((len(pos), len(pos[0]),3))
    for i in range(len(vel)):
        m_v[i] = vel[i]*m/tot_mass
        m_p[i] = pos[i]*m/tot_mass
    #m_v = vel*m/tot_mass
    #m_p = pos*m/tot_mass
    print("after m_v and m_p calculation for c.o.m")
    print(psutil.virtual_memory())
    return m_p.reshape(-1,len(mass),3).sum(1).reshape(len(pos),-1,3), m_v.reshape(-1,len(mass),3).sum(1).reshape(len(vel),-1,3)

# calculate principal moments of inertia from Moment of Inertia Tensor (MIT)
def prinaxes(MI,x):
    print(psutil.virtual_memory())
    v,w = np.linalg.eig(MI.reshape(-1,3,3)) # calculating eigen values and eigen vectors for the MIT inputs
    a1 = np.argsort(-v) # sorting eigen values in desending order and storing the indices from largest to smallest eigen values for each MIT
    a3 = np.zeros((len(w),3,3)) # list to store eigen vectors (column wise) corresponding to each eigen values
    a4 = np.zeros((len(v),3,3)) # list to store eigen values for each MIT in descending order
    a4inv = np.zeros((len(v),3,3))
    for i in range(len(a1)):
        a3[i] = (w[i][:,a1[i]]) # storing eigen vectors corresponding to decreasing eigen values
        a4[i] = np.diag((v[i][a1[i]])) # storing Principal Moments of Inertia (PMIT) which is the diagonal matrix of eigen values
        a4inv[i] = np.linalg.inv(a4[i])
    a3 = a3.reshape(x,-1,3)
    a4 = a4.reshape(x,-1,3)
    a4inv = a4inv.reshape(x,-1,3)
    #a4inv = np.linalg.inv(a4.reshape(-1,3,3)).reshape(x,-1,3) # calculating inverse of PMIT
    v = []
    w = []
    print("after computing MI data using numpy")
    print(psutil.virtual_memory())
    return a4, a4inv, a3

# function to allign and orient principal axes consistently 
def mod_prinaxes(MI):
    ctr = np.zeros((int(len(MI[0])/3),3))
    ctc = 0
    for i in range(0,len(MI[0]),3):
        for j in range(1,len(MI)):
            a = np.zeros((3,1))
            a[0] = np.dot(MI[j-1,i:i+3,0],MI[j,i:i+3,0])
            a[1] = np.dot(MI[j-1,i:i+3,1],MI[j,i:i+3,1])
            a[2] = np.dot(MI[j-1,i:i+3,2],MI[j,i:i+3,2])
            b = np.abs(a)
            c = np.argsort(-b.T).T
            if (j==1):
                ctr[ctc] = c.T
                ctc = ctc+1
            if (a[int(c[0])]*b[int(c[0])] < 0):
                MI[j,i:i+3,int(c[0])] = -MI[j,i:i+3,int(c[0])]
            if (a[int(c[1])]*b[int(c[1])] < 0):
                MI[j,i:i+3,int(c[1])] = -MI[j,i:i+3,int(c[1])]
            if (int(np.abs(c[0]-c[1])) == 1):
                MI[j,i:i+3,int(c[2])] = np.cross(MI[j,i:i+3,min(int(c[0]),int(c[1]))], MI[j,i:i+3,max(int(c[0]),int(c[1]))])
            else:
                MI[j,i:i+3,int(c[2])] = np.cross(MI[j,i:i+3,max(int(c[0]),int(c[1]))], MI[j,i:i+3,min(int(c[0]),int(c[1]))])
    for i in range(len(ctr)):
        if(np.linalg.det(MI[0,3*i:3*i+3,:]) < 0):
            MI[0,3*i:3*i+3,int(ctr[i][2])] = -MI[0,3*i:3*i+3,int(ctr[i][2])]
    return MI 

# calculating MIT from conventional formula
def moment_of_inertia(mass, vec):
    m = np.tile(mass.transpose(), int(len(vec[0])/len(mass)))
    xx = m*(vec[:,:,1]**2 + vec[:,:,2]**2) # I_xx
    xxu = xx.reshape(-1,len(mass)).sum(1).reshape(len(vec),-1)
    yy = m*(vec[:,:,0]**2 + vec[:,:,2]**2) # I_yy
    yyu = yy.reshape(-1,len(mass)).sum(1).reshape(len(vec),-1)
    zz = m*(vec[:,:,0]**2 + vec[:,:,1]**2) # I_zz
    zzu = zz.reshape(-1,len(mass)).sum(1).reshape(len(vec),-1)
    xy = -m*(vec[:,:,0]*vec[:,:,1]) # I_xy = I_yx
    xyu = xy.reshape(-1,len(mass)).sum(1).reshape(len(vec),-1)
    xz = -m*(vec[:,:,0]*vec[:,:,2]) # I_xz = I_zx
    xzu = xz.reshape(-1,len(mass)).sum(1).reshape(len(vec),-1)
    yz = -m*(vec[:,:,1]*vec[:,:,2])# I_yz = I_zy
    yzu = yz.reshape(-1,len(mass)).sum(1).reshape(len(vec),-1)
    MIl = [] # list to store MIT
    for i in range(len(xxu)):
        MIlts = []
        for j in range(len(xxu[0])):
            k = np.array([[xxu[i,j],xyu[i,j],xzu[i,j]],[xyu[i,j],yyu[i,j],yzu[i,j]],[xzu[i,j],yzu[i,j],zzu[i,j]]])
            MIlts.append(k)
        MIats = np.array(MIlts).reshape(len(xxu[0])*3,3)
        MIl.append(MIats)
    MIa = np.array(MIl)
    MIainv = np.linalg.inv(MIa.reshape(-1,3,3)).reshape(len(vec),-1,3) # calculating inverse of MIT 
    return MIa, MIainv

# calculating data along the principal axes of rotation (Xn = X*P', where Xn -> transformed X and P' -> the transformtion matrix with each columns being the principal axes of rotation in our case)
def data_along_prinaxes(pos, prinaxes):
    pos = pos.reshape(-1,4,3)
    prinaxes = prinaxes.reshape(-1,3,3)
    pr_pos = np.zeros((len(pos),4,3))
    print("after preallocation of pr_axis transformed rel_pos")
    print(psutil.virtual_memory())
    for i in range(len(pos)):
        pr_pos[i] = np.matmul(pos[i],prinaxes[i])
    return pr_pos.reshape(len(pos),-1,3)

# calculating relative distance and velocities for each atom of a molecule w.r.t. its c.o.m. parameters
def rel_vector(comp, pos):
    rel_dis = pos - comp # relative position
    return rel_dis

# calculating angular momentum (m*(r x v))
def angular_momentum(mass, vec, vel):
    m = np.tile(mass.transpose(), int(len(vel[0])/len(mass))).transpose()
    m_v = np.zeros((len(vel),len(vel[0]),3))
    for i in range(len(vel)):
        k = vel[i]*m
        m_v[i] = k
    #m_v = vel*m
    print("after calculating mass*velocity array")
    print(psutil.virtual_memory())
    print("storage for m_v in gb = %f"%(m_v.nbytes/(1024*1024*1024)))
    L_i = np.zeros((len(vel),len(vel[0]),3))
    print("after preallocation of L(ang. mom.) array")
    print(psutil.virtual_memory())
    for i in range(len(vel)):
        L_i[i] = np.cross(vec[i], m_v[i])
    return L_i.reshape(-1,len(mass),3).sum(1).reshape(len(vel),-1,3) 

# calculating angular velocity from angular momentum (MITinv * L)
def rot_angular_velocity(MIinv, L, k):
    w = np.matmul(MIinv.reshape(-1,3,3),L.reshape(-1,3,1)).reshape(k,-1,1) 
    return w  

# PMIT from literature for SPC:
def PMIT_SPC(mass):
    lpos = np.array([[[0,0,0],[1,0,0],[-0.3333132476,0.9428161427,0]]]) # dummy coordinates set according to data for Tip4p mwater model (r(OH) = 0.9572 angs, HOH angle = 104.52 deg) 
    tot_mass = np.sum(mass) # total mass of a molecule 
    #part to calculate c.o.m.
    m_p = lpos*mass/tot_mass 
    lcom = m_p.reshape(-1,len(mass),3).sum(1).reshape(len(lpos),-1,3) 
    lcom = np.repeat(lcom, repeats=len(mass), axis=1) # c.o.m data
    rpos = (lpos-lcom)*(1e-10) # relative positions of atoms from c.o.m.
    lmit, lmiti = moment_of_inertia(mass,rpos) # MIT tensor
    plmit, plmiti, pla = prinaxes(lmit,len(lpos)) # principal MIT
    return plmit, plmiti

# PMIT from literature for TIP4P:
def PMIT_T4P(mass):
    lpos = np.array([[[0,0,0],[0.9572,0,0],[-0.23999,0.92663,0],[3,4,5]]]) # dummy coordinates set according to data for Tip4p mwater model (r(OH) = 0.9572 angs, HOH angle = 104.52 deg) 
    tot_mass = np.sum(mass) # total mass of a molecule 
    #part to calculate c.o.m.
    m_p = lpos*mass/tot_mass 
    lcom = m_p.reshape(-1,len(mass),3).sum(1).reshape(len(lpos),-1,3) 
    lcom = np.repeat(lcom, repeats=len(mass), axis=1) # c.o.m data
    rpos = (lpos-lcom)*(1e-10) # relative positions of atoms from c.o.m.
    lmit, lmiti = moment_of_inertia(mass,rpos) # MIT tensor
    plmit, plmiti, pla = prinaxes(lmit,len(lpos)) # principal MIT
    return plmit, plmiti


# calculating spectral density
def vel_fft(vel):
    n = len(vel)
    c = vel
    c = np.power(np.absolute(fft(c)), 2)/len(c)
    return c

# calculating mass weighted autocorrelation values
def ini_dos(vel, m):
    corr = [[],[],[]]
    for dim in range(3):
        ma = m[:,dim]
        for atom in range(vel.shape[1]):
            corr[dim].append(vel_fft(vel[:, atom, dim]))
            corr[dim][atom] = corr[dim][atom]*ma[atom]
    corr = np.array(corr, dtype=np.float_)
    return corr

# calculating DoS 
def density_of_states(corr, dt):
    dos = np.sum(corr, axis=(0,1))
    ff = fftfreq(len(dos), d = dt)
    return dos[:len(dos)//2], ff[:len(dos)//2]

def integ_DoS(y, ff, dt):
    y = y*(beta)*2*dt
    integ = integrate.simps(y.T,ff[:len(y)].T)
    return y, np.float(integ)

# all required box frame data
def box_frame_data(pos, vel, mass, k):
    c_pos, c_vel = com_data(pos, vel, mass) # c.o.m. positions and velocity 
    com_pos = np.repeat(c_pos, repeats=k, axis=1)
    com_vel = np.repeat(c_vel, repeats=k, axis=1)
    rel_dis = rel_vector(com_pos, pos)
    c_pos = []
    c_vel = []
    return rel_dis, com_vel 

# all principal axes frame data
def prinaxes_data(rel_dis, vel, np_prinaxes, mass, np_MIinv):
    #rel_dispr = data_along_prinaxes(rel_dis, np_prinaxes) # translating relative positions and velocities along principal axes frame
    #MI, MIinv =  moment_of_inertia(mass, rel_dispr) # MI calculated in pricipal axes frame (should be a diagonal matrix)
    #print("storage for MIinv in gb = %f"%((MIinv.nbytes)/(1024*1024*1024)))
    L = angular_momentum(mass, rel_dis, vel) # calculating angular momnetum 
    print("storage for L in gb = %f"%((L.nbytes)/(1024*1024*1024)))
    Lt = ((np.transpose(np_prinaxes.reshape(-1,3,3), (0,2,1)))@(np.transpose(L.reshape(-1,1,3), (0,2,1)))).reshape(len(vel),-1,3) # calculating angular momentum in principal axes frame
    print("storage for Lt in gb = %f"%((Lt.nbytes)/(1024*1024*1024)))
    L = []
    print("before calculating angular velocity")
    print(psutil.virtual_memory())
    w_np = rot_angular_velocity(np_MIinv, Lt.reshape(len(vel),-1), len(vel)) # calculating angular velocity
    print("before calculating angular velocity")
    print(psutil.virtual_memory())
    aMIT, aMITi = PMIT_SPC(mass)
    aMIT = aMIT.reshape(-1,3)
    w_np = w_np.reshape(len(vel),-1,3)
    return aMIT, w_np  

# dos data
def dos_data(x, y, dt):
    dos_pr, ff1 = density_of_states(ini_dos(x,y), dt) # calculating density of states with data in principal axes frame
    dos_pr = dos_pr.reshape(len(dos_pr),1)
    dos_req, tdos = integ_DoS(dos_pr, ff1, dt) # calculating total density of states    
    return dos_req, tdos, ff1

def auto_corr_fft(vel):
    n = len(vel)
    #print(n)
    p = vel.astype("complex")
    #p = np.conj(p)*complex(0,1)
    norm_fac = np.arange(n, 0, -1, dtype=np.float_)
    c = np.hstack((vel,np.zeros(n)))
    #c = vel
    c = c.astype("complex")
    c = np.real(ifft(np.power(np.absolute(fft(c)), 2)))
    c[:n] = c[:n]/norm_fac[:n]
    #c[n:] = c[n:]/norm_fac[::-1]
    #print(len(c))
    #return c[:n]
    return c[:n]

# calculating mass weighted autocorrelation values
def corr_3D(vel, m):
    corr = [[],[],[]]
    for dim in range(3):
        ma = m[:,dim]
        for atom in range(vel.shape[1]):
            corr[dim].append(auto_corr_fft(vel[:, atom, dim]))
            corr[dim][atom] = corr[dim][atom]*ma[atom]
    corr = np.array(corr, dtype=np.float_)
    return corr

# calculating total autocorrelation 
def total_autocorr(corr):
    total_ac = np.sum(corr, axis=(0,1))
    return total_ac

# calculating DoS from total auto-correlation values
def density_of_states2(ac, dt):
    length = len(ac)
    yfft = np.empty(length, dtype = np.float_)
    yfft = fft(ac) 
    yfft = np.real(yfft)
    ff = fftfreq(len(yfft), d = dt)
    return yfft[:length//2], ff[:length//2]

# calculating Total DoS
def integ_DoS2(y, ff, dt):
    y = y*(beta)*2*dt
    integ = integrate.simps(y.T,ff[:len(y)].T)
    return y, np.float(integ)

def dos_data2(x, y, dt):
    dos_pr, ff2 = density_of_states2(total_autocorr(corr_3D(x,y)), dt) # calculating density of states with data in principal axes frame
    #dos_pr, ff2 = density_of_states2(np.hstack((total_autocorr(corr_3D(x,y)),total_autocorr(corr_3D(x,y))[::-1])), dt)
    dos_pr = dos_pr.reshape(len(dos_pr),1)
    dos_req, tdos = integ_DoS2(dos_pr, ff2, dt) # calculating total density of states    
    return dos_req, tdos, ff2

# all universe data
def extract_data_from_universe(u1, u2, x, y):
    vel, time = extract_velocity_and_time(u1, x) # extracting velocities and time
    print("vel read")
    pos, md_MI = extract_pos_com_MI_prinaxes(u2, len(vel), len(y), y)
    print("pos read")
    print(md_MI.shape)
    time = time.reshape(len(time),1)
    mass = extract_atom_weight(y[0].atoms)
    tot_mass = float(sum(mass))
    print(tot_mass)
    N = len(y) # storing number of molecules  
    dt = np.float((time[1]-time[0])) # length of each time step
    #ff = np.arange(len(time), dtype=np.float_)/(len(time)*dt) # storing frequency
    #ff = ff.reshape(len(time),1)
    return N, md_MI, vel, time, pos, mass, dt, tot_mass

#--------------------------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    time_start = perf_counter() # starting time counter

    analysis_tpr_file = sys.argv[1] # system argument to specify analysis.tpr file
    trr_file = sys.argv[2] # system argument to specify .trr file
    #xtc_file = sys.argv[3] # system argument to specify .xtc file
    runnum = sys.argv[3]

    print("files read")
   
    # system information 
    u1 = mda.Universe(analysis_tpr_file,trr_file) # creating universe with velocities
    print(u1)
    #u2 = mda.Universe(analysis_tpr_file,xtc_file) # creating universe with positions
    #print(u2)
    wat_atoms = u1.select_atoms("resname SOL")
    print(wat_atoms)
    vol = 2.48262*2.48262*2.48262*(1e-27) # volume of the system
    p1 = u1.trajectory[0:len(u1.trajectory):10]
    print(p1)
    #p2 = u2.trajectory[0:len(u2.trajectory):10]
    #print(p2)
    N, md_MI, vel, time, pos, mass, dt, tot_mass = extract_data_from_universe(p1, p1, wat_atoms, u1.residues)
    print(N)
    den = N/vol # calculating density
    print(vel.shape)
    print("storage for vel in gb = %f"%((vel.nbytes)/(1024*1024*1024)))
    print("storage for md_MI in gb = %f"%((md_MI.nbytes)/(1024*1024*1024)))
    np_MI, np_MIinv, np_prinaxes = prinaxes(md_MI,len(vel))
    print(np_prinaxes.shape)
    np_prinaxes_mod = mod_prinaxes(np.copy(np_prinaxes))
    print("storage for np_MI in gb = %f"%((np_MI.nbytes)/(1024*1024*1024)))
    print("storage for np_prinaxes in gb = %f"%((np_prinaxes.nbytes)/(1024*1024*1024)))
    rel_dis, com_vel = box_frame_data(pos, vel, mass, len(u1.residues[0].atoms))
    print("storage for rel_dis in gb = %f"%((rel_dis.nbytes)/(1024*1024*1024)))
    print("storage for com_vel in gb = %f"%((com_vel.nbytes)/(1024*1024*1024)))
    #mass_PMIT = np.append(mass,0).reshape(-1,1)
    MIT, w_npt = prinaxes_data(rel_dis, vel, np_prinaxes_mod, mass, np_MIinv)
    print("storage for w_npt in gb = %f"%((w_npt.nbytes)/(1024*1024*1024)))
    m = np.repeat(np.tile(mass.transpose(), int(len(vel[0])/len(mass))).transpose(), repeats=3, axis=1)
    print("memory details before correlation analysis")
    print(psutil.virtual_memory())
    MITd = np.diagonal(MIT).reshape(-1,3)
    mit = np.repeat(MITd, repeats=len(w_npt[0]), axis=0)
    tot_dos, tot_dos_val, ff = dos_data(vel, m, dt) # Computing total density of states
    print(tot_dos_val) 
    tr_dos, tr_dos_val, ff = dos_data(com_vel, m, dt) # Computing translational density of states
    print(tr_dos_val)
    rot_dos, rot_dos_val, ff = dos_data(w_npt, mit, dt) # Computing rotational density of states
    print(rot_dos_val)
    outd1 = open("tot_dos_bulk_water.txt", "w") 
    outd2 = open("tr_dos_bulk_water.txt", "w")
    outd3 = open("rot_dos_bulk_water.txt", "w")
    outd1.truncate(0) 
    outd2.truncate(0)
    outd3.truncate(0)    
    for i in range(len(tot_dos)): 
        print(np.float(ff[i]*1e-12),np.float(tot_dos[i]*1e12), file = outd1) 
    outd1.close() 
    for i in range(len(tr_dos)): 
        print(np.float(ff[i]*1e-12),np.float(tr_dos[i]*1e12), file = outd2) 
    outd2.close() 
    for i in range(len(rot_dos)): 
        print(np.float(ff[i]*1e-12),np.float(rot_dos[i]*1e12), file = outd3) 
    outd3.close()  


    mvac_with_IFFT_FFT_tvel = corr_3D(vel,m)
    mvac_with_IFFT_FFT_trvel = corr_3D(com_vel,m)
    mvac_with_IFFT_FFT_rotvel = corr_3D(w_npt,mit)
    mvac_axes= np.sum(mvac_with_IFFT_FFT_rotvel, axis=(1))
    mvac = np.sum(mvac_with_IFFT_FFT_rotvel, axis=(0,1))
    outd1 = open("rot_axis1_vacf.txt", "w") 
    outd2 = open("rot_axis2_vacf.txt", "w")
    outd3 = open("rot_axis3_vacf.txt", "w")
    outd4 = open("rot_sum_all_axes.txt", "w")
    outd1.truncate(0) 
    outd2.truncate(0)
    outd3.truncate(0) 
    outd4.truncate(0)
    for i in range(len(mvac_axes[0,:])): 
        print(np.float(time[i]*1e12),np.float(mvac_axes[0,:][i]*(Na/1e3)), file = outd1) 
    outd1.close() 
    for i in range(len(mvac_axes[1,:])): 
        print(np.float(time[i]*1e12),np.float(mvac_axes[1,:][i]*(Na/1e3)), file = outd2) 
    outd2.close() 
    for i in range(len(mvac_axes[2,:])): 
        print(np.float(time[i]*1e12),np.float(mvac_axes[2,:][i]*(Na/1e3)), file = outd3) 
    outd3.close() 
    for i in range(len(mvac)): 
        print(np.float(time[i]*1e12),np.float(mvac[i]*(Na/1e3)), file = outd4) 
    outd4.close() 

    time_stop = perf_counter() # stopping time counter
    run_time = time_stop - time_start # run_time
    print(run_time)

