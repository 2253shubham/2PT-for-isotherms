# python ../../dos_H2O.py ../../2PT_H2O_files/350/RUN2/analysis.tpr ../../2PT_H2O_files/350/RUN2/prod1.trr ../../2PT_H2O_files/350/RUN2/prod1.xtc RUN2-output.txt "-760565.0" "RUN2" ../ref_data/avg_ref_tr_dos.txt ../ref_data/avg_ref_data.txt

import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
import sys
import math
import scipy.constants
from scipy.fftpack import fft, ifft
import MDAnalysis as mda
import sys
from scipy import integrate
from scipy import optimize
from time import perf_counter
from scipy.optimize import curve_fit
from numba import jit
import psutil

# Defining Consatants
kb = scipy.constants.k  # Boltzmann Constant k
h = scipy.constants.Planck  # Planck's constant
Na = scipy.constants.Avogadro  # Avogadro's number
T = 300  # Temperature, will make sure to make this as user input
beta = 1 / (kb * T)  # beta


def aid(x):
    # This function returns the memory
    # block address of an array.
    return x.__array_interface__["data"][0]


# extracting atom weight and storing them in SI units (kg/particle)
def extract_atom_weight(b):
    mass_dict = dict({i + 1: b[i].mass * (1e-3) / Na for i in range(len(b))})
    mass = np.array([x for x in mass_dict.values()])
    return mass.reshape(len(b), 1)  # reshaping them in 4x1 array


# extacting velocity and time and storing them in SI units (ms-1, s)


def extract_velocity_and_time(univ, x):
    print(psutil.virtual_memory())
    vel = np.zeros((len(univ), len(x), 3))
    time = np.zeros((len(vel), 1))
    print("after preallocation of vel and time arrays")
    print(psutil.virtual_memory())
    for i in range(len(univ)):
        univ[i]
        vel[i] = x.velocities
        time[i] = univ[i].time
    print("after extracting vel and time data")
    print(psutil.virtual_memory())
    return vel * 1e2, (time.reshape(-1)) * (
        1e-12
    )  # vel -> (40001x100x3), time -> (40001x1)


# extracting positions (m), c.o.m (m), MI (kg*m^2), MI inverse, and principal moments of inertia


def extract_pos_com_MI_prinaxes(univ, x, y, red):
    print(psutil.virtual_memory())
    l = np.zeros(
        (len(univ), len(red) * 3, 3)
    )  # list to store moment of inertia of all molecules at all time steps)
    pos = np.zeros(
        (len(univ), len(red.atoms), 3)
    )  # array to store positions of all atoms
    print("after preallocation of pos and MI arrays")
    print(psutil.virtual_memory())
    for i in range(len(univ)):
        univ[i]
        q = np.zeros(
            (len(red), 3, 3)
        )  # array to store moment of inertia of all molecules at one time step
        for j in range(len(red)):
            q[j] = red[j].atoms.moment_of_inertia()  # moment of inertia
        l[i] = q.reshape(-1, 3)
        pos[i] = red.atoms.positions
    print("after extracting pos and MI data")
    print(psutil.virtual_memory())
    return pos * (1e-10), l * (1e-20) / (1e3 * Na)


# calculate c.o.m. information from extracted positions and velocities and store them in SI units
def com_data(pos, vel, mass):
    tot_mass = np.sum(mass)  # total mass of a molecule
    m = np.tile(mass.transpose(), int(len(vel[0]) / len(mass))).transpose()
    m_v = np.zeros((len(vel), len(vel[0]), 3))
    m_p = np.zeros((len(pos), len(pos[0]), 3))
    for i in range(len(vel)):
        m_v[i] = vel[i] * m / tot_mass
        m_p[i] = pos[i] * m / tot_mass
    # m_v = vel*m/tot_mass
    # m_p = pos*m/tot_mass
    print("after m_v and m_p calculation for c.o.m")
    print(psutil.virtual_memory())
    return m_p.reshape(-1, len(mass), 3).sum(1).reshape(len(pos), -1, 3), m_v.reshape(
        -1, len(mass), 3
    ).sum(1).reshape(len(vel), -1, 3)


# calculate principal moments of inertia from Moment of Inertia Tensor (MIT)


def prinaxes(MI, x):
    print(psutil.virtual_memory())
    v, w = np.linalg.eig(
        MI.reshape(-1, 3, 3)
    )  # calculating eigen values and eigen vectors for the MIT inputs
    a1 = np.argsort(
        -v
    )  # sorting eigen values in desending order and storing the indices from largest to smallest eigen values for each MIT
    a3 = np.zeros(
        (len(w), 3, 3)
    )  # list to store eigen vectors (column wise) corresponding to each eigen values
    a4 = np.zeros(
        (len(v), 3, 3)
    )  # list to store eigen values for each MIT in descending order
    for i in range(len(a1)):
        a3[i] = w[i][
            :, a1[i]
        ]  # storing eigen vectors corresponding to decreasing eigen values
        a4[i] = np.diag(
            (v[i][a1[i]])
        )  # storing Principal Moments of Inertia (PMIT) which is the diagonal matrix of eigen values
    a3 = a3.reshape(x, -1, 3)  # shape -> 40001x75x3
    a4 = a4.reshape(x, -1, 3)  # shape -> 40001x75x3
    a4inv = np.linalg.inv(a4.reshape(-1, 3, 3)).reshape(
        x, -1, 3
    )  # calculating inverse of PMIT
    v = []
    w = []
    print("after computing MI data using numpy")
    print(psutil.virtual_memory())
    return a4, a4inv, a3


# calculating data along the principal axes of rotation (Xn = X*P', where Xn -> transformed X and P' -> the transformtion matrix with each columns being the principal axes of rotation in our case)


def data_along_prinaxes(pos, prinaxes):
    pos = pos.reshape(-1, 4, 3)
    prinaxes = prinaxes.reshape(-1, 3, 3)
    pr_pos = np.zeros((len(pos), 4, 3))
    print("after preallocation of pr_axis transformed rel_pos")
    print(psutil.virtual_memory())
    for i in range(len(pos)):
        pr_pos[i] = np.matmul(pos[i], prinaxes[i])
    return pr_pos.reshape(len(pos), -1, 3)


"""
def data_along_prinaxes(pos, prinaxes):
    pr_pos = np.matmul(pos.reshape(-1,4,3),prinaxes.reshape(-1,3,3)).reshape(len(pos),-1,3)
    return pr_pos
"""


# calculating rotational velocity from difference between total velocity and centre of mass velocity
def rot_velocity(vel1, vel2):
    rvel = vel1 - vel2
    return rvel


# calculating relative distance and velocities for each atom of a molecule w.r.t. its c.o.m. parameters
def rel_vector(comp, pos):
    rel_dis = pos - comp  # relative position
    return rel_dis


# calculating MIT from conventional formula


def Ixx_cal(vec, x, y, m, v):
    cal = np.zeros((len(vec), int(len(vec[0]) / v)))
    print("after preallocation of sub MI array")
    print(psutil.virtual_memory())
    for i in range(len(vec)):
        cal[i] = np.sum(
            (m * (vec[i][:, x] ** 2 + vec[i][:, y] ** 2)).reshape(-1, v), axis=1
        )
    return cal


def Ixy_cal(vec, x, y, m, v):
    cal = np.zeros((len(vec), int(len(vec[0]) / v)))
    print("after preallocation of sub MI array")
    print(psutil.virtual_memory())
    for i in range(len(vec)):
        cal[i] = np.sum((-m * (vec[i][:, x] * vec[i][:, y])).reshape(-1, v), axis=1)
    return cal


# calculating MIT from conventional formula
def moment_of_inertia(mass, vec):
    m = np.tile(mass.transpose(), int(len(vec[0]) / len(mass)))
    xxu = Ixx_cal(vec, 1, 2, m, len(mass))  # I_xx
    yyu = Ixx_cal(vec, 0, 2, m, len(mass))  # I_yy
    zzu = Ixx_cal(vec, 0, 1, m, len(mass))  # I_zz
    xyu = Ixy_cal(vec, 0, 1, m, len(mass))  # I_xy = I_yx
    xzu = Ixy_cal(vec, 0, 2, m, len(mass))  # I_xz = I_zx
    yzu = Ixy_cal(vec, 1, 2, m, len(mass))  # I_yz = I_zy
    MI = np.zeros((len(xxu), len(xxu[0]) * 3, 3))
    MIinv = np.zeros((len(xxu), len(xxu[0]) * 3, 3))
    print("after preallocation of MI and MIinv array")
    print(psutil.virtual_memory())
    for i in range(len(xxu)):
        MIts = np.zeros((len(xxu[0]), 3, 3))
        MItsinv = np.zeros((len(xxu[0]), 3, 3))
        for j in range(len(xxu[0])):
            k = np.array(
                [
                    [xxu[i, j], xyu[i, j], xzu[i, j]],
                    [xyu[i, j], yyu[i, j], yzu[i, j]],
                    [xzu[i, j], yzu[i, j], zzu[i, j]],
                ]
            )
            MIts[j] = k
            MItsinv[j] = np.linalg.inv(k)
        MIts = MIts.reshape(len(xxu[0]) * 3, 3)
        MItsinv = MItsinv.reshape(len(xxu[0]) * 3, 3)
        MI[i] = MIts
        MIinv[i] = MItsinv
    print("after computation of MI and MIinv array")
    print(psutil.virtual_memory())
    return MI, MIinv


"""
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
    yz = -m*(vec[:,:,1]*vec[:,:,2]) # I_yz = I_zy
    yzu = yz.reshape(-1,len(mass)).sum(1).reshape(len(vec),-1)
    MI = np.zeros((len(xxu),len(xxu[0])*3,3))
    for i in range(len(xxu)):
        MIts = np.zeros((len(xxu[0]),3,3))
        for j in range(len(xxu[0])):
            k = np.array([[xxu[i,j],xyu[i,j],xzu[i,j]],[xyu[i,j],yyu[i,j],yzu[i,j]],[xzu[i,j],yzu[i,j],zzu[i,j]]])
            MIts[j] = k
        MIts = MIts.reshape(len(xxu[0])*3,3)
        MI[i] = MIts
    print("storage for MI in gb = %f"%(MI.nbytes/(1024*1024*1024)))
    MIinv = np.linalg.inv(MI.reshape(-1,3,3)).reshape(len(vec),-1,3) # calculating inverse of MIT 
    return MI, MIinv
"""


# calculating angular momentum (m*(r x v))
def angular_momentum(mass, vec, vel):
    m = np.tile(mass.transpose(), int(len(vel[0]) / len(mass))).transpose()
    m_v = np.zeros((len(vel), len(vel[0]), 3))
    for i in range(len(vel)):
        k = vel[i] * m
        m_v[i] = k
    # m_v = vel*m
    print("after calculating mass*velocity array")
    print(psutil.virtual_memory())
    print("storage for m_v in gb = %f" % (m_v.nbytes / (1024 * 1024 * 1024)))
    L_i = np.zeros((len(vel), len(vel[0]), 3))
    print("after preallocation of L(ang. mom.) array")
    print(psutil.virtual_memory())
    for i in range(len(vel)):
        L_i[i] = np.cross(vec[i], m_v[i])
    return L_i.reshape(-1, len(mass), 3).sum(1).reshape(len(vel), -1, 3)


# @jit(nopython=True)
"""
def angular_momentum(mass, vec, vel):
    m = np.tile(mass.transpose(), int(len(vel[0])/len(mass))).transpose()
    m_v = vel*m
    print("storage for m_v in gb = %f"%(m_v.nbytes/(1024*1024*1024)))
    L_i = np.cross(vec, m_v)
    return L_i.reshape(-1,len(mass),3).sum(1).reshape(len(vel),-1,3) 
"""


# calculating angular velocity from angular momentum (MITinv * L)
def rot_angular_velocity(MIinvr, L, k):
    w = np.matmul(MIinvr.reshape(-1, 3, 3), L.reshape(-1, 3, 1)).reshape(k, -1, 1)
    return w


# calculating rotational velocity from angular velocity (w x r)
def rot_from_ang_velocity(w, d, x, y):
    v = np.cross(w.reshape(-1, 1, 3), d.reshape(-1, x, 3)).reshape(y, -1, 3)
    return v


# PMIT from literature:
def PMIT(mass):
    lpos = np.array(
        [[[0, 0, 0], [0.9572, 0, 0], [-0.23999, 0.92663, 0], [3, 4, 5]]]
    )  # dummy coordinates set according to data for Tip4p mwater model (r(OH) = 0.9572 angs, HOH angle = 104.52 deg)
    tot_mass = np.sum(mass)  # total mass of a molecule
    # part to calculate c.o.m.
    m_p = lpos * mass / tot_mass
    lcom = m_p.reshape(-1, len(mass), 3).sum(1).reshape(len(lpos), -1, 3)
    lcom = np.repeat(lcom, repeats=len(mass), axis=1)  # c.o.m data
    rpos = (lpos - lcom) * (1e-10)  # relative positions of atoms from c.o.m.
    lmit, lmiti = moment_of_inertia(mass, rpos)  # MIT tensor
    plmit, plmiti, pla = prinaxes(lmit, len(lpos))  # principal MIT
    return plmit, plmiti


# calculating autocorrelation values
def auto_corr_fft(vel):
    n = len(vel)
    norm_fac = np.arange(n, 0, -1, dtype=np.float_)
    c = np.hstack((vel, np.zeros(n)))
    c = np.real(ifft(np.power(np.absolute(fft(c)), 2)))
    c = c[:n] / norm_fac
    return c[: n // 2]


# calculating mass weighted autocorrelation values
def corr_3D(vel, m):
    corr = [[], [], []]
    for dim in range(3):
        ma = m[:, dim]
        for atom in range(vel.shape[1]):
            corr[dim].append(auto_corr_fft(vel[:, atom, dim]))
            corr[dim][atom] = corr[dim][atom] * ma[atom]
    corr = np.array(corr, dtype=np.float_)
    return corr


# calculating total autocorrelation
def total_autocorr(corr):
    total_ac = np.sum(corr, axis=(0, 1))
    return total_ac


# calculating DoS from total auto-correlation values
def density_of_states(ac):
    # ac = (ac/ac[0])*3*N*kb*T
    ac = np.hstack((ac, ac[::-1]))
    length = len(ac)
    yfft = np.empty(length, dtype=np.float_)
    yfft = fft(ac)
    yfft = np.real(yfft)
    # yfft = signal.savgol_filter(yn,21,5)
    return yfft


# calculating Total DoS
def integ_DoS(y, ff, dt):
    y = y * (beta) * 2 * dt
    integ = integrate.simps(y[: (int)((len(y)) / 2)].T, ff[: (int)((len(y)) / 2)].T)
    return y, np.float(integ)


# all required box frame data
def box_frame_data(pos, vel, mass, k):
    c_pos, c_vel = com_data(pos, vel, mass)  # c.o.m. positions and velocity
    com_pos = np.repeat(c_pos, repeats=k, axis=1)
    com_vel = np.repeat(c_vel, repeats=k, axis=1)
    rel_dis = rel_vector(com_pos, pos)
    c_pos = []
    c_vel = []
    return rel_dis, com_vel


# all principal axes frame data
def prinaxes_data(rel_dis, vel, np_prinaxes, mass, np_MIinv):
    # rel_dispr = data_along_prinaxes(rel_dis, np_prinaxes) # translating relative positions and velocities along principal axes frame
    # MI, MIinv =  moment_of_inertia(mass, rel_dispr) # MI calculated in pricipal axes frame (should be a diagonal matrix)
    # print("storage for MIinv in gb = %f"%((MIinv.nbytes)/(1024*1024*1024)))
    L = angular_momentum(mass, rel_dis, vel)  # calculating angular momnetum
    print("storage for L in gb = %f" % ((L.nbytes) / (1024 * 1024 * 1024)))
    Lt = (
        (np.transpose(np_prinaxes.reshape(-1, 3, 3), (0, 2, 1)))
        @ (np.transpose(L.reshape(-1, 1, 3), (0, 2, 1)))
    ).reshape(
        len(vel), -1, 3
    )  # calculating angular momentum in principal axes frame
    print("storage for Lt in gb = %f" % ((Lt.nbytes) / (1024 * 1024 * 1024)))
    L = []
    print("before calculating angular velocity")
    print(psutil.virtual_memory())
    w_np = rot_angular_velocity(
        np_MIinv, Lt.reshape(len(vel), -1), len(vel)
    )  # calculating angular velocity
    print("before calculating angular velocity")
    print(psutil.virtual_memory())
    aMIT, aMITi = PMIT(mass)
    aMIT = aMIT.reshape(-1, 3)
    w_np = w_np.reshape(len(vel), -1, 3)
    return aMIT, w_np


# dos data
def dos_data(x, y, ff, dt):
    dos_pr = density_of_states(
        total_autocorr(corr_3D(x, y))
    )  # calculating density of states with data in principal axes frame
    dos_pr = dos_pr.reshape(len(dos_pr), 1)
    dos_req, tdos = integ_DoS(
        dos_pr, ff[: len(ff) - 1], dt
    )  # calculating total density of states
    return dos_req, tdos


# all universe data
def extract_data_from_universe(u1, u2, x, y):
    vel, time = extract_velocity_and_time(u1, x)  # extracting velocities and time
    print("vel read")
    pos, md_MI = extract_pos_com_MI_prinaxes(u2, len(vel), len(y), y)
    print("pos read")
    print(md_MI.shape)
    time = time.reshape(len(time), 1)
    mass = extract_atom_weight(y[0].atoms)
    tot_mass = float(sum(mass))
    print(tot_mass)
    N = len(y)  # storing number of molecules
    dt = np.float((time[1] - time[0]))  # length of each time step
    ff = np.arange(len(time), dtype=np.float_) / (len(time) * dt)  # storing frequency
    ff = ff.reshape(len(time), 1)
    return N, md_MI, vel, time, pos, mass, ff, dt, tot_mass


# rotational fluidicity
def rot_hs_fluidicity(m, den, y0, T, k):
    delta = (
        math.sqrt((np.pi * kb * T) / m)
        * ((6 / np.pi) ** (2 / 3))
        * (den ** (1 / 3))
        * (2 * y0)
        / (9 * k)
    )  # normalized diffusivity constant

    def fluidicity(f):
        return (
            2 * (delta ** (0 - 9 / 2)) * (f ** (15 / 2))
            - (delta ** (0 - 3 / 2)) * (f ** (7 / 2))
            - 6 * (f ** (5)) * (delta ** (0 - 3))
            + 2 * f
            + 6 * (delta ** (0 - 3 / 2)) * (f ** (5 / 2))
            - 2
        )

    f_value = optimize.root_scalar(fluidicity, bracket=[0, 1], method="brentq")
    frac = f_value.root
    return frac


# rotational dos
def rotational_dos(fl, tdos, y, ff):
    p = (np.pi * y[0] * ff) / (2 * fl * tdos)
    yg = y[0] / (1 + (p**2))  # calculating gas DOS function vs frequency
    ys = y - yg  # calculating solid DOS function vs frequency
    return yg, ys


# hard sphere rotational entropy
def hs_rot_entropy_wf():
    r = (
        (
            (((np.pi) * (np.e) ** 3) * (T**3))
            / (13.706211107457026 * 20.998505222742857 * 39.467682045442515)
        )
        ** (0.5)
    ) * 0.5
    si = kb * np.log(r)  # entropy value per molecule basis
    return si / (3 * kb)


# entropy
def entropy(W_Shs, ff, y, ys, yg):
    ygm = W_Shs * yg[1 : (int)((len(y)) / 2)]
    W_Scho = 1 - np.log(beta * h * ff[1 : (int)((len(y)) / 2)])
    ysm = W_Scho * ys[1 : (int)((len(y)) / 2)]
    S_gas = (integrate.simps(ygm.T, ff[1 : (int)((len(y)) / 2)].T)) * kb  # Gas Entropy
    S_cho1 = (
        integrate.simps(ysm.T, ff[1 : (int)((len(y)) / 2)].T)
    ) * kb  # Solid Entropy (assuming classical harmonic oscillator)
    S_cho = S_cho1 + S_gas  # Entropy (assuming classical harmonic oscillator)
    return np.float(S_cho)


# internal energy
def energy(W_Ehs, ff, y, ys, yg):
    E_cho = (
        integrate.simps(
            yg[1 : (int)((len(y)) / 2)].T * W_Ehs, ff[1 : (int)((len(y)) / 2)].T
        )
        + integrate.simps(ys[1 : (int)((len(y)) / 2)].T, ff[1 : (int)((len(y)) / 2)].T)
    ) * (
        beta ** (-1)
    )  # Energy (assuming classical harmonic oscillator)
    return np.float(E_cho)


# helmholtz energy
def helmholtz(W_Ahs, ff, y, ys, yg):
    W_Acho = np.log(beta * h * ff[1 : (int)((len(y)) / 2)])
    ysm = W_Acho * ys[1 : (int)((len(y)) / 2)]
    ygm = W_Ahs * yg[1 : (int)((len(y)) / 2)]
    A_cho = (
        integrate.simps(ygm.T, ff[1 : (int)((len(y)) / 2)].T)
        + (integrate.simps(ysm.T, ff[1 : (int)((len(y)) / 2)].T))
    ) * beta ** (
        -1
    )  # Helmholtz free energy (assuming classical harmonic oscillator)
    return np.float(A_cho)


# --------------------------------------------------------------------------------------------------------------------

if __name__ == "__main__":
    time_start = perf_counter()  # starting time counter

    analysis_tpr_file = sys.argv[1]  # system argument to specify analysis.tpr file
    trr_file = sys.argv[2]  # system argument to specify .trr file
    xtc_file = sys.argv[3]  # system argument to specify .xtc file
    outfile = sys.argv[4]
    enerfile = sys.argv[5]
    runnum = sys.argv[6]
    reffile = sys.argv[
        7
    ]  # system argument to specify translational gas reference dos file
    refpropfile = sys.argv[
        8
    ]  # system argument to specify translational reference state properties file

    print("files read")

    # system information
    md_energy = np.float(enerfile)
    print(outfile)
    print(md_energy)
    u1 = mda.Universe(analysis_tpr_file, trr_file)  # creating universe with velocities
    print(u1)
    u2 = mda.Universe(analysis_tpr_file, xtc_file)  # creating universe with positions
    print(u2)
    wat_atoms = u1.select_atoms("resname WAT")
    print(wat_atoms)
    vol = 4.00440 * 3.97980 * 4.01490 * (1e-27)  # volume of the system
    den = 1 / vol  # calculating density
    p1 = u1.trajectory[0 : len(u1.trajectory) : 10]
    print(p1)
    p2 = u2.trajectory[0 : len(u2.trajectory) : 10]
    print(p2)
    N, md_MI, vel, time, pos, mass, ff, dt, tot_mass = extract_data_from_universe(
        p1, p2, wat_atoms, u2.residues
    )
    print(N)
    print(vel.shape)
    print("storage for vel in gb = %f" % ((vel.nbytes) / (1024 * 1024 * 1024)))
    print("storage for md_MI in gb = %f" % ((md_MI.nbytes) / (1024 * 1024 * 1024)))
    np_MI, np_MIinv, np_prinaxes = prinaxes(md_MI, len(vel))
    print(np_prinaxes.shape)
    print("storage for np_MI in gb = %f" % ((np_MI.nbytes) / (1024 * 1024 * 1024)))
    print(
        "storage for np_prinaxes in gb = %f"
        % ((np_prinaxes.nbytes) / (1024 * 1024 * 1024))
    )
    rel_dis, com_vel = box_frame_data(pos, vel, mass, len(u1.residues[0].atoms))
    print("storage for rel_dis in gb = %f" % ((rel_dis.nbytes) / (1024 * 1024 * 1024)))
    print("storage for com_vel in gb = %f" % ((com_vel.nbytes) / (1024 * 1024 * 1024)))
    MIT, w_npt = prinaxes_data(rel_dis, vel, np_prinaxes, mass, np_MIinv)
    print("storage for w_npt in gb = %f" % ((w_npt.nbytes) / (1024 * 1024 * 1024)))
    m = np.repeat(
        np.tile(mass.transpose(), int(len(vel[0]) / len(mass))).transpose(),
        repeats=3,
        axis=1,
    )
    MITd = np.diagonal(MIT).reshape(-1, 3)
    mit = np.repeat(MITd, repeats=len(w_npt[0]), axis=0)
    tot_corr = total_autocorr(corr_3D(vel, m))
    tot_dos, tot_dos_val = dos_data(vel, m, ff, dt)  # Computing total density of states
    outc = open("tot_corr_per_10_for_" + runnum + ".txt", "w")
    outd = open("dos_per_10_for_" + runnum + ".txt", "w")
    outc.truncate(0)
    outd.truncate(0)
    for i in range(len(tot_corr)):
        print(np.float(time[i]), np.float(tot_corr[i]), file=outc)
    outc.close()
    for i in range(len(tot_dos)):
        print(np.float(ff[i]), np.float(tot_dos[i]), file=outd)
    outd.close()
    print(tot_dos_val)
    tr_dos, tr_dos_val = dos_data(
        com_vel, m, ff, dt
    )  # Computing translational density of states
    print(tr_dos_val)
    rot_dos, rot_dos_val = dos_data(
        w_npt, mit, ff, dt
    )  # Computing rotational density of states
    print(rot_dos_val)

    # Translational gas reference system data
    tr_dos_ref = np.loadtxt(fname=reffile)
    dataref = np.loadtxt(fname=refpropfile)
    tr_dif_sp_ref = np.float(
        dataref[0]
    )  # value of self-diffusion coefficient for reference state
    tr_e_gas_ref = dataref[1]  # value of reference gas md energy
    tr_s_gas_ref = dataref[2]  # value of reference gas entropy
    a = dataref[4]  # value of curve parameter for Ideal Helmholtz energy
    b = dataref[5]  # value of curve parameter for Ideal Helmholtz energy
    c = dataref[6]  # value of curve parameter for Ideal Helmholtz energy
    d = dataref[7]  # value of curve parameter for Ideal Helmholtz energy
    e = dataref[8]  # value of curve parameter for Ideal Helmholtz energy

    # Total system data
    # md_energy = (-244.025/14.3771)*1e3/Na + 3*kb*T # total md energy
    md_energy = (
        np.float(enerfile) * 1e3 / Na + 750629 * 1e3 / Na
    )  # total md energy*1e3/Na + 750629*1e3/Na
    tot_u_ex_ref = dataref[3]  # value of reference excess gibbs free energy
    tot_u_ex = N * tot_u_ex_ref  # excess gibbs free energy

    # Rotational system data
    rot_dif_sp = np.float(
        rot_dos[0] * (T * kb) / (12 * tot_mass * 1)
    )  # Rotational Self diffusion coefficient (p)
    rot_fl = rot_hs_fluidicity(tot_mass, den, rot_dos[0], T, N)  # Rotational fluidicity
    print(rot_fl)
    gas_rot_dos, solid_rot_dos = rotational_dos(
        rot_fl, rot_dos_val, rot_dos, ff[0 : len(ff) - 1]
    )  # gas and solid rotational dos
    W_Shs = hs_rot_entropy_wf()  # hs approx. rotational entropy weighing function
    W_Ahs = 0.5 - W_Shs  # hs approx. rotational helmholtz energy weighing function
    W_Ehs = 0.5  # hs approx. rotational internal energy weighing function
    rot_s = entropy(
        W_Shs, ff, rot_dos, solid_rot_dos, gas_rot_dos
    )  # rotational part entropy
    rot_e = energy(
        W_Ehs, ff, rot_dos, solid_rot_dos, gas_rot_dos
    )  # rotational part internal energy
    rot_h = helmholtz(
        W_Ahs, ff, rot_dos, solid_rot_dos, gas_rot_dos
    )  # rotational part helmholtz energy
    rot_ln_q = -rot_h * beta  # ln of rotational partition function
    print(rot_s)

    # Translational system data
    tr_dif_sp = np.float(
        tr_dos[0] * (T * kb) / (12 * tot_mass * N)
    )  # Translational Self diffusion coefficient (p)
    tr_fl = tr_dif_sp / tr_dif_sp_ref  # fluidicity
    gas_tr_dos = (
        tr_fl * tr_dos_ref.reshape(-1, 1) * N
    )  # calculating gas DOS function vs frequency
    solid_tr_dos = tr_dos - gas_tr_dos  # calculating solid DOS function vs frequency
    tr_u_ex = tot_u_ex  # excess translational gibbs energy
    tr_h_gas = tr_fl * tr_u_ex + (
        a
        + b * (N * tr_fl)
        + c * ((N * tr_fl) ** 2)
        + d * ((N * tr_fl) ** 3)
        + e * ((N * tr_fl) ** 4)
    )
    tr_h = tr_h_gas + helmholtz(
        0, ff, tr_dos, solid_tr_dos, gas_tr_dos
    )  # translational part entropy
    tr_e = N * tr_fl * tr_e_gas_ref + energy(
        0, ff, tr_dos, solid_tr_dos, gas_tr_dos
    )  # translational part internal energy
    tr_s = (N * tr_fl * tr_e_gas_ref - tr_h_gas) / T + entropy(
        0, ff, tr_dos, solid_tr_dos, gas_tr_dos
    )  # translational part helmholtz energy
    tr_ln_q = -tr_h * beta  # ln of translational partition function

    # Total system data
    reference = (
        md_energy
        - ((3 * N) * (2 - 0.5 * rot_fl - tr_fl)) / (beta)
        - N * tr_fl * tr_e_gas_ref
    )  # reference energy
    tot_e = reference + tr_e + rot_e  # total internal energy
    tot_h = reference + tr_h + rot_h  # total helmholtz free energy
    tot_s = tr_s + rot_s  # total entropy
    tot_ln_q = -tot_h * beta  # ln of total partition function

    time_stop = perf_counter()  # stopping time counter
    run_time = time_stop - time_start  # run_time
    print(run_time)

    out = open(outfile, "w")
    out.truncate(0)
    print("Data collected for N = " + str(N) + " molecules system", file=out)
    print("All values are in their SI units \n", file=out)
    print("\n Total system data :\n", file=out)
    print("Total Density of states in the system = ", tot_dos_val, file=out)
    print(
        "MD energy of the system scales with N (normalized by beta) = ",
        md_energy / (kb * T),
        file=out,
    )
    print(
        "Total Helmholtz free energy of the system scales with N (normalized by beta) = ",
        tot_h / (kb * T),
        file=out,
    )
    print(
        "Total Internal energy of the system scales with N (normalized by beta) = ",
        tot_e / (kb * T),
        file=out,
    )
    print("Total Entropy (scaled by kb) = ", tot_s / kb, file=out)
    print(
        "Reference potential energy of the system scales with N (normalized by beta) = ",
        reference / (kb * T),
        file=out,
    )
    print(
        "ln of Total Partition function of the system (assuming classical harmonic oscillator for solid part of the system) scales with N = ",
        tot_ln_q,
        file=out,
    )
    print("\n\n Rotational System Data: \n", file=out)
    print("Total Rotational Density of states in the system = ", rot_dos_val, file=out)
    print("Rotational Fluidicity / Fraction of gas = ", rot_fl, file=out)
    print(
        "Total number of gas Rotational DOS in the system = ",
        rot_fl * rot_dos_val,
        file=out,
    )
    print(
        "Total number of solid Rotational DOS in the system = ",
        (1 - rot_fl) * rot_dos_val,
        file=out,
    )
    print(
        "Total Rotational Helmholtz free energy of the system scales with N (normalized by beta) = ",
        rot_h / (kb * T),
        file=out,
    )
    print(
        "Total Rotational energy of the system scales with N (normalized by beta) = ",
        rot_e / (kb * T),
        file=out,
    )
    print("Total Rotational Entropy (scaled by kb) = ", rot_s / kb, file=out)
    print(
        "ln of Rotational Partition function of the system (assuming classical harmonic oscillator for solid part of the system) scales with N = ",
        rot_ln_q,
        file=out,
    )
    print("\n\n Translational System Data: \n", file=out)
    print(
        "Total Translational Density of states in the system = ", tr_dos_val, file=out
    )
    print("Translational Fluidicity / Fraction of gas = ", tr_fl, file=out)
    print(
        "Total number of gas Translational DOS in the system = ",
        tr_fl * tr_dos_val,
        file=out,
    )
    print(
        "Total number of solid Translational DOS in the system = ",
        (1 - tr_fl) * tr_dos_val,
        file=out,
    )
    print(
        "Total Translational Helmholtz free energy of the system scales with N (normalized by beta) = ",
        tr_h / (kb * T),
        file=out,
    )
    print(
        "Total Translational energy of the system scales with N (normalized by beta) = ",
        tr_e / (kb * T),
        file=out,
    )
    print("Total Translational Entropy (scaled by kb) = ", tr_s / kb, file=out)
    print(
        "ln of Translational Partition function of the system (assuming classical harmonic oscillator for solid part of the system) scales with N = ",
        tr_ln_q,
        file=out,
    )
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
