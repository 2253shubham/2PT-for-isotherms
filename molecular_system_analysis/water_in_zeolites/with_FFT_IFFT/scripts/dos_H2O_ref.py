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

# Defining Consatants
kb = scipy.constants.k  # Boltzmann Constant k
h = scipy.constants.Planck  # Planck's constant
Na = scipy.constants.Avogadro  # Avogadro's number
T = 300  # Temperature, will make sure to make this as user input
beta = 1 / (kb * T)  # beta


# extracting atom weight and storing them in SI units (kg/particle)
def extract_atom_weight(b):
    mass_dict = dict({i + 1: b[i].mass * (1e-3) / Na for i in range(len(b))})
    mass = np.array([x for x in mass_dict.values()])
    return mass.reshape(len(b), 1)  # reshaping them in 4x1 array


# extacting velocity and time and storing them in SI units (ms-1, s)
def extract_velocity_and_time(univ, x):
    vel = []
    time = []
    for i in univ:
        vel.append(x.velocities)
        time.append(i.time)
    return np.array(vel, dtype=np.float_) * 1e2, np.array(time, dtype=np.float_) * (
        1e-12
    )  # vel -> (40001x100x3), time -> (40001x1)


# extracting positions (m), c.o.m (m), MI (kg*m^2), MI inverse, and principal moments of inertia
def extract_pos_com_MI_prinaxes(univ, x, y, red):
    l = []  # list to store moment of inertia of all molecules at all time steps)
    pos = []  # list to store positions of all atoms
    for i in univ:
        q = []  # list to store moment of inertia of all molecules at one time step
        for j in red:
            q.append(j.atoms.moment_of_inertia())  # moment of inertia
        l.append(q)
        pos.append(red.atoms.positions)
    posa = np.array(pos, dtype=np.float_) * (
        1e-10
    )  # array to store positions in SI units
    la = (
        np.array(l, dtype=np.float_).reshape(x, y * 3, 3) * (1e-20) / (1e3 * Na)
    )  # array to store moment of inertia in SI units
    return posa, la


# calculate c.o.m. information from extracted positions and velocities and store them in SI units
def com_data(pos, vel, mass):
    tot_mass = np.sum(mass)  # total mass of a molecule
    m = np.tile(mass.transpose(), int(len(vel[0]) / len(mass))).transpose()
    m_v = vel * m / tot_mass
    m_p = pos * m / tot_mass
    return m_p.reshape(-1, len(mass), 3).sum(1).reshape(len(pos), -1, 3), m_v.reshape(
        -1, len(mass), 3
    ).sum(1).reshape(len(vel), -1, 3)


# calculate principal moments of inertia from Moment of Inertia Tensor (MIT)
def prinaxes(MI, x):
    v, w = np.linalg.eig(
        MI.reshape(-1, 3, 3)
    )  # calculating eigen values and eigen vectors for the MIT inputs
    a1 = np.argsort(
        -v
    )  # sorting eigen values in desending order and storing the indices from largest to smallest eigen values for each MIT
    a3 = (
        []
    )  # list to store eigen vectors (column wise) corresponding to each eigen values
    a4 = []  # list to store eigen values for each MIT in descending order
    for i in range(len(a1)):
        a3.append(
            w[i][:, a1[i]]
        )  # storing eigen vectors corresponding to decreasing eigen values
        a4.append(v[i][a1[i]])  # storing eigen values in descending order
    a3 = np.array(a3).reshape(x, -1, 3)  # shape -> 40001x75x3
    a4 = np.array(a4)
    diag = []
    for i in a4:
        diag.append(
            np.diag(i)
        )  # storing Principal Moments of Inertia (PMIT) which is the diagonal matrix of eigen values
    diaga = np.array(diag).reshape(x, -1, 3)  # shape -> 40001x75x3
    diagainv = np.linalg.inv(diaga.reshape(-1, 3, 3)).reshape(
        x, -1, 3
    )  # calculating inverse of PMIT
    return diaga, diagainv, a3


# calculating data along the principal axes of rotation (Xn = X*P', where Xn -> transformed X and P' -> the transformtion matrix with each columns being the principal axes of rotation in our case)
def data_along_prinaxes(pos, prinaxes):
    pr_pos = np.matmul(pos.reshape(-1, 4, 3), prinaxes.reshape(-1, 3, 3)).reshape(
        len(pos), -1, 3
    )
    return pr_pos


# calculating rotational velocity from difference between total velocity and centre of mass velocity
def rot_velocity(vel1, vel2):
    rvel = vel1 - vel2
    return rvel


# calculating relative distance and velocities for each atom of a molecule w.r.t. its c.o.m. parameters
def rel_vector(comp, pos, comv, vel):
    rel_dis = pos - comp  # relative position
    rel_vel = vel - comv  # relative velocities
    return rel_dis, rel_vel


# calculating MIT from conventional formula
def moment_of_inertia(mass, vec):
    m = np.tile(mass.transpose(), int(len(vec[0]) / len(mass)))
    xx = m * (vec[:, :, 1] ** 2 + vec[:, :, 2] ** 2)  # I_xx
    xxu = xx.reshape(-1, len(mass)).sum(1).reshape(len(vec), -1)
    yy = m * (vec[:, :, 0] ** 2 + vec[:, :, 2] ** 2)  # I_yy
    yyu = yy.reshape(-1, len(mass)).sum(1).reshape(len(vec), -1)
    zz = m * (vec[:, :, 0] ** 2 + vec[:, :, 1] ** 2)  # I_zz
    zzu = zz.reshape(-1, len(mass)).sum(1).reshape(len(vec), -1)
    xy = -m * (vec[:, :, 0] * vec[:, :, 1])  # I_xy = I_yx
    xyu = xy.reshape(-1, len(mass)).sum(1).reshape(len(vec), -1)
    xz = -m * (vec[:, :, 0] * vec[:, :, 2])  # I_xz = I_zx
    xzu = xz.reshape(-1, len(mass)).sum(1).reshape(len(vec), -1)
    yz = -m * (vec[:, :, 1] * vec[:, :, 2])  # I_yz = I_zy
    yzu = yz.reshape(-1, len(mass)).sum(1).reshape(len(vec), -1)
    MIl = []  # list to store MIT
    for i in range(len(xxu)):
        MIlts = []
        for j in range(len(xxu[0])):
            k = np.array(
                [
                    [xxu[i, j], xyu[i, j], xzu[i, j]],
                    [xyu[i, j], yyu[i, j], yzu[i, j]],
                    [xzu[i, j], yzu[i, j], zzu[i, j]],
                ]
            )
            MIlts.append(k)
        MIats = np.array(MIlts).reshape(len(xxu[0]) * 3, 3)
        MIl.append(MIats)
    MIa = np.array(MIl)
    MIainv = np.linalg.inv(MIa.reshape(-1, 3, 3)).reshape(
        len(vec), -1, 3
    )  # calculating inverse of MIT
    return MIa, MIainv


# calculating angular momentum (m*(r x v))
def angular_momentum(mass, vec, vel):
    m = np.tile(mass.transpose(), int(len(vel[0]) / len(mass))).transpose()
    m_v = vel * m
    L_i = np.cross(vec, m_v)
    return (
        L_i.reshape(-1, len(mass), 3).sum(1).reshape(len(vel), -1, 3)
    )  # shape -> 40001x25x3


# calculating angular velocity from angular momentum (MITinv * L)
def rot_angular_velocity(MIinv, L, k):
    w = np.matmul(MIinv.reshape(-1, 3, 3), L.reshape(-1, 3, 1)).reshape(
        k, -1, 1
    )  # shape -> 40001x75x1
    return w


# calculating rotational velocity from angular velocity (w x r)
def rot_from_ang_velocity(w, d, x, y):
    v = np.cross(w.reshape(-1, 1, 3), d.reshape(-1, x, 3)).reshape(
        y, -1, 3
    )  # shape -> 40001x100x3
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
def integ_DoS(y, ff, dt, k):
    y = y * (beta) * 2 * dt / k
    integ = integrate.simps(y[: (int)((len(y)) / 2)].T, ff[: (int)((len(y)) / 2)].T)
    return y, np.float(integ)


# all required box frame data
def box_frame_data(pos, vel, mass, k):
    c_pos, c_vel = com_data(pos, vel, mass)  # c.o.m. positions and velocity
    com_pos = np.repeat(c_pos, repeats=k, axis=1)
    com_vel = np.repeat(c_vel, repeats=k, axis=1)
    rot_vel = rot_velocity(vel, com_vel)  # calculating rotational velocity
    rel_dis, rel_vel = rel_vector(com_pos, pos, com_vel, vel)
    return rel_dis, com_vel


# all principal axes frame data
def prinaxes_data(rel_dis, vel, np_prinaxes, mass):
    rel_dispr = data_along_prinaxes(
        rel_dis, np_prinaxes
    )  # translating relative positions and velocities along principal axes frame
    MI, MIinv = moment_of_inertia(
        mass, rel_dispr
    )  # MI calculated in pricipal axes frame (should be a diagonal matrix)
    L = angular_momentum(mass, rel_dis, vel)  # calculating angular momnetum
    Lt = (
        (np.transpose(np_prinaxes.reshape(-1, 3, 3), (0, 2, 1)))
        @ (np.transpose(L.reshape(-1, 1, 3), (0, 2, 1)))
    ).reshape(
        len(vel), -1, 3
    )  # calculating angular momnetum in principal axes frame
    w_np = rot_angular_velocity(
        MIinv, Lt.reshape(len(vel), -1), len(vel)
    )  # calculating angular velocity
    aMIT, aMITi = PMIT(mass)
    aMIT = aMIT.reshape(-1, 3)
    w_npt = w_np.reshape(len(vel), -1, 3)
    return aMIT, w_npt


# dos data
def dos_data(x, y, ff, dt, k):
    dos_pr = density_of_states(
        total_autocorr(corr_3D(x, y))
    )  # calculating density of states with data in principal axes frame
    dos_pr = dos_pr.reshape(len(dos_pr), 1)
    dos_req, tdos = integ_DoS(
        dos_pr, ff[: len(ff) - 1], dt, k
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
def rot_hs_fluidicity(m, den, y0, T):
    delta = (
        math.sqrt((np.pi * kb * T) / m)
        * ((6 / np.pi) ** (2 / 3))
        * (den ** (1 / 3))
        * (2 * y0)
        / 9
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
def ref_rotational_dos(fl, tdos, y, ff):
    p = (np.pi * y[0] * ff) / (2 * fl * tdos)
    yg = y[0] / (1 + (p**2))  # calculating gas DOS function vs frequency
    ys = y - yg  # calculating solid DOS function vs frequency
    return yg, ys


# hard sphere rotational entropy
def hs_rot_entropy_wf(tot_mass, den, rot_fl, k):
    r = (
        (
            (((np.pi) * (np.e) ** 3) * (T**3))
            / (13.706211107457026 * 20.998505222742857 * 39.467682045442515)
        )
        ** (0.5)
    ) * 0.5
    si = kb * np.log(r)  # entropy value per molecule basis
    return si / (3 * kb)


# total rotational entropy
def rot_entropy(W_Shs, ff, y, ys, yg):
    ygm = W_Shs * yg[1 : (int)((len(y)) / 2)]
    W_Scho = 1 - np.log(beta * h * ff[1 : (int)((len(y)) / 2)])
    ysm = W_Scho * ys[1 : (int)((len(y)) / 2)]
    S_gas = (integrate.simps(ygm.T, ff[1 : (int)((len(y)) / 2)].T)) * kb  # Gas Entropy
    S_cho1 = (
        integrate.simps(ysm.T, ff[1 : (int)((len(y)) / 2)].T)
    ) * kb  # Solid Entropy (assuming classical harmonic oscillator)
    S_cho = S_cho1 + S_gas  # Entropy (assuming classical harmonic oscillator)
    return np.float(S_cho)


# rotational internal energy (assumption --> reference energy = 0)
def rot_energy(W_Ehs, ff, y, ys, yg):
    E_cho = (
        integrate.simps(
            yg[1 : (int)((len(y)) / 2)].T * W_Ehs, ff[1 : (int)((len(y)) / 2)].T
        )
        + integrate.simps(ys[1 : (int)((len(y)) / 2)].T, ff[1 : (int)((len(y)) / 2)].T)
    ) * (
        beta ** (-1)
    )  # Energy (assuming classical harmonic oscillator)
    return np.float(E_cho)


# rotational helmholtz energy (assumption --> reference energy = 0)
def rot_helmholtz(W_Ahs, ff, y, ys, yg):
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


def trans_iag_helmholtz(vol, m, T):
    A_ig = []
    Nr = np.arange(0, 200, 1)
    for i in Nr:
        x = -kb * T * i * (np.log(vol * ((2 * np.pi * m * kb * T) / (h**2)) ** (1.5)))
        for j in range(i):
            x = x + kb * T * np.log(j + 1)
        A_ig.append(x)
        x = 0

    def func(x, a, b, c, d, e):
        return a + b * x + c * (x**2) + d * (x**3) + e * (x**4)

    popt, pconv = curve_fit(func, Nr, A_ig)
    a = popt[0]
    b = popt[1]
    c = popt[2]
    d = popt[3]
    e = popt[4]

    return a, b, c, d, e


# --------------------------------------------------------------------------------------------------------------------

if __name__ == "__main__":
    time_start = perf_counter()  # starting time counter

    analysis_tpr_file = sys.argv[1]  # system argument to specify analysis.tpr file
    trr_file = sys.argv[2]  # system argument to specify .trr file
    xtc_file = sys.argv[3]  # system argument to specify .xtc file
    outfile1 = sys.argv[4]
    outfile2 = sys.argv[5]
    enerfile = sys.argv[6]
    runnum = sys.argv[7]

    print("files read")

    # system information
    md_energy = np.float(enerfile)
    print(outfile1)
    print(md_energy)
    u1_ref = mda.Universe(
        analysis_tpr_file, trr_file
    )  # creating universe with velocities
    print(u1_ref)
    u2_ref = mda.Universe(
        analysis_tpr_file, xtc_file
    )  # creating universe with positions
    print(u2_ref)
    wat_atoms = u1_ref.select_atoms("resname WAT")
    print(wat_atoms)
    vol = 4.00440 * 3.97980 * 4.01490 * (1e-27)  # volume of the system
    den = 1 / vol  # calculating density
    p1 = u1_ref.trajectory[0 : len(u1_ref.trajectory) : 10]
    print(p1)
    p2 = u2_ref.trajectory[0 : len(u2_ref.trajectory) : 10]
    print(p2)
    N_ref, md_MI_ref, vel_ref, time, pos_ref, mass_ref, ff, dt, tot_mass = (
        extract_data_from_universe(p1, p2, wat_atoms, u2_ref.residues)
    )
    print(N_ref)
    print(vel_ref.shape)
    np_MI_ref, np_MIinv_ref, np_prinaxes_ref = prinaxes(md_MI_ref, len(vel_ref))
    print(np_prinaxes_ref.shape)
    rel_dis_ref, com_vel_ref = box_frame_data(
        pos_ref, vel_ref, mass_ref, len(u1_ref.residues[0].atoms)
    )
    MIT_ref, w_npt_ref = prinaxes_data(rel_dis_ref, vel_ref, np_prinaxes_ref, mass_ref)
    m_ref = np.repeat(
        np.tile(mass_ref.transpose(), int(len(vel_ref[0]) / len(mass_ref))).transpose(),
        repeats=3,
        axis=1,
    )
    MITd_ref = np.diagonal(MIT_ref).reshape(-1, 3)
    mit_ref = np.repeat(MITd_ref, repeats=len(w_npt_ref[0]), axis=0)
    tot_corr = total_autocorr(corr_3D(vel_ref, m_ref))
    tot_dos_ref, tot_dos_val_ref = dos_data(
        vel_ref, m_ref, ff, dt, N_ref
    )  # Computing total density of states
    outc = open("tot_corr_per_10_" + runnum + ".txt", "w")
    outd = open("dos_per_10_" + runnum + ".txt", "w")
    outc.truncate(0)
    outd.truncate(0)
    for i in range(len(tot_corr)):
        print(np.float(time[i]), np.float(tot_corr[i]), file=outc)
    outc.close()
    for i in range(len(tot_dos_ref)):
        print(np.float(ff[i]), np.float(tot_dos_ref[i]), file=outd)
    outd.close()
    print(tot_dos_val_ref)
    tr_dos_ref, tr_dos_val_ref = dos_data(
        com_vel_ref, m_ref, ff, dt, N_ref
    )  # Computing translational density of states
    print(tr_dos_val_ref)
    rot_dos_ref, rot_dos_val_ref = dos_data(
        w_npt_ref, mit_ref, ff, dt, N_ref
    )  # Computing rotational density of states
    print(rot_dos_val_ref)

    # copying translational part dos to a file
    out = open(outfile1, "w")
    out.truncate(0)
    for i in range(len(tr_dos_ref)):
        print(np.float(tr_dos_ref[i]), file=out)
    out.close()

    # Total system data
    # md_energy = (-244.025/14.3771)*1e3/Na + 3*kb*T # total md energy
    md_energy = (
        np.float(enerfile) * 1e3 / Na + 750629 * 1e3 / Na
    )  # total md energy*1e3/Na + 750629*1e3/Na
    u_ex = -6.64901 * (1e3) / Na  # excess gibbs free energy
    r = (((2 * np.pi * tot_mass * kb * T) / (h**2)) ** (1.5)) * (vol / 1)
    u_ig = -kb * T * np.log(r)  # ideal gibbs free energy
    u = u_ig + u_ex  # total gibbs free energy
    s = -(u - md_energy) / T  # total entropy

    # Rotational system data
    rot_dif_sp_ref = np.float(
        rot_dos_ref[0] * (T * kb) / (12 * tot_mass * 1)
    )  # Rotational Self diffusion coefficient (p)
    rot_fl = rot_hs_fluidicity(
        tot_mass, den, rot_dos_ref[0], T
    )  # Rotational fluidicity
    print(rot_fl)
    gas_rot_dos_ref, solid_rot_dos_ref = ref_rotational_dos(
        rot_fl, rot_dos_val_ref, rot_dos_ref, ff[0 : len(ff) - 1]
    )  # gas and solid rotational dos
    W_Shs = hs_rot_entropy_wf(
        tot_mass, den, rot_fl, rot_dos_ref[0]
    )  # hs approx. rotational entropy weighing function
    W_Ahs = 0.5 - W_Shs  # hs approx. rotational helmholtz energy weighing function
    W_Ehs = 0.5  # hs approx. rotational internal energy weighing function
    rot_s = rot_entropy(
        W_Shs, ff, rot_dos_ref, solid_rot_dos_ref, gas_rot_dos_ref
    )  # rotational part entropy
    rot_e = rot_energy(
        W_Ehs, ff, rot_dos_ref, solid_rot_dos_ref, gas_rot_dos_ref
    )  # rotational part internal energy
    rot_h = rot_helmholtz(
        W_Ahs, ff, rot_dos_ref, solid_rot_dos_ref, gas_rot_dos_ref
    )  # rotational part helmholtz energy
    print(rot_s)

    # Translational system data
    reference = md_energy - 1.5 * (1 - rot_fl) * (beta**-1)
    tr_dif_sp_ref = np.float(
        tr_dos_ref[0] * (T * kb) / (12 * tot_mass * 1)
    )  # Translational Self diffusion coefficient (p)
    tr_s = s - rot_s  # translational part entropy
    tr_e = md_energy - rot_e - reference  # translational part internal energy
    tr_h = u - rot_h - reference  # translational part helmholtz energy
    a, b, c, d, e = trans_iag_helmholtz(vol, tot_mass, T)

    time_stop = perf_counter()  # stopping time counter
    run_time = time_stop - time_start  # run_time

    out = open(outfile2, "w")
    out.truncate(0)
    print(tr_dif_sp_ref, file=out)
    print(tr_e, file=out)
    print(tr_s, file=out)
    print(u_ex - rot_h, file=out)
    print(a, file=out)
    print(b, file=out)
    print(c, file=out)
    print(d, file=out)
    print(e, file=out)
    out.close()

    print(run_time)

"""
    plt.figure(1)
    plt.plot(ff[:len(ff)-1], tot_dos_ref, "k--", label="total DoS")
    plt.plot(ff[:len(ff)-1], tr_dos_ref, "r--", label="translational DoS")
    plt.plot(ff[:len(ff)-1], rot_dos_ref, "b--", label="rotational DoS")
    plt.xlabel("freq in s-1")
    plt.ylabel("S(v) in s")
    plt.legend()
    plt.show()

    plt.figure(2)
    plt.plot(ff[:len(ff)-1], rot_dos_ref, "k--", label="total rotational DoS") 
    plt.plot(ff[:len(ff)-1], solid_rot_dos_ref, "r--", label="rotational solid DoS") 
    plt.plot(ff[:len(ff)-1], gas_rot_dos_ref, "b--", label="rotational gas DoS") 
    plt.xlabel("freq in s-1") 
    plt.ylabel("S(v) in s") 
    plt.legend() 
    plt.show()  
"""
