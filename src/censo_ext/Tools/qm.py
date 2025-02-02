#!/usr/bin/env python3
import numpy as np
import os
import argparse
# from scipy.linalg import eigh
# from scipy.sparse import kron, csc_matrix, csr_matrix, lil_matrix, bmat
from icecream import ic


def qm_parameter(v: list[float], J: np.ndarray):

    sigma_x = np.array([[0, 1 / 2], [1 / 2, 0]])
    sigma_y = np.array([[0, -1j / 2], [1j / 2, 0]])
    sigma_z = np.array([[1 / 2, 0], [0, -1 / 2]])
    unit = np.array([[1, 0], [0, 1]])

    nspins: int = len(v)
    L = np.empty((3, nspins, 2 ** nspins, 2 ** nspins), dtype=np.complex128)
    for n in range(nspins):
        Lx_current = 1
        Ly_current = 1
        Lz_current = 1

        for k in range(nspins):
            if k == n:
                Lx_current = np.kron(Lx_current, sigma_x)
                Ly_current = np.kron(Ly_current, sigma_y)
                Lz_current = np.kron(Lz_current, sigma_z)
            else:
                Lx_current = np.kron(Lx_current, unit)
                Ly_current = np.kron(Ly_current, unit)
                Lz_current = np.kron(Lz_current, unit)

        L[0][n] = Lx_current
        L[1][n] = Ly_current
        L[2][n] = Lz_current

    L_T = L.transpose(1, 0, 2, 3)
    Lproduct = np.tensordot(L_T, L, axes=((1, 3), (0, 2))).swapaxes(1, 2)

    Lz = L[2]  # array of Lz operators
    H = np.tensordot(v, Lz, axes=1)
    # print(" H : ", H)

    J = np.array(J)  # convert to numpy array first
    scalars = 0.5 * J
    H += np.tensordot(scalars, Lproduct, axes=2)

    # print(" H : ", H)
    n = 2 ** nspins
    T = np.zeros((n, n))
    idx_number = 0
    for i in range(n - 1):
        for j in range(i + 1, n):
            if bin(i ^ j).count('1') == 1:
                T[i, j] = 1
    T += T.T
    return H, T


def qm_full(v: list[float], J: np.ndarray, nIntergals: int, args: argparse.Namespace) -> list:
    '''
    v, J : unit Hz 
    nIntergals : the total numbers of Intensities
    '''
    ic(v, J)
    nspins: int = len(v)
    H, T = qm_parameter(v, J)

    # ic(T)
    E, V = np.linalg.eigh(H)
    ic(E, V)
    V = V.real
    I = np.square(V.T.dot(T.dot(V)))
    # symmetry makes it possible to use only one half of the matrix for faster calculation
    I_upper = np.triu(I)
    ic(I)
    # print(" E : \n", E[:, np.newaxis]-[0, 0, 0, 0])
    E_matrix = np.abs(E[:, np.newaxis] - E)
    # print(" E_matrix : \n", E_matrix)
    E_upper = np.triu(E_matrix)
    combo = np.stack([E_upper, I_upper])
    iv = combo.reshape(2, I.shape[0] ** 2).T
    # an arbitrary cutoff where peaks below this intensity are filtered out of the solution
    peaklist = iv[iv[:, 1] >= args.cutoff]
    # peaklist.T[0] = peaklist.T[0]/args.mf
    ic(peaklist)
    from nmrsim.math import normalize_peaklist
    normalized_plist = normalize_peaklist(peaklist, nIntergals)
    print_plot(normalized_plist, dpi=1000,
               nIntergals=nIntergals, args=args, Active_range=10)
    return normalized_plist


def qm_partial(v: list[float], J: np.ndarray, idx0_nspins, nIntergals, args: argparse.Namespace) -> list:
    '''
    v, J : unit Hz 
    nIntergals : the total numbers of Intensities 
    idx0_nspins : for each spin serial numbers in nspins from 0 to n-1 
    '''
    nspins = len(v)
    H, T = qm_parameter(v, J)

    n = 2 ** nspins
    F = np.zeros((n, n), dtype=int)
    idx = int(2**(nspins-idx0_nspins-1))
    # idx = ~int(2**idx0_nspins)+1
    for i in range(n - 1):
        for j in range(i + 1, n):
            if bin(i ^ j).count('1') == 1:
                if bin((i & idx) ^ (j & idx)).count('1') == 1:
                    F[i][j] = 1
    F += F.T
    # print(" F : \n", F)
    E, V = np.linalg.eigh(H)
    V = V.real
    # print(" E : \n", E)
    # print(" I : \n", I)
    # symmetry makes it possible to use only one half of the matrix for faster calculation
    I = np.square(V.T.dot(T.dot(V)))
    # print(" I : \n",I)
    IF = np.square(V.T.dot(F.dot(V)))
    # print(" IF : \n",IF)
    # print(" I*IF : \n",I*IF)
    I_upper = np.triu(I*IF)
    # print(" I : \n", I)
    # print(" F : \n", F)

    E_matrix = np.abs(E[:, np.newaxis] - E)
    # print(" E_matrix : \n", E_matrix)

    E_upper = np.triu(E_matrix)

    combo = np.stack([E_upper, I_upper])
    iv = combo.reshape(2, I.shape[0] ** 2).T
    # print(" iv \n", iv)

    peaklist = iv[iv[:, 1] >= args.cutoff]
    # peaklist.T[0] = peaklist.T[0]/args.mf
    # print(" peaklist: \n", peaklist)
    # print(" iv[:,1]: ", iv[:, 1])
    from nmrsim.math import normalize_peaklist
    normalized_plist = normalize_peaklist(peaklist, nIntergals)
    # print(normalized_plist)
    return normalized_plist


def print_plot(inpeaklist: list, dpi: int, nIntergals: int, args: argparse.Namespace, Active_range: int, hidden=True) -> np.ndarray:
    '''
    peaklist :  Chemical Shift(Hz), Intensities
    dpi: 10000 for Hydrogen, 500 for Carbon
    nIntergals : total numbers of Intensities
    '''
    from nmrsim.math import normalize_peaklist
    peaklist: np.ndarray = np.array(inpeaklist)
    peaklist.T[0] = peaklist.T[0] / args.mf
    # print(" peaklist : \n", peaklist)
    normalized_plist: list[tuple[float, float]
                           ] = normalize_peaklist(peaklist, nIntergals)
    # print(" normalized_plist : ", normalized_plist)
    if args.start == None:
        args.start = (peaklist.T)[0].min() - Active_range * 0.1
    if args.end == None:
        args.end = (peaklist.T)[0].max() + Active_range * 0.1

    y_max_normalized = np.array(normalized_plist).max()*2

    args.start = round(args.start, 4)
    args.end = round(args.end, 4)

    lw = args.lw * 2 / 1000
    lw_points: int = int((args.end - args.start) * dpi)+1

    from nmrsim.plt import mplplot
    if hidden == False:
        if input("    Do you want to show matplotlib results?    ") in ("y", "yes"):
            x, y = mplplot(normalized_plist, w=lw, y_max=y_max_normalized, y_min=-
                           y_max_normalized*0.01, limits=(args.start, args.end), points=lw_points)
        else:
            x, y = mpl_plot(normalized_plist, w=lw, y_max=y_max_normalized, y_min=-
                            y_max_normalized*0.01, limits=(args.start, args.end), points=lw_points, hidden=True)
        print("    Plot is saved to {} !".format(args.out))
        np.savetxt(args.out, np.vstack((x, y)).T, fmt='%2.5f %12.5e')
        print("    All Done!")
        return np.vstack((x, y))
    elif hidden == True:
        x, y = mpl_plot(normalized_plist, w=lw, y_max=y_max_normalized, y_min=-
                        y_max_normalized*0.01, limits=(args.start, args.end), points=lw_points, hidden=True)
        if not args.bobyqa:
            np.savetxt(args.out, np.vstack((x, y)).T, fmt='%2.5f %12.5e')
        return np.vstack((x, y))
    else:
        print("something wrong in your print_plot method hidden setting")
        ic()
        exit(1)


def mpl_plot(peaklist, w=1, y_min=-0.01, y_max=1, points=800, limits=None, hidden=False):
    """
    Modification by Vitamin Cheng
    Only for no show picture  for hidden == True 

    """
    from nmrsim.plt import low_high, add_lorentzians
    if hidden == True:
        peaklist.sort()
        if limits:
            l_limit, r_limit = low_high(limits)
        else:
            l_limit = peaklist[0][0] - 50
            r_limit = peaklist[-1][0] + 50
        x = np.linspace(float(l_limit), float(r_limit), points)
        y = add_lorentzians(x, peaklist, w)
        return x, y
    else:
        print("Please use the nmrsim mplplot ")
        ic()
        exit(0)


def qm_base(v: list, J: np.ndarray, nIntergals, idx0_nspins, args: argparse.Namespace) -> list:
    """ QM_Base 
    v           : list[float] unit Hz
    J           : np.ndarray  unit Hz
     nIntergals : the total numbers of Intensities 
    idx0_npsins : idx numbers in AB quartet
    Returns     :
    peaklist    : list[(peak,intensity)] 
    """
    peaklist: list = []
    if len(v) > 1:
        peaklist = qm_partial(v=v, J=J, idx0_nspins=idx0_nspins,
                              nIntergals=nIntergals, args=args)
    elif len(v) == 1:
        import math
        peaklist = [(math.fabs(v[0]), 1.00000)]
    else:
        print("something wrong in your qm_Base cal.")
    return peaklist


def qm_multiplet(v: list[float], nIntergals, J: list) -> list:
    '''
    v : Hz 
    J : Hz
    nIntergals : the total numbers of Intensities 
    '''
    from nmrsim import Multiplet
    td = Multiplet(v, nIntergals, J)
    # ic(td.peaklist())
    return td.peaklist()


if __name__ == "__main__":

    x: dict = {"out": "output.dat", "start": -
               0.5, "end": 10.5, "lw": 1, "mf": 500.0, "cutoff": 0.001, "debug": False}

    v: list = [964, 2775.76, 2768.20, 928]

    J: np.ndarray = np.array([[0.0,   0.0,   0.0,   0.0],
                             [0.0,   0.0, 16.97,   0.0],
                              [0.0, 16.97,   0.0,   7.0],
                              [0.0,   0.0,   7.0,   0.0]])

    R_peak: list = qm_full(v=v, J=J, nIntergals=1,
                           args=argparse.Namespace(**x))
    print_plot(inpeaklist=R_peak, dpi=10000, nIntergals=2,
               Active_range=10, args=argparse.Namespace(**x), hidden=False)
