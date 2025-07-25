#!/usr/bin/env python3
import numpy as np
import os
import argparse
from icecream import ic


def qm_parameter(v: list[float], J: np.ndarray):

    sigma_x: np.ndarray = np.array([[0, 1 / 2], [1 / 2, 0]])
    sigma_y: np.ndarray = np.array([[0, -1j / 2], [1j / 2, 0]])
    sigma_z: np.ndarray = np.array([[1 / 2, 0], [0, -1 / 2]])
    unit: np.ndarray = np.array([[1, 0], [0, 1]])

    nspins: int = len(v)
    L: np.ndarray = np.empty(
        (3, nspins, 2 ** nspins, 2 ** nspins), dtype=np.complex128)
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
    H: np.ndarray = np.tensordot(v, Lz, axes=1)
    # ic(H)

    J = np.array(J)  # convert to numpy array first
    scalars = 0.5 * J
    H += np.tensordot(scalars, Lproduct, axes=2)

    # ic(H)
    n: int = 2 ** nspins
    T: np.ndarray = np.zeros((n, n))
    for i in range(n - 1):
        for j in range(i + 1, n):
            if bin(i ^ j).count('1') == 1:
                T[i, j] = 1
    T += T.T
    # ic(T)
    return H, T


def qm_full(v: list[float], J: np.ndarray, nIntergals: int, args: argparse.Namespace) -> list:
    '''
    v, J : unit Hz 
    nIntergals : the total numbers of Intensities
    '''
    nspins: int = len(v)
    if J.shape != (nspins, nspins):
        raise ValueError("Your JCoupl is Error")
        os._exit(0)

    H, T = qm_parameter(v, J)

    # ic(H)
    E, V = np.linalg.eigh(H)
    # ic(E, V)
    V = V.real
    I = np.square(V.T.dot(T.dot(V)))

    # symmetry makes it possible to use only one half of the matrix for faster calculation
    I_upper = np.triu(I)
    # ic(I)
    E_matrix = np.abs(E[:, np.newaxis] - E)
    # ic(E_matrix)
    E_upper = np.triu(E_matrix)
    combo = np.stack([E_upper, I_upper])
    iv = combo.reshape(2, I.shape[0] ** 2).T

    # an arbitrary cutoff where peaks below this intensity are filtered out of the solution
    peaklist = iv[iv[:, 1] >= args.cutoff]
    # ic(peaklist)
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
    nspins: int = len(v)
    if J.shape != (nspins, nspins):
        raise ValueError("Your JCoupl is Error")
        os._exit(0)
    if idx0_nspins >= nspins:
        raise ValueError("Your idx0_nspins is Error")
        os._exit(0)

    H, T = qm_parameter(v, J)

    n: int = 2 ** nspins
    F: np.ndarray = np.zeros((n, n), dtype=np.float64)
    idx: int = int(2**(nspins-idx0_nspins-1))
    # idx = ~int(2**idx0_nspins)+1
    for i in range(n - 1):
        for j in range(i + 1, n):
            if bin((i & idx) ^ (j & idx)).count('1') == 1:
                # if bin(i ^ j).count('1') == 1 and bin((i & idx) ^ (j & idx)).count('1') == 1:
                F[i][j] = 1
    F += F.T
    F = F*T
    # ic(F)

    E, V = np.linalg.eigh(H)
    V = V.real
    # ic(E,V)

    # symmetry makes it possible to use only one half of the matrix for faster calculation
    I: np.ndarray = np.square(V.T.dot(T.dot(V)))
    # ic(I)
    IF: np.ndarray = np.square(V.T.dot(F.dot(V)))
    I_upper: np.ndarray = np.triu(I*IF)
    # ic(IF)
    # ic(I*IF)

    E_matrix: np.ndarray = np.abs(E[:, np.newaxis] - E)
    # ic(E_matrix)

    E_upper: np.ndarray = np.triu(E_matrix)

    combo: np.ndarray = np.stack([E_upper, I_upper])
    iv: np.ndarray = combo.reshape(2, I.shape[0] ** 2).T
    # ic(iv)
    thr = np.max(iv[:, 1])*args.cutoff
    peaklist: np.ndarray = iv[iv[:, 1] >= thr]
    # ic(peaklist)
    # ic(iv[:, 1])
    from nmrsim.math import normalize_peaklist
    normalized_plist = normalize_peaklist(peaklist, nIntergals)
    # ic(normalized_plist)
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
    # ic(peaklist)
    normalized_plist: list[tuple[float, float]
                           ] = normalize_peaklist(peaklist, nIntergals)
    # ic(normalized_plist)
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
    # ic(v, J)
    if len(v) > 1:
        peaklist = qm_partial(v=v, J=J, idx0_nspins=idx0_nspins,
                              nIntergals=nIntergals, args=args)
    elif len(v) == 1:
        import math
        peaklist = [(np.fabs(v[0]), np.float64(1.00000))]
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
               0.5, "end": 10.5, "lw": 0.1, "mf": 500.0, "cutoff": 0.001, "debug": False, "bobyqa": False}
    # v: positive or negative of the frequency is the same of spectra
    # J: only one AB quartet, positive or negative of the J Coupling constant the spectra is the same

    # v: list = [1100, 1200, 1900, 2500]
    # J: np.ndarray = np.array([[0.0,   -16.0,   0.0,   4.0],
    #                         [-16.0,   0.0, 2.0,   4.0],
    #                          [0.0, 2.0,   0.0,   0.0],
    #                          [4.0,   4.0,   0.0,   0.0]])

    v: list = [480, 645, 645, 645, 645, 480]
    J: np.ndarray = np.array([[0.00000,      7.12744,      7.12267,     -0.22011,   -0.21844,     -0.02230],
                              [7.12744,      0.00000,    -13.32467,
                                  6.11500,    6.85267,     -0.21511],
                              [7.12267,    -13.32467,      0.00000,
                                  6.81333,    6.12033,     -0.21878],
                              [-0.22011,      6.11500,      6.81333,
                                  0.00000, -13.32433,      7.12589],
                              [-0.21844,      6.85267,      6.12033,    -
                                  13.32433,   0.00000,      7.12811],
                              [-0.02230,     -0.21511,     -0.21878,      7.12589,   7.12811,      0.00000]])

    ic(v)
    ic(J)
    R_peak: list = qm_full(v=v, J=J, nIntergals=1,
                           args=argparse.Namespace(**x))
    ic(R_peak)
    ic(len(R_peak))
    print_plot(inpeaklist=R_peak, dpi=10000, nIntergals=2,
               Active_range=10, args=argparse.Namespace(**x), hidden=False)

    R_peaks: list = []
    for idx in range(len(v)):

        # ic(v[idx])
        R_peaks += qm_partial(v=v, J=J, idx0_nspins=idx, nIntergals=1,
                              args=argparse.Namespace(**x))

    ic(R_peaks)
    ic(len(R_peaks))
    print_plot(inpeaklist=R_peaks, dpi=10000, nIntergals=2,
               Active_range=10, args=argparse.Namespace(**x), hidden=False)
