#!/usr/bin/env python
import numpy as np
import numpy.typing as npt
import argparse
from icecream import ic


def qm_parameter(v: list[float], J: npt.NDArray[np.float64]) -> tuple[npt.NDArray[np.complex128], npt.NDArray[np.float64]]:
    """
    Calculate the Hamiltonian and transition matrix for a spin system.

    This function constructs the angular momentum operators (Lx, Ly, Lz) for each spin
    and builds the total Hamiltonian H from the Zeeman terms (v) and dipolar coupling terms (J).

    Args:
        v (list[float]): List of resonance frequencies in Hz for each spin.
        J (npt.NDArray[np.float64]): Dipolar coupling matrix (Hz) with shape (nspins, nspins).

    Returns:
        tuple[npt.NDArray[np.complex128], npt.NDArray[np.float64]]: 
        - H: The total Hamiltonian matrix (complex128)
        - T: Transition matrix for intensity calculations (float64)
    """
    sigma_x: npt.NDArray[np.complex128] = np.array([[0, 1 / 2], [1 / 2, 0]])
    sigma_y: npt.NDArray[np.complex128] = np.array([[0, -1j / 2], [1j / 2, 0]])
    sigma_z: npt.NDArray[np.complex128] = np.array([[1 / 2, 0], [0, -1 / 2]])
    unit: npt.NDArray[np.complex128] = np.array([[1, 0], [0, 1]])

    nspins: int = len(v)
    L: npt.NDArray[np.complex128] = np.empty(
        (3, nspins, 2 ** nspins, 2 ** nspins), dtype=np.complex128)
    for n in range(nspins):
        Lx_current: npt.NDArray[np.complex128] = np.array([1])
        Ly_current: npt.NDArray[np.complex128] = np.array([1])
        Lz_current: npt.NDArray[np.complex128] = np.array([1])

        for k in range(nspins):
            if k == n:
                Lx_current = np.kron(Lx_current, sigma_x).astype(np.complex128)
                Ly_current = np.kron(Ly_current, sigma_y).astype(np.complex128)
                Lz_current = np.kron(Lz_current, sigma_z).astype(np.complex128)
            else:
                Lx_current = np.kron(Lx_current, unit).astype(np.complex128)
                Ly_current = np.kron(Ly_current, unit).astype(np.complex128)
                Lz_current = np.kron(Lz_current, unit).astype(np.complex128)

        L[0][n] = Lx_current
        L[1][n] = Ly_current
        L[2][n] = Lz_current

    L_T: npt.NDArray[np.complex128] = L.transpose(1, 0, 2, 3)
    Lproduct: npt.NDArray[np.complex128] = np.tensordot(
        L_T, L, axes=((1, 3), (0, 2))).swapaxes(1, 2).astype(np.complex128)

    Lz = L[2]  # array of Lz operators
    H: npt.NDArray[np.complex128] = np.tensordot(
        v, Lz, axes=1).astype(np.complex128)
    # ic(H)

    # J = np.array(J)  # convert to numpy array first
    scalars: npt.NDArray[np.float64] = 0.5 * J
    H += np.tensordot(scalars, Lproduct, axes=2)

    # ic(H)
    n: int = 2 ** nspins
    T: npt.NDArray[np.float64] = np.zeros((n, n), dtype=np.float64)
    for i in range(n - 1):
        for j in range(i + 1, n):
            if bin(i ^ j).count('1') == 1:
                T[i, j] = 1
    T += T.T
    return H, T


def qm_full(v: list[float], J: npt.NDArray[np.float64], nIntergals: int, args: argparse.Namespace) -> list[tuple[float, float]]:
    """
    Calculate full spin system spectrum using quantum mechanical approach.

    This function computes the complete energy eigenvalues and eigenvectors for a 
    spin system, calculates intensities based on transition matrix, and normalizes
    the resulting peaklist.

    Args:
        v (list[float]): List of resonance frequencies in Hz for each spin.
        J (npt.NDArray[np.float64]): Dipolar coupling matrix (Hz) with shape (nspins, nspins).
        nIntergals (int): The total number of intensities to generate.
        args (argparse.Namespace): Command line arguments containing plotting parameters.

    Returns:
        list[tuple[float, float]]: Normalized peaklist with (frequency, intensity) tuples.
    """
    nspins: int = len(v)
    if J.shape != (nspins, nspins):
        raise ValueError("Your JCoup is Error")

    H: npt.NDArray[np.complex128]
    T: npt.NDArray[np.float64]
    H, T = qm_parameter(v, J)

    # ic(H)
    E: npt.NDArray[np.float64]
    V: npt.NDArray[np.complex128 | np.float64]

    E, V = np.linalg.eigh(H)
    # ic(E, V)
    V = V.real
    I_np: npt.NDArray[np.float64] = np.square(V.T.dot(T.dot(V)))

    # symmetry makes it possible to use only one half of the matrix for faster calculation
    I_upper: npt.NDArray[np.float64] = np.triu(I_np)
    # ic(I)
    E_matrix: npt.NDArray[np.float64] = np.abs(E[:, np.newaxis] - E)
    # ic(E_matrix)
    E_upper: npt.NDArray[np.float64] = np.triu(E_matrix)
    combo: npt.NDArray[np.float64] = np.stack([E_upper, I_upper])
    iv: npt.NDArray[np.float64] = combo.reshape(2, I_np.shape[0] ** 2).T

    # an arbitrary cutoff where peaks below this intensity are filtered out of the solution
    peaklist: npt.NDArray[np.float64] = iv[iv[:, 1] >= args.cutoff]
    # ic(peaklist)
    from nmrsim.math import normalize_peaklist
    Normal_plist: list[tuple[float, float]
                       ] = normalize_peaklist(peaklist, nIntergals)
    # print_plot(Normal_plist, dpi=1000,
    #           nIntergals=nIntergals, args=args, Active_range=10)
    return Normal_plist


def qm_partial(v: list[float], J: npt.NDArray[np.float64], idx0_nspins, nIntergals, args: argparse.Namespace) -> list[tuple[float, float]]:
    """
    Calculate partial spin system spectrum for a specific spin.

    This function computes the spectrum contribution from a single spin (idx0_nspins) 
    by restricting transitions to only those involving that spin.

    Args:
        v (list[float]): List of resonance frequencies in Hz for each spin.
        J (npt.NDArray[np.float64]): Dipolar coupling matrix (Hz) with shape (nspins, nspins).
        idx0_nspins (int): Index of the spin to calculate spectrum for (0-based).
        nIntergals (int): The total number of intensities to generate.
        args (argparse.Namespace): Command line arguments containing plotting parameters.

    Returns:
        list[tuple[float, float]]: Normalized peaklist with (frequency, intensity) tuples.
    """
    nspins: int = len(v)
    if J.shape != (nspins, nspins):
        raise ValueError("Your JCoup is Error")
    if idx0_nspins >= nspins:
        raise ValueError("Your idx0_nspins is Error")

    H: npt.NDArray[np.complex128]
    T: npt.NDArray[np.float64]
    H, T = qm_parameter(v, J)

    n: int = 2 ** nspins
    F: npt.NDArray[np.float64] = np.zeros((n, n), dtype=np.float64)
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
    E: npt.NDArray[np.float64]
    V: npt.NDArray[np.complex128 | np.float64]
    E, V = np.linalg.eigh(H)
    V = V.real
    # ic(E,V)

    # symmetry makes it possible to use only one half of the matrix for faster calculation
    I_np: npt.NDArray[np.float64] = np.square(V.T.dot(T.dot(V)))
    # ic(I)
    IF: npt.NDArray[np.float64] = np.square(V.T.dot(F.dot(V)))
    I_upper: npt.NDArray[np.float64] = np.triu(I_np*IF)
    # ic(IF)
    # ic(I*IF)

    E_matrix: npt.NDArray[np.float64] = np.abs(E[:, np.newaxis] - E)
    # ic(E_matrix)

    E_upper: npt.NDArray[np.float64] = np.triu(E_matrix)

    combo: npt.NDArray[np.float64] = np.stack([E_upper, I_upper])
    iv: npt.NDArray[np.float64] = combo.reshape(2, I_np.shape[0] ** 2).T
    # ic(iv)
    thr: np.float64 = np.max(iv[:, 1])*args.cutoff
    peaklist: npt.NDArray[np.float64] = iv[iv[:, 1] >= thr]
    # ic(peaklist)
    # ic(iv[:, 1])
    from nmrsim.math import normalize_peaklist
    Normal_plist: list[tuple[float, float]] = np.array(normalize_peaklist(
        peaklist, nIntergals)).tolist()
    return Normal_plist


def print_plot(in_plist: list[tuple[float, float]], dpi: int, nIntergals: int, args: argparse.Namespace, Active_range: int) -> npt.NDArray:
    """
    Generate and save a plot of the NMR spectrum.

    This function creates a matplotlib plot of the peaklist with specified parameters
    and saves it to a file. It also outputs the data to a text file.

    Args:
        in_plist (list[tuple[float, float]]): List of (frequency, intensity) tuples.
        dpi (int): Plot resolution (10000 for Hydrogen, 500 for Carbon).
        nIntergals (int): Total number of intensities to generate.
        args (argparse.Namespace): Command line arguments containing plotting parameters.
        Active_range (int): Range around spectrum to display.

    Returns:
        npt.NDArray: Array containing x and y coordinates of the plot data.
    """
    from nmrsim.math import normalize_peaklist
    plist: npt.NDArray[np.float64] = np.array(in_plist)
    plist.T[0] = plist.T[0] / args.mf
    # ic(peaklist)
    Normal_plist: list[tuple[float, float]
                       ] = normalize_peaklist(plist, nIntergals)
    # ic(normalized_plist)
    if not args.start:
        args.start = (plist.T)[0].min() - Active_range * 0.1
    if not args.end:
        args.end = (plist.T)[0].max() + Active_range * 0.1

    Normal_y_max: float = np.array(Normal_plist).max()*2

    args.start = round(args.start, 4)
    args.end = round(args.end, 4)

    lw: float = args.lw * 2 / 1000
    lw_points: int = int((args.end - args.start) * dpi)+1

    x: npt.NDArray[np.float64]
    y: npt.NDArray[np.float64]
    x, y = mpl_plot(Normal_plist, w=lw, y_max=Normal_y_max, y_min=-
                    Normal_y_max*0.01, limits=(args.start, args.end), points=lw_points)
    if not args.bobyqa:
        np.savetxt(args.out, np.vstack((x, y)).T, fmt='%2.5f %12.5e')
        print(f" the spectra is saved to : {args.out}")
    return np.vstack((x, y))


def mpl_plot(plist, w=1.0, y_min=-0.01, y_max=1.0, points=800, limits=None) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """
    Generate a plot using lorentzian lineshape for NMR spectrum.

    This function creates an NMR spectrum plot using lorentzian line shapes.
    It's designed to be used internally by print_plot when hidden mode is enabled.

    Args:
        plist (list[tuple[float, float]]): List of (frequency, intensity) tuples.
        w (float): Lorentzian width parameter (default: 1.0).
        y_min (float): Minimum y-axis value for plot (default: -0.01).
        y_max (float): Maximum y-axis value for plot (default: 1.0).
        points (int): Number of points to generate for the curve (default: 800).
        limits (tuple, optional): x-axis limits as (min, max) tuple.

    Returns:
        tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]: 
        - x: Array of x-coordinates
        - y: Array of y-coordinates
    """
    from nmrsim.plt import low_high, add_lorentzians
    plist.sort()
    if limits:
        l_limit, r_limit = low_high(limits)
    else:
        l_limit = plist[0][0] - 50
        r_limit = plist[-1][0] + 50
    x: npt.NDArray[np.float64] = np.linspace(
        float(l_limit), float(r_limit), points).astype(np.float64)
    y: npt.NDArray[np.float64] = add_lorentzians(x, plist, w)
    return x, y


def qm_base(v: list[float], J: npt.NDArray[np.float64], nIntergals, idx0_nspins, args: argparse.Namespace) -> list[tuple[float, float]]:
    """
    Base quantum mechanical calculation function for spin systems.

    This function serves as the main interface for quantum mechanical calculations,
    handling both single spin and multi-spin cases appropriately.

    Args:
        v (list[float]): List of resonance frequencies in Hz for each spin.
        J (npt.NDArray[np.float64]): Dipolar coupling matrix (Hz) with shape (nspins, nspins).
        nIntergals (int): The total number of intensities to generate.
        idx0_nspins (int): Index of the spin to calculate spectrum for (0-based), 
                          used in partial calculations.
        args (argparse.Namespace): Command line arguments containing plotting parameters.

    Returns:
        list[tuple[float, float]]: Normalized peaklist with (frequency, intensity) tuples.
    """
    plist: list[tuple[float, float]] = []
    # ic(v, J)
    if len(v) > 1:
        plist = qm_partial(v=v, J=J, idx0_nspins=idx0_nspins,
                           nIntergals=nIntergals, args=args)
    elif len(v) == 1:
        import math
        plist = [(math.fabs(v[0]), float(1.00000))]
    else:
        print("something wrong in your qm_Base cal.")
    return plist


def qm_multiplet(v: float | int, nIntergals, J: list[tuple[float, int]]) -> list[tuple[float, float]]:
    """
    Calculate multiplet spectrum using nmrsim library.

    This function generates a multiplet spectrum for a single spin system 
    with specified coupling constants and number of peaks.

    Args:
        v (float | int): Chemical shift in Hz.
        nIntergals (int): The total number of intensities to generate.
        J (list[tuple[float, int]]): List of (coupling_constant, multiplicity) tuples.

    Returns:
        list[tuple[float, float]]: Normalized peaklist with (frequency, intensity) tuples.
    """
    from nmrsim import Multiplet
    td = Multiplet(v, nIntergals, J)
    return td.peaklist()


if __name__ == "__main__":

    x: dict = {"out": "output.dat", "start": -
               0.5, "end": 10.5, "lw": 0.1, "mf": 500.0, "cutoff": 0.001, "debug": False, "bobyqa": False}
    # v: positive or negative of the frequency is the same of spectra
    # J: only one AB quartet, positive or negative of the J Coupling constant the spectra is the same

    # v: list = [1100, 1200, 1900, 2500]
    # J: npt.NDArray[np.float64] = np.array([[0.0,   -16.0,   0.0,   4.0],
    #                         [-16.0,   0.0, 2.0,   4.0],
    #                          [0.0, 2.0,   0.0,   0.0],
    #                          [4.0,   4.0,   0.0,   0.0]])

    v: list[float] = [480, 645, 645, 645, 645, 480]
    J: npt.NDArray[np.float64] = np.array([[0.00000,      7.12744,      7.12267,     -0.22011,   -0.21844,     -0.02230],
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
    ic(len(R_peak))
    print_plot(in_plist=R_peak, dpi=10000, nIntergals=2,
               Active_range=10, args=argparse.Namespace(**x))

    R_peaks: list = []
    for idx in range(len(v)):

        R_peaks += qm_partial(v=v, J=J, idx0_nspins=idx, nIntergals=1,
                              args=argparse.Namespace(**x))

    ic(R_peaks)
    ic(len(R_peaks))
    print_plot(in_plist=R_peaks, dpi=10000, nIntergals=2,
               Active_range=10, args=argparse.Namespace(**x))
