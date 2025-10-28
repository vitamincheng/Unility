#!/usr/bin/env python
import numpy as np
import numpy.typing as npt
from numba import jit
import argparse
from icecream import ic
from cachier import cachier

from censo_ext.anmr import Anmr
type cplex = npt.NDArray[np.complex128]


@cachier(separate_files=True)
def Pauil_matrix(nspins: int) -> tuple[npt.NDArray[np.complex128], npt.NDArray[np.complex128]]:
    """
    Create Pauli matrices for a given number of spins.

    This function generates the standard Pauli matrices (sigma_x, sigma_y, sigma_z)
    scaled by 1/2, which are fundamental operators in quantum mechanics for describing
    spin-1/2 particles.

    Args:
        nspins: Number of spins in the system

    Returns:
        Tuple containing the three Pauli matrices (sigma_x, sigma_y, sigma_z)
    """

    sigma_x: cplex = np.array([[0, 1 / 2], [1 / 2, 0]])
    sigma_y: cplex = np.array([[0, -1j / 2], [1j / 2, 0]])
    sigma_z: cplex = np.array([[1 / 2, 0], [0, -1 / 2]])
    unit: cplex = np.array([[1, 0], [0, 1]])

    L: cplex = np.empty(
        (3, nspins, 2 ** nspins, 2 ** nspins), dtype=np.complex128)
    for n in range(nspins):
        Lx_current: cplex = np.array([1])
        Ly_current: cplex = np.array([1])
        Lz_current: cplex = np.array([1])

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

    L_T: cplex = L.transpose(1, 0, 2, 3)
    Lproduct: cplex = np.tensordot(
        L_T, L, axes=((1, 3), (0, 2))).swapaxes(1, 2).astype(np.complex128)

    return L, Lproduct


@cachier(separate_files=True)
def F_matrix(nspins: int, idx0_nspins: int) -> npt.NDArray[np.uint8]:
    """
    Generate interaction matrix F for spin systems.

    This function creates a matrix that represents interactions between spins
    in a quantum system, typically used in the context of quantum Hamiltonians.

    Args:
        nspins: Number of spins in the system
        idx0_nspins: Index parameter used in bit manipulation for interaction calculation

    Returns:
        Matrix F representing spin interactions
    """

    n: int = 2 ** nspins
    F: npt.NDArray[np.uint8] = np.zeros((n, n), dtype=np.uint8)
    idx: int = int(2**(nspins-idx0_nspins-1))
    # idx = ~int(2**idx0_nspins)+1
    for i in range(n - 1):
        for j in range(i + 1, n):
            if bin((i & idx) ^ (j & idx)).count('1') == 1:
                # if bin(i ^ j).count('1') == 1 and bin((i & idx) ^ (j & idx)).count('1') == 1:
                F[i][j] = 1
    return F


@cachier(separate_files=True)
def T_matrix(nspins: int) -> npt.NDArray[np.uint8]:
    """
    Generate transition matrix T for spin systems.

    This function creates a binary matrix that indicates which states can transition
    to each other in a quantum spin system, where transitions occur when only one
    spin flips (Hamming distance of 1).

    Args:
        nspins: Number of spins in the system

    Returns:
        Matrix T representing possible spin transitions
    """

    n: int = 2 ** nspins
    T: npt.NDArray[np.uint8] = np.zeros((n, n), dtype=np.uint8)
    for i in range(n - 1):
        for j in range(i + 1, n):
            if bin(i ^ j).count('1') == 1:
                T[i, j] = 1
    T += T.T
    return T


def qm_parameter(v: list[float], J: npt.NDArray[np.float64]) -> tuple[npt.NDArray[np.complex128], npt.NDArray[np.uint8]]:
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

    L, Lproduct = Pauil_matrix(len(v))
    T: npt.NDArray[np.uint8] = T_matrix(len(v))

    Lz = L[2]  # array of Lz operators
    H: npt.NDArray[np.complex128] = np.tensordot(
        v, Lz, axes=1).astype(np.complex128)
    # ic(H)

    scalars: npt.NDArray[np.float64] = 0.5 * J
    H += np.tensordot(scalars, Lproduct, axes=2)

    return H, T


def qm_full(v: list[float], J: npt.NDArray[np.float64], args: argparse.Namespace) -> list[tuple[float, float]]:
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

    H, T = qm_parameter(v, J)

    E: npt.NDArray[np.float64]
    V: npt.NDArray[np.complex128 | np.float64]

    E, V = np.linalg.eigh(H)
    if args.verbose:
        ic(H)
        ic(T)
        ic(E, V)
        np.savetxt("Hamiltonian.out", H, fmt="%6.2f")
        np.savetxt("eigenValue.out", E.real, fmt="%6.2f")
        np.savetxt("eigenVector.out", V.real, fmt="%6.2f")
    V = V.real
    I_np: npt.NDArray[np.float64] = np.square(V.T.dot(T.dot(V)))

    # symmetry makes it possible to use only one half of the matrix for faster calculation
    I_upper: npt.NDArray[np.float64] = np.triu(I_np)
    E_matrix: npt.NDArray[np.float64] = np.abs(E[:, np.newaxis] - E)
    E_upper: npt.NDArray[np.float64] = np.triu(E_matrix)
    combo: npt.NDArray[np.float64] = np.stack([E_upper, I_upper])
    iv: npt.NDArray[np.float64] = combo.reshape(2, I_np.shape[0] ** 2).T

    # an arbitrary cutoff where peaks below this intensity are filtered out of the solution
    peaklist: npt.NDArray[np.float64] = iv[iv[:, 1] >= args.cutoff]
    if args.verbose:
        ic(I_upper)
        ic(E_matrix)
        np.savetxt("E_matrix.out", E_matrix, fmt="%7.2f")
        np.savetxt("I_matrix.out", I_np, fmt="%7.2f")
        ic(peaklist)

    freq, intensit = [x for x, y in peaklist], [y for x, y in peaklist]
    return list(zip(freq, intensit))


def qm_partial(v: list[float], J: npt.NDArray[np.float64], idx0_nspins, args: argparse.Namespace) -> list[tuple[float, float]]:
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

    H, T = qm_parameter(v, J)
    F: npt.NDArray[np.uint8] = F_matrix(nspins, idx0_nspins)

    F += F.T
    F = F*T
    E: npt.NDArray[np.float64]
    V: npt.NDArray[np.complex128 | np.float64]
    E, V = np.linalg.eigh(H)
    V = V.real
    if args.verbose:
        ic(F)
        ic(E, V)

    # symmetry makes it possible to use only one half of the matrix for faster calculation
    I_np: npt.NDArray[np.float64] = np.square(V.T.dot(T.dot(V)))
    IF: npt.NDArray[np.float64] = np.square(V.T.dot(F.dot(V)))
    I_upper: npt.NDArray[np.float64] = np.triu(I_np*IF)
    if args.verbose:
        ic(I_np)
        ic(IF)
        ic(I_np*IF)

    E_matrix: npt.NDArray[np.float64] = np.abs(E[:, np.newaxis] - E)

    E_upper: npt.NDArray[np.float64] = np.triu(E_matrix)

    combo: npt.NDArray[np.float64] = np.stack([E_upper, I_upper])
    iv: npt.NDArray[np.float64] = combo.reshape(2, I_np.shape[0] ** 2).T
    thr: np.float64 = np.max(iv[:, 1])*args.cutoff
    peaklist: npt.NDArray[np.float64] = iv[iv[:, 1] >= thr]
    if args.verbose:
        ic(E_matrix)
        ic(iv)
        ic(peaklist)
        ic(iv[:, 1])

    freq, intensit = [x for x, y in peaklist], [y for x, y in peaklist]
    return list(zip(freq, intensit))


def print_plot(inAnmr: Anmr, in_plist: list[tuple[float, float]], dpi: int,
               args: argparse.Namespace, Active_range: int) -> npt.NDArray[np.float64]:
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
    plist: npt.NDArray[np.float64] = np.array(in_plist)
    plist.T[0] = plist.T[0] / args.mf
    a, b = inAnmr.get_Anmrrc_linear()
    plist.T[0] = a*plist.T[0]+b
    Normal_plist = plist.tolist()
    if args.verbose:
        ic(plist)
        ic(Normal_plist)
    if not args.start:
        args.start = (plist.T)[0].min() - Active_range
    if not args.end:
        args.end = (plist.T)[0].max() + Active_range

    limits: tuple[float, float] = round(args.start, 4), round(args.end, 4)

    lw: float = args.lw * 2 / 1000
    lw_points: int = int((args.end - args.start) * dpi)+1

    xy_curve: tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]] = \
        mpl_plot(Normal_plist, lw=lw, limits=limits, lw_points=lw_points)
    from censo_ext.Tools.utility import save_simulation_spectra_file
    save_simulation_spectra_file(args.out, np.vstack(xy_curve).T)
    return np.vstack(xy_curve)


def mpl_plot(plist: list[tuple[float, float]], limits: tuple[float, float], lw=1.0, lw_points=200_000) \
        -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
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
    plist.sort()
    if limits:
        l_limit, r_limit = min(limits), max(limits)
    else:
        l_limit: float = plist[0][0] - 50
        r_limit: float = plist[-1][0] + 50
    x: npt.NDArray[np.float64] = np.linspace(
        float(l_limit), float(r_limit), lw_points).astype(np.float64)
    y: npt.NDArray[np.float64] = add_lorentzians(x, plist, lw)
    return x, y


def add_lorentzians(linspace: npt.NDArray[np.float64], plist: list[tuple[float, float]], lw: float) -> npt.NDArray[np.float64]:

    for freq, intensit in plist:
        try:
            result += lorentz(linspace, freq, intensit, lw)  # type: ignore # nopep8
        except NameError:
            result = lorentz(linspace, plist[0][0], plist[0][1], lw)
    return result  # type: ignore #nopep8


@jit
def lorentz(linspace: npt.NDArray[np.float64], freq: float, Intensity: float, lw: float) -> npt.NDArray[np.float64]:
    scaling_factor: float = 0.5 / lw
    return scaling_factor * Intensity * ((0.5 * lw) ** 2 / ((0.5 * lw) ** 2 + (linspace - freq) ** 2))


def qm_base(v: list[float], J: npt.NDArray[np.float64], idx0_nspins, args: argparse.Namespace) -> list[tuple[float, float]]:
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
    if args.verbose:
        ic(v, J)
    if len(v) > 1:
        plist = qm_partial(v=v, J=J, idx0_nspins=idx0_nspins, args=args)
    elif len(v) == 1:
        plist = [(np.fabs(v[0]), float(1.0))]
    else:
        print("something wrong in your qm_Base cal.")
    return plist


def qm_multiplet(v: float | int, nIntergals, J: list[tuple[float, int]], delta: list[float]) -> list[tuple[float, float]]:
    """
    Calculate multiplet spectrum 

    This function generates a multiplet spectrum for a single spin system
    with specified coupling constants and number of peaks.

    Args:
        v (float | int): Chemical shift in Hz.
        nIntergals (int): The total number of intensities to generate.
        J (list[tuple[float, int]]): List of (coupling_constant, multiplicity) tuples.

    Returns:
        list[tuple[float, float]]: Normalized peaklist with (frequency, intensity) tuples.
    """
    return Multiplet(v, nIntergals, J, delta).peaklist()


class Multiplet:

    def __init__(self, v: float, nIntergals: int, J: list[tuple[float, int]], delta: list[float], w=0.5):
        self.v: float = v
        self.nIntergals: int = nIntergals
        self.J: list[tuple[float, int]] = J
        self.delta: list[float] = delta
        self.w: float = w
        self._peaklist: list = multiplet((v, nIntergals), J, delta)

    def _refresh(self) -> None:
        self._peaklist = multiplet(
            (self.v, self.nIntergals), self.J, self.delta)

    def peaklist(self) -> list:
        self._refresh()
        return self._peaklist


def multiplet(signal: tuple[float, int], JCoups: list[tuple[float, int]], delta: list[float]) -> list[tuple[float, float]]:
    res: list = [signal]
    for idx, JCoup in enumerate(JCoups):
        for _ in range(JCoup[1]):
            res = _doublet(res, JCoup[0], delta[idx])
    return reduce_peaks(res)


def reduce_peaks(plist_: list[tuple[float, float]], tolerance=0.02) -> list[tuple[float, float]]:
    res: list[tuple[float, float]] = []
    work: list[tuple[float, float]] = []  # an accumulator of peaks to be added
    plist: list[tuple[float, float]] = sorted(plist_)
    for peak in plist:
        if not work:
            work.append(peak)
            continue
        if peak[0] - work[-1][0] <= tolerance:
            work.append(peak)
            continue
        else:
            res.append(add_peaks(work))
            work = [peak]
    if work:  # process any remaining work after for loop
        res.append(add_peaks(work))
    return res


def add_peaks(plist: list[tuple[float, float]]) -> tuple[float, float]:
    v_total = 0
    i_total = 0
    for v, i in plist:
        v_total += v
        i_total += i
    return v_total / len(plist), i_total


def _doublet(plist: list[tuple[float, int]], JCoups, delta) -> list[tuple[float, float]]:
    # see http://www.ebyte.it/library/docs/kts/KTS_isoAB_Geometry.html
    # if c is positive, peaks must be the left of doublet is more low and the right is more high
    # if c is negative, peaks must be the left of doublet is more high and the right is more low
    #
    if (JCoups+delta) == 0:
        _k = 0
    else:
        _k = JCoups / (JCoups+delta)

    k_small: float = 1 - _k
    k_large: float = 1 + _k

    res: list = []
    for v, intensit in plist:
        # the left of doublet if J is positive
        res.append((v + JCoups / 2, intensit / 2 * k_small))
        # the right of doublet if J is positive
        res.append((v - JCoups / 2, intensit / 2 * k_large))
    return res
