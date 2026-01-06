#!/usr/bin/env python
from censo_ext.Tools.Pauli import PauliComposer
from scipy.sparse import csr_matrix
import numpy as np
import numpy.typing as npt
from numba import njit
import argparse
from icecream import ic
from cachier import cachier
from censo_ext.anmr import Anmr
type cplex = npt.NDArray[np.complex64]
type np_uint = npt.NDArray[np.uint8]
type np_float = npt.NDArray[np.float64]


@cachier(separate_files=True)
def F_matrix(nspins: int, idx0_nspins: int) -> cplex:
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
    F: np_uint = np.zeros((n, n), dtype=np.uint8)
    idx: int = int(2**(nspins-idx0_nspins-1))
    # idx = ~int(2**idx0_nspins)+1
    for i in range(n - 1):
        for j in range(i + 1, n):
            if bin((i & idx) ^ (j & idx)).count('1') == 1:
                # if bin(i ^ j).count('1') == 1 and bin((i & idx) ^ (j & idx)).count('1') == 1:
                F[i][j] = 1
    return F.astype(np.complex64)


@cachier(separate_files=True)
def T_matrix(nspins: int) -> cplex:
    total: csr_matrix = Pauli_matrix("X", nspins, 1)
    for x in range(2, nspins + 1):
        total += Pauli_matrix("X", nspins, x)
    return (total.toarray()*2).astype(np.complex64)


@cachier(separate_files=True)
def Pauli_matrix(axis: str, nspins: int, idx1: int) -> csr_matrix:
    if not (axis == "X" or "Y" or "Z"):
        print(f"  {axis}, Something wrong in your Pauli matrix")
        print("  Close and Exit the program !!!")
        exit(0)

    if nspins < 1 or idx1 < 1 or nspins < idx1:
        print("  Something wrong in your npsins")
        print("  Close and Exit the program !!!")
        exit(0)

    sym: str = ("I"*(idx1-1) + str(axis)*1 + "I"*(nspins-idx1))
    return PauliComposer(sym).to_sparse()*0.5


def Hami_Zeeman(freq: list[float]) -> csr_matrix:
    nspins: int = len(freq)
    res: csr_matrix = Pauli_matrix("Z", nspins, 1)*freq[0]
    for idx1 in range(2, nspins + 1):
        res += Pauli_matrix("Z", nspins, idx1) * freq[idx1-1]
    return res


def Hami_JCoups(JCoups: npt.NDArray) -> csr_matrix:
    nspins: int = len(JCoups[0])
    res: csr_matrix = Lproduct_ij(nspins, 1, 1)*JCoups[0][0]
    for idx1_i in range(1, nspins+1):
        for idx1_j in range(idx1_i, nspins+1):
            # ic(J[i-1][j-1])
            res += Lproduct_ij(nspins, idx1_i, idx1_j) * \
                JCoups[idx1_i-1][idx1_j-1]
    return res


@cachier(separate_files=True)
def Lproduct_ij(nspins: int, idx1_i: int, idx1_j: int) -> csr_matrix:
    res: csr_matrix = Pauli_matrix(
        "X", nspins, idx1_i)*Pauli_matrix("X", nspins, idx1_j)
    res += Pauli_matrix("Y", nspins, idx1_i)*Pauli_matrix("Y", nspins, idx1_j)
    res += Pauli_matrix("Z", nspins, idx1_i)*Pauli_matrix("Z", nspins, idx1_j)
    return res


def Hamiltonian(freq: list[float], JCoups: np_float) -> csr_matrix:
    Hami: csr_matrix = Hami_Zeeman(freq)
    Hami += Hami_JCoups(JCoups)
    return Hami


def qm_full(freq: list[float], JCoups: np_float, _cutoff: float, _verbose: bool) -> list[tuple[float, float]]:
    """
    Calculate full spin system spectrum using quantum mechanical approach.

    This function computes the complete energy eigenvalues and eigenvectors for a
    spin system, calculates intensities based on transition matrix, and normalizes
    the resulting peaklist.

    Args:
        freq (list[float]): List of resonance frequencies in Hz for each spin.
        JCoups (npt.NDArray[np.float64]): Dipolar coupling matrix (Hz) with shape (nspins, nspins).
        nIntergals (int): The total number of intensities to generate.
        args (argparse.Namespace): Command line arguments containing plotting parameters.

    Returns:
        list[tuple[float, float]]: Normalized peaklist with (frequency, intensity) tuples.
    """
    nspins: int = len(freq)
    if JCoups.shape != (nspins, nspins):
        raise ValueError("Your JCoups is Error")

    H: csr_matrix = Hamiltonian(freq, JCoups)
    Energy: np_float
    Vector: cplex | np_float
    Energy, Vector = np.linalg.eigh(H.toarray())

    if _verbose:
        ic(H.toarray())
        ic(Energy, Vector)
        np.savetxt("Hamiltonian.out", H.toarray(), fmt="%6.2f")
        np.savetxt("eigenValue.out", Energy.real, fmt="%6.2f")
        np.savetxt("eigenVector.out", Vector.real, fmt="%6.2f")
    Vector = Vector.real  # type: ignore

    T: cplex = T_matrix(nspins)
    I_np: np_float = np.square(Vector.T.dot(T.dot(Vector)))

    # symmetry makes it possible to use only one half of the matrix for faster calculation
    I_upper: np_float = np.triu(I_np)
    E_matrix: np_float = np.abs(Energy[:, np.newaxis] - Energy)
    E_upper: np_float = np.triu(E_matrix)
    combo: np_float = np.stack([E_upper, I_upper])
    iv: np_float = combo.reshape(2, I_np.shape[0] ** 2).T

    # an arbitrary cutoff where peaks below this intensity are filtered out of the solution
    peaklist: np_float = iv[iv[:, 1] >= _cutoff]
    if _verbose:
        ic(I_upper)
        ic(E_matrix)
        np.savetxt("E_matrix.out", E_matrix, fmt="%7.2f")
        np.savetxt("I_matrix.out", I_np, fmt="%7.2f")
        ic(peaklist)

    freq, intensit = [x for x, y in peaklist], [y for x, y in peaklist]
    return list(zip(freq, intensit))


def qm_partial(freq: list[float], JCoups: np_float, idx0_nspins: int, _cutoff: float, _verbose: bool) -> list[tuple[float, float]]:
    """
    Calculate partial spin system spectrum for a specific spin.

    This function computes the spectrum contribution from a single spin (idx0_nspins)
    by restricting transitions to only those involving that spin.

    Args:
        freq (list[float]): List of resonance frequencies in Hz for each spin.
        JCoups (npt.NDArray[np.float64]): Dipolar coupling matrix (Hz) with shape (nspins, nspins).
        idx0_nspins (int): Index of the spin to calculate spectrum for (0-based).
        nIntergals (int): The total number of intensities to generate.
        args (argparse.Namespace): Command line arguments containing plotting parameters.

    Returns:
        list[tuple[float, float]]: Normalized peaklist with (frequency, intensity) tuples.
    """
    nspins: int = len(freq)
    if JCoups.shape != (nspins, nspins):
        raise ValueError("Your JCoups is Error")
    if idx0_nspins >= nspins:
        raise ValueError("Your idx0_nspins is Error")

    H: csr_matrix = Hamiltonian(freq, JCoups)

    Energy: np_float
    Vector: cplex | np_float

    Energy, Vector = np.linalg.eigh(H.toarray())

    Vector = Vector.real  # type: ignore
    if _verbose:
        ic(H.toarray())
        ic(Energy, Vector)
    F: cplex = F_matrix(nspins, idx0_nspins)
    T: cplex = T_matrix(nspins)
    F += F.T
    F = F*T

    # symmetry makes it possible to use only one half of the matrix for faster calculation
    I_np: np_float = np.square(Vector.T.dot(T.dot(Vector)))
    IF: np_float = np.square(Vector.T.dot(F.dot(Vector)))
    I_upper: np_float = np.triu(I_np*IF)
    if _verbose:
        ic(I_np)
        ic(IF)
        ic(I_np*IF)

    E_matrix: np_float = np.abs(Energy[:, np.newaxis] - Energy)

    E_upper: np_float = np.triu(E_matrix)

    combo: np_float = np.stack([E_upper, I_upper])
    iv: np_float = combo.reshape(2, I_np.shape[0] ** 2).T
    thr = np.max(iv[:, 1]) * _cutoff
    peaklist: np_float = iv[iv[:, 1] >= thr]
    if _verbose:
        ic(E_matrix)
        ic(iv)
        ic(peaklist)
        ic(iv[:, 1])

    freq, intensit = [x for x, y in peaklist], [y for x, y in peaklist]
    return list(zip(freq, intensit))


def print_plot(inAnmr: Anmr, in_plist: list[tuple[float, float]], dpi: int,
               args: argparse.Namespace, Active_range: int) -> np_float:
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
    plist: np_float = np.array(in_plist)
    plist.T[0] = plist.T[0] / args.mf
    a, b = inAnmr.get_Anmrrc_linear()
    plist.T[0] = a * plist.T[0] + b
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

    xy_curve: tuple[np_float, np_float] = \
        mpl_plot(Normal_plist, lw=lw, limits=limits, lw_points=lw_points)
    from censo_ext.Tools.utility import save_simulation_spectra_file
    save_simulation_spectra_file(args.out, np.vstack(xy_curve).T)
    return np.vstack(xy_curve)


def mpl_plot(plist: list[tuple[float, float]], limits: tuple[float, float], lw: float = 1.0, lw_points: int = 200_000) \
        -> tuple[np_float, np_float]:
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
    x: np_float = np.linspace(
        float(l_limit), float(r_limit), lw_points).astype(np.float64)
    y: np_float = add_lorentzians(x, plist, lw)
    return x, y


def add_lorentzians(linspace: np_float, plist: list[tuple[float, float]], lw: float) -> np_float:

    for freq, intensit in plist:
        try:
            result += lorentz(linspace, freq, intensit, lw)  # type: ignore # nopep8
        except NameError:
            result = lorentz(linspace, plist[0][0], plist[0][1], lw)
    return result  # type: ignore #nopep8


@njit
def lorentz(linspace: np_float, freq: float, Intensity: float, lw: float) -> np_float:
    scaling_factor: float = 0.5 / lw
    return scaling_factor * Intensity * ((0.5 * lw) ** 2 / ((0.5 * lw) ** 2 + (linspace - freq) ** 2))


def qm_base(freq: list[float], JCoups: np_float, idx0_nspins: int, _cutoff: float, _verbose: bool) -> list[tuple[float, float]]:
    """
    Base quantum mechanical calculation function for spin systems.

    This function serves as the main interface for quantum mechanical calculations,
    handling both single spin and multi-spin cases appropriately.

    Args:
        freq (list[float]): List of resonance frequencies in Hz for each spin.
        JCoups (npt.NDArray[np.float64]): Dipolar coupling matrix (Hz) with shape (nspins, nspins).
        nIntergals (int): The total number of intensities to generate.
        idx0_nspins (int): Index of the spin to calculate spectrum for (0-based),
                          used in partial calculations.
        args (argparse.Namespace): Command line arguments containing plotting parameters.

    Returns:
        list[tuple[float, float]]: Normalized peaklist with (frequency, intensity) tuples.
    """
    plist: list[tuple[float, float]] = []
    if _verbose:
        ic(freq, JCoups)
    if len(freq) > 1:
        plist = qm_partial(freq=freq, JCoups=JCoups, idx0_nspins=idx0_nspins,
                           _cutoff=_cutoff, _verbose=_verbose)
    elif len(freq) == 1:
        plist = [(np.fabs(freq[0]), float(1.0))]
    else:
        print("something wrong in your qm_Base cal.")
    return plist


def qm_multiplet(freq: float | int, nIntergals: int, JCoups: list[tuple[float, int]], delta: list[float], _verbose: bool) -> list[tuple[float, float]]:
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
    if _verbose:
        ic(freq, nIntergals, JCoups, delta)
    return Multiplet(freq, nIntergals, JCoups, delta, _verbose).peaklist()


class Multiplet:

    def __init__(self, freq: float, nIntergals: int, JCoups: list[tuple[float, int]], delta: list[float], _verbose: bool, w: float = 0.5) -> None:
        self.v: float = freq
        self.nIntergals: int = nIntergals
        self.J: list[tuple[float, int]] = JCoups
        self.delta: list[float] = delta
        self.w: float = w
        self._peaklist: list = multiplet((freq, nIntergals), JCoups, delta)
        self._verbose: bool = _verbose

    def _refresh(self) -> None:
        self._peaklist = multiplet(
            (self.v, self.nIntergals), self.J, self.delta)

    def peaklist(self) -> list:
        self._refresh()
        if self._verbose:
            ic(self._peaklist)
        return self._peaklist


def multiplet(signal: tuple[float, int], JCoups: list[tuple[float, int]], delta: list[float]) -> list[tuple[float, float]]:
    res: list = [signal]
    for idx, JCoup in enumerate(JCoups):
        for _ in range(JCoup[1]):
            res = _doublet(res, JCoup[0], delta[idx])
    return reduce_peaks(res)


def reduce_peaks(plist_: list[tuple[float, float]], tolerance: float = 0.02) -> list[tuple[float, float]]:
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
    freq_total = 0
    intensit_total = 0
    for freq, intensit in plist:
        freq_total += freq
        intensit_total += intensit
    return freq_total / len(plist), intensit_total


def _doublet(plist: list[tuple[float, int]], JCoups: float, delta: float) -> list[tuple[float, float]]:
    # see http://www.ebyte.it/library/docs/kts/KTS_isoAB_Geometry.html
    # if c is positive, peaks must be the left of doublet is more low and the right is more high
    # if c is negative, peaks must be the left of doublet is more high and the right is more low
    #
    if (JCoups+delta) == 0:
        _k = 0
    else:
        _k: float = JCoups / (JCoups+delta)

    k_small: float = 1 - _k
    k_large: float = 1 + _k

    res: list[tuple[float, float]] = []
    for freq, intensit in plist:
        # the left of doublet if J is positive
        res.append((freq + JCoups / 2, intensit / 2 * k_small))
        # the right of doublet if J is positive
        res.append((freq - JCoups / 2, intensit / 2 * k_large))
    return res

# @cachier(separate_files=True)
# def Pauil_matrix(nspins: int) -> tuple[cplex, cplex]:
#    """
#    Create Pauli matrices for a given number of spins.
#
#    This function generates the standard Pauli matrices (sigma_x, sigma_y, sigma_z)
#    scaled by 1/2, which are fundamental operators in quantum mechanics for describing
#    spin-1/2 particles.
#
#    Args:
#        nspins: Number of spins in the system
#
#    Returns:
#        Tuple containing the three Pauli matrices (sigma_x, sigma_y, sigma_z)
#    """
#
#    sigma_x: cplex = np.array([[0, 1 / 2], [1 / 2, 0]])
#    sigma_y: cplex = np.array([[0, -1j / 2], [1j / 2, 0]])
#    sigma_z: cplex = np.array([[1 / 2, 0], [0, -1 / 2]])
#    unit: cplex = np.array([[1, 0], [0, 1]])
#
#    L: cplex = np.empty(
#        (3, nspins, 2 ** nspins, 2 ** nspins), dtype=np.complex64)
#    for n in range(nspins):
#        Lx_current: cplex = np.array([1])
#        Ly_current: cplex = np.array([1])
#        Lz_current: cplex = np.array([1])
#
#        for k in range(nspins):
#            if k == n:
#                Lx_current = np.kron(Lx_current, sigma_x).astype(np.complex64)
#                Ly_current = np.kron(Ly_current, sigma_y).astype(np.complex64)
#                Lz_current = np.kron(Lz_current, sigma_z).astype(np.complex64)
#            else:
#                Lx_current = np.kron(Lx_current, unit).astype(np.complex64)
#                Ly_current = np.kron(Ly_current, unit).astype(np.complex64)
#                Lz_current = np.kron(Lz_current, unit).astype(np.complex64)
#
#        L[0][n] = Lx_current
#        L[1][n] = Ly_current
#        L[2][n] = Lz_current
#
#    L_T: cplex = L.transpose(1, 0, 2, 3)
#    Lproduct: cplex = np.tensordot(
#        L_T, L, axes=((1, 3), (0, 2))).swapaxes(1, 2).astype(np.complex64)
#
#    return L, Lproduct

# @cachier(separate_files=True)
# def T_matrix1(nspins: int) -> np_uint:
#    """
#    Generate transition matrix T for spin systems.
#
#    This function creates a binary matrix that indicates which states can transition
#    to each other in a quantum spin system, where transitions occur when only one
#    spin flips (Hamming distance of 1).
#
#    Args:
  #      nspins: Number of spins in the system
#
#    Returns:
  #      Matrix T representing possible spin transitions
#    """
#
#    n: int = 2 ** nspins
#    T: np_uint = np.zeros((n, n), dtype=np.uint8)
#    for i in range(n - 1):
  #      for j in range(i + 1, n):
    #      if bin(i ^ j).count('1') == 1:
    #      T[i, j] = 1
#    T += T.T
#    return T

# def qm_parameter1(v: list[float], J: np_float) -> tuple[cplex, np_uint]:
#    """
#    Calculate the Hamiltonian and transition matrix for a spin system.
#
#    This function constructs the angular momentum operators (Lx, Ly, Lz) for each spin
#    and builds the total Hamiltonian H from the Zeeman terms (v) and dipolar coupling terms (J).
#
#    Args:
#        v (list[float]): List of resonance frequencies in Hz for each spin.
#        J (npt.NDArray[np.float64]): Dipolar coupling matrix (Hz) with shape (nspins, nspins).
#
#    Returns:
#        tuple[npt.NDArray[np.complex128], npt.NDArray[np.float64]]:
#        - H: The total Hamiltonian matrix (complex128)
#        - T: Transition matrix for intensity calculations (float64)
#    """
#
#    L, Lproduct = Pauli_matrix(len(v))
#    T: np_uint = T_matrix(len(v))
#
#    Lz = L[2]  # array of Lz operators
#    H: cplex = np.tensordot(
#        v, Lz, axes=1).astype(np.complex128)
#    # ic(H)
#
#    scalars: np_float = 0.5 * J
#    H += np.tensordot(scalars, Lproduct, axes=2)
#
#    return H, T
