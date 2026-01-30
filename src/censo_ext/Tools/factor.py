#!/usr/bin/env python
from typing import Literal
from censo_ext.Tools.xyzfile import GeometryXYZs
import numpy as np
import numpy.typing as npt
from censo_ext.Tools.calculate_rmsd import cal_RMSD_xyz
from censo_ext.Tools.utility import AtomID


def method_factor_analysis(xyzFile: GeometryXYZs, _factor) -> tuple[list[AtomID], dict[AtomID, float]]:
    """
    Performs factor analysis on a set of geometries to identify atoms with high and low structural variability.

    This function calculates the standard deviation of atomic positions across multiple conformations
    to determine which atoms exhibit significant structural variation (major factors) versus those with
    minimal variation (minor factors).

    Args:
        args: Command-line arguments containing file path and factor threshold.
            Expected attributes:
            - file: Path to the XYZ file containing geometries
            - factor: Threshold multiplier for determining major vs minor factors

    Returns:
        tuple[list[int],dict[int,float]]:
            A tuple containing:
            - a list of atom indices identified as having low structural variability (minor factors)
            - and a dictionary mapping atom indices to their corresponding standard deviations

    Example:
        >>> args = argparse.Namespace(file="geometries.xyz", factor=1.5)
        >>> minor_factors, std_dict = method_factor_analysis(args)
        >>> print(f"Minor factors: {minor_factors}")
        >>> print(f"Standard deviations: {std_dict}")
    """

    coord: list[list[float]] = []

    # For idxElement for the data of first xyzFile
    tmp, _ = cal_RMSD_xyz(xyzFile, 1, 1, _remove_idx=None,
                          _add_idx=None, _bond_broken=None, _ignore_Hydrogen=True)
    _Element: list[AtomID] = list(tmp.keys())

    # Get variance of coord square of all xyzFile
    for idx0 in range(len(xyzFile)):
        coord_square, _ = cal_RMSD_xyz(
            xyzFile, 1, idx0+1, _remove_idx=None, _add_idx=None, _bond_broken=None, _ignore_Hydrogen=True)
        var: list[float] = list(coord_square.values())
        coord.append(var)

    atomIDs_std: dict[AtomID, float] = dict(
        zip(_Element, np.std(np.array(coord).T, axis=1).astype(float)))
    average_std: np.float64 = np.float64(
        np.average(np.array(list(atomIDs_std.values()))))

    print(" ========== Factor Analysis Processing ========== ")
    print("\n Average of STD      : ", end="")
    print(f"{average_std:>12.8f}")
    print(f" Threshold {
          _factor:3.2f} *STD : {average_std*_factor:>12.8}", "\n")
    print(
        f" Atom        STD     Major ( >STD)  or    Low ( <{_factor:>5.5f} STD)")
    _MajorFactor: list[AtomID] = []
    _MinorFactor: list[AtomID] = []

    for atomID, atomID_std in atomIDs_std.items():
        if (atomID_std >= average_std):
            print(f"{int(atomID):>5d} {atomID_std:>10.5f}     Major factor")
            _MajorFactor.append(AtomID(int(atomID)))
        elif (atomID_std <= average_std*_factor):
            print(f"{int(atomID):>5d} {atomID_std:>10.5f}",
                  " "*23, "Low factor")
            _MinorFactor.append(AtomID(int(atomID)))
        else:
            print(f"{atomID:>5d} {atomID_std:>10.5f}")

    print(f"\n Major Factor List: {_MajorFactor}")
    print(" ========== Finished ==========")

    return _MinorFactor, atomIDs_std


def method_factor_opt(xyzFile: GeometryXYZs, _lowFactor: list[AtomID], table_std: dict[AtomID, float]) -> tuple[Literal[True], list[int], float] | Literal[False]:
    """
    Optimizes the location of a broken bond based on factor analysis results.

    This function analyzes pairs of atoms that show low structural variability and
    determines which bond breaking configuration yields the most balanced distribution
    of structural deviations across conformations.

    Args:
        args: Command-line arguments containing file path.
            Expected attributes:
            - file: Path to the XYZ file containing geometries
            - low_factor: List of atom indices with low structural variability (from factor analysis)
            - Table_S: Dictionary mapping atom indices to their standard deviations from factor analysis

    Returns:
        tuple[Literal[True],list[int],float]|Literal[False]:
            A tuple containing:
            - Boolean indicating success (True if ratio > 0.60, False otherwise)
            - List of two atom indices representing the optimized broken bond location
            - The calculated ratio value for the optimal configuration

    Example:
        >>> args = argparse.Namespace(file="geometries.xyz")
        >>> low_factors = [1, 2, 3]
        >>> std_dict = {1: 0.5, 2: 0.3, 3: 0.7}
        >>> success, bond_location, ratio = method_factor_opt(args, low_factors, std_dict)
        >>> if success:
        ...     print(f"Optimal bond location: {bond_location}")
        ...     print(f"Ratio: {ratio}")

    Note:
        The function prints detailed information about the optimization process,
        including atom indices, deviation sizes, and standard deviation ratios.
        If the maximum ratio is below 0.60, it returns False indicating
        the configuration is not recommended for use.
    """

    print(" ")
    print(" ========== Optimized Broken-bond Location Process ==========")
    from censo_ext.Tools.topo import Topo
    Bonding_LowFactor: list[npt.NDArray[np.int64]] = []
    _topo: Topo = Topo(xyzFile)

    for atomID in _lowFactor:
        Bonding_LowFactor.append(
            np.array(_topo.method_bonding(_bonding=atomID, _print=False)))

    Pair_LowFactor: list[list[int]] = []

    for idx0, atomID in enumerate(_lowFactor):
        for idy0, _ in enumerate(Bonding_LowFactor[idx0]):
            Pair_LowFactor.append([atomID, int(Bonding_LowFactor[idx0][idy0])])

    for x in Pair_LowFactor:
        if x[0] > x[1]:
            x[0], x[1] = x[1], x[0]

    unique_PairLowFactor: list[list[int]] = [
        list(t) for t in set(tuple(x) for x in Pair_LowFactor)]

    nCONFs: int = len(list(table_std.keys()))
    atomIDs_std: dict[AtomID, float] = table_std

    idx_ratio: list[list[int]] = []
    Ratio: list[float] = []
    for x in unique_PairLowFactor:

        atomIDs_L: list[AtomID] = _topo.method_broken_bond(
            _bond_broken=(x[0], x[1]), _print=False)
        atomIDs_R: list[AtomID] = _topo.method_broken_bond(
            _bond_broken=(x[1], x[0]), _print=False)

        # total std of Left fragment of inputted data
        tSTD_L: float = float(0.0)
        # total std of Right fragment of inputted data
        tSTD_R: float = float(0.0)

        if len(atomIDs_L) < 1 or len(atomIDs_R) < 1:
            raise ValueError("something wrong in your List_STD ")

        elif len(atomIDs_L) < (nCONFs-2) and len(atomIDs_R) < (nCONFs-2):

            print(f" Index of atoms :      {x[0]:4d}   vs {x[1]:4d}")
            print(f" Sizes of deviation :  {int(len(atomIDs_L)): 4d}   vs {int(len(atomIDs_R)): 4d}")  # nopep8

            for y in atomIDs_L:
                tSTD_L += float(atomIDs_std[y])
            for y in atomIDs_R:
                tSTD_R += float(atomIDs_std[y])

            print(f" STD :           {tSTD_L:10.5f}   vs {tSTD_R: 10.5f}")  # nopep8
            print(f" STD/STD =    {(tSTD_L/tSTD_R):10.7f}")
            if tSTD_L/tSTD_R < 1:
                print(f" RATIO   =    {(tSTD_L/tSTD_R):10.7f}")
                idx_ratio.append([x[1], x[0]])
                Ratio.append(tSTD_L/tSTD_R)
            else:
                print(f" RATIO   =    {(tSTD_R/tSTD_L):10.7f}")
                idx_ratio.append([x[0], x[1]])
                Ratio.append(tSTD_R/tSTD_L)
            print("")

    print("")

    print(f" The Optimized Broken-bond location :  {idx_ratio[Ratio.index(max(Ratio))]}")  # nopep8
    print(f" The max ratio location :  {max(Ratio)}")
    if max(Ratio) <= 0.60:
        print(" Ratio is below 0.60")
        print(" It is not a good choice as broken-bond location.")
        print(" Not recommended to use it.")
        return False
    print(" ========== Finished ========== ")
    return True, idx_ratio[Ratio.index(max(Ratio))], max(Ratio)


def idx_3atom_opt(xyzFile: GeometryXYZs) -> tuple[AtomID, AtomID, AtomID]:
    from censo_ext.Tools.factor import method_factor_analysis
    _LowFactor: list[AtomID]
    _Deviation: dict[AtomID, float]
    _factor = 0.50
    while (True):
        _LowFactor, _Deviation = method_factor_analysis(
            xyzFile, _factor=_factor)
        if len(_LowFactor) >= 1:
            break
        else:
            _factor = _factor * 1.1

    _Bonding: list[list[AtomID]] = []
    from censo_ext.Tools.topo import Topo
    _topo = Topo(xyzFile)
    for x in _LowFactor:
        _Bonding.append(_topo.method_bonding(_bonding=x, _print=False))
        # print(x, _Bonding)

    _3AtomID: list[list[AtomID]] = []
    for idx0, x in enumerate(_LowFactor):
        # total numbers >=3 or >2 (one of total numbers is )
        if len(_Bonding[idx0]) > 1:
            tmp: list[AtomID] = []
            tmp.append(x)
            for y in _Bonding[idx0]:
                tmp.append(y)
            _3AtomID.append(tmp)
    # print(_3AtomID)
    # print(_Bonding)
    if len(_3AtomID) == 0:
        # print(_LowFactor[0])
        # print(_Bonding[0][0])
        a_bonding = _topo.method_bonding(_bonding=_Bonding[0][0], _print=False)
        # print(a_bonding)
        a_bonding.remove(_LowFactor[0])
        # print(a_bonding[0])
        result = (_LowFactor[0], _Bonding[0][0], a_bonding[0])
        print("")
        print(f" 3 atom idx of lowest total factor {result}")  # nopep8
        print("")
        return result

    from itertools import combinations
    Combined_3AtomID: list[tuple[AtomID, AtomID, AtomID]] = []
    for x in _3AtomID:
        for y in list(combinations(x, 3)):
            Combined_3AtomID.append(y)

    idx1_Atoms: list[AtomID] = list(_Deviation.keys())
    STD_Atoms: list[float] = list(_Deviation.values())

    intp_minTotalDev: int = 0
    minTotalDev: float = 100
    for idx0, x in enumerate(Combined_3AtomID):
        TotalDevAtoms: float = 0.0
        for y in x:
            TotalDevAtoms += (STD_Atoms[idx1_Atoms.index(AtomID(y))])
        if minTotalDev > TotalDevAtoms:
            minTotalDev = TotalDevAtoms
            intp_minTotalDev = idx0
    print("")
    print(f" 3 atom idx of lowest total factor {Combined_3AtomID[intp_minTotalDev]}")  # nopep8
    print("")
    return (Combined_3AtomID[intp_minTotalDev])
