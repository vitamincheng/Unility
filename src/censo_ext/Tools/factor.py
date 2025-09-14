#!/usr/bin/env python
from typing import Literal
from censo_ext.Tools.xyzfile import GeometryXYZs
import argparse
import numpy as np
import numpy.typing as npt
from censo_ext.Tools.calculate_rmsd import cal_RMSD_xyz


def method_factor_analysis(args) -> tuple[list[int], dict[int, float]]:
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

    xyzFile: GeometryXYZs = GeometryXYZs(args.file)
    xyzFile.method_read_xyz()
    args_x: dict = {"remove_idx": None, "add_idx": None,
                    "bond_broken": None, "ignore_Hydrogen": True, "debug": False, }
    coord: list[list[float]] = []
    idxElement: list[int] = []
    for idx0 in range(len(xyzFile)):
        coord_square, _ = cal_RMSD_xyz(
            xyzFile, 1, idx0+1, args=argparse.Namespace(**args_x))
        var: list[float] = list(coord_square.values())
        if idx0 == 0:
            idxElement = list(coord_square.keys())
        coord.append(var)

    idx_dev: dict[int, float] = dict(
        zip(idxElement, np.std(np.array(coord).T, axis=1).astype(float)))
    avSTD: np.float64 = np.float64(
        np.average(np.array(list(idx_dev.values()))))

    print(" ========== Factor Analysis Processing ========== ")
    print("\n Average of STD      : ", end="")
    print(f"{avSTD:>12.8f}")
    print(f" Threshold {
          args.factor:3.2f} *STD : {avSTD*args.factor:>12.8}", "\n")
    print(f" Atom        STD     Major (>STD)  or    Low (<({args.factor}STD)")
    idx1_MajorFactor: list[int] = []
    idx1_MinorFactor: list[int] = []

    for idx, x in idx_dev.items():
        if (x >= avSTD):
            print(f"{int(idx):>5d} {x:>10.5f}     Major factor")
            idx1_MajorFactor.append(int(idx))
        elif (x <= avSTD*args.factor):
            print(f"{int(idx):>5d} {x:>10.5f}", " "*23, "Low factor")
            idx1_MinorFactor.append(int(idx))
        else:
            print(f"{idx:>5d} {x:>10.5f}")

    print(f"\n Major Factor List: {idx1_MajorFactor}")
    print(" ========== Finished ==========")
    return idx1_MinorFactor, idx_dev


def method_factor_opt(args, low_factor: list[int], Table_S: dict[int, float]) -> tuple[Literal[True], list[int], float] | Literal[False]:
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
    bonding_LowFactor: list[npt.NDArray[np.int64]] = []
    for idx in low_factor:
        args_x: dict = {"file": args.file, "bonding": idx,
                        "print": False, "debug": False}
        Sts_topo: Topo = Topo(args_x["file"])
        bonding_LowFactor.append(
            np.array(Sts_topo.method_bonding(argparse.Namespace(**args_x))))

    PairLowFactor: list[list[int]] = []

    for idx0, x in enumerate(low_factor):
        for idy0, y in enumerate(bonding_LowFactor[idx0]):
            PairLowFactor.append([x, int(bonding_LowFactor[idx0][idy0])])

    for x in PairLowFactor:
        if x[0] > x[1]:
            x[0], x[1] = x[1], x[0]

    unique_PairLowFactor: list[list[int]] = [
        list(t) for t in set(tuple(x) for x in PairLowFactor)]

    nConfs: int = len(list(Table_S.keys()))
    idx_STD: dict[int, float] = Table_S

    idx_ratio: list[list[int]] = []
    Ratio: list[float] = []
    for x in unique_PairLowFactor:

        args_x = {"file": args.file, "bond_broken": [
            x[0], x[1]], "print": False, "debug": False}
        Sts_topo: Topo = Topo(args_x["file"])
        idxSTD_L: list[int] = Sts_topo.method_broken_bond(
            argparse.Namespace(**args_x))
        args_x = {"file": args.file, "bond_broken": [
            x[1], x[0]], "print": False, "debug": False}
        idxSTD_R: list[int] = Sts_topo.method_broken_bond(
            argparse.Namespace(**args_x))

        tSTD_L: float = float(0.0)  # total STD Left Data
        tSTD_R: float = float(0.0)  # Total STD Right Data

        if len(idxSTD_L) < 1 or len(idxSTD_R) < 1:
            raise ValueError("something wrong in your List_STD ")

        elif len(idxSTD_L) < (nConfs-2) and len(idxSTD_R) < (nConfs-2):

            print(f" Index of atoms :      {x[0]:4d}   vs {x[1]:4d}")
            print(f" Sizes of deviation :  {int(len(idxSTD_L)): 4d}   vs {int(len(idxSTD_R)): 4d}")  # nopep8

            for y in idxSTD_L:
                tSTD_L += float(idx_STD[y])
            for y in idxSTD_R:
                tSTD_R += float(idx_STD[y])

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
