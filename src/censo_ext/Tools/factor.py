#!/usr/bin/env python3
from typing import Literal
from censo_ext.Tools.xyzfile import GeometryXYZs
import argparse
import numpy as np
import numpy.typing as npt
from icecream import ic
from censo_ext.Tools.calculate_rmsd import cal_rmsd_xyz


def method_factor_analysis(args) -> tuple[list[int], dict[int, float]]:

    xyzfile: GeometryXYZs = GeometryXYZs(args.file)
    xyzfile.method_read_xyz()
    args_x: dict = {"remove_idx": None, "add_idx": None,
                    "bond_broken": None, "ignore_Hydrogen": True, "debug": False, }
    coord: list[list[float]] = []
    idx_element: list[int] = []
    for idx in range(len(xyzfile)):
        coord_square, _ = cal_rmsd_xyz(
            xyzfile, 1, idx+1, args=argparse.Namespace(**args_x))
        var: list[float] = list(coord_square.values())
        if idx == 0:
            idx_element = list(coord_square.keys())
        coord.append(var)

    if len(idx_element) != 0:
        idx_STD: list[int] = idx_element
    else:
        print("idx_coord is not exist and so quit this program !!!")
        ic()
        exit(1)

    dict_idx_STD: dict[int, float] = dict(
        zip(idx_STD, np.std(np.array(coord).T, axis=1).astype(float)))
    avg_STD: np.float64 = np.float64(
        np.average(np.array(list(dict_idx_STD.values()))))

    print(" ========== Factor Analysis Processing ========== ")
    print("\n Average of STD      : ", end="")
    print(f"{avg_STD:>12.8f}")
    print(f" Threshold {
          args.factor:3.2f} *STD : {avg_STD*args.factor:>12.8}", "\n")
    print(f" Atom        STD     Major (>STD)  or    Low (<({args.factor}STD)")
    idx1_major_factor: list[int] = []
    idx1_minor_factor: list[int] = []

    for idx, x in dict_idx_STD.items():
        if (x >= avg_STD):
            print(f"{int(idx):>5d} {x:>10.5f}     Major factor")
            idx1_major_factor.append(int(idx))
        elif (x <= avg_STD*args.factor):
            print(f"{int(idx):>5d} {x:>10.5f}", " "*23, "Low factor")
            idx1_minor_factor.append(int(idx))
        else:
            print(f"{idx:>5d} {x:>10.5f}")

    print("\n Major Factor List: ", idx1_major_factor)
    print(" ========== Finished ==========")
    return idx1_minor_factor, dict_idx_STD


def method_factor_opt(args, low_factor: list[int], Table_S: dict[int, float]) -> tuple[Literal[True], list[int], float] | Literal[False]:
    """
    Optimize the broken-bond location based on calculated STD values.
    return: A tuple indicating whether the optimization is successful, the optimized broken-bond location, and the ratio.
    """

    print(" ")
    print(" ========== Optimized Broken-Bond Location Process ==========")
    from censo_ext.Tools.topo import Topo
    bonding_low_factor: list[npt.NDArray] = []
    for idx in (low_factor):
        args_x: dict = {"file": args.file, "bonding": idx,
                        "print": False, "debug": False}
        Sts_topo: Topo = Topo(args_x["file"])
        bonding_low_factor.append(
            np.array(Sts_topo.method_bonding(argparse.Namespace(**args_x))))

    Pair_low_factor: list[list[int]] = []

    for idx, x in enumerate(low_factor):
        for idy, y in enumerate(bonding_low_factor[idx]):
            Pair_low_factor.append([x, int(bonding_low_factor[idx][idy])])

    for idx, x in enumerate(Pair_low_factor):
        if x[0] > x[1]:
            x[0], x[1] = x[1], x[0]

    unique_pair_low_factor: list[list[int]] = [
        list(t) for t in set(tuple(x) for x in Pair_low_factor)]

    nCONFS: int = len(list(Table_S.keys()))
    idx_STD: dict[int, float] = Table_S

    idx_ratio: list[list[int]] = []
    list_ratio: list[float] = []
    for idx, x in enumerate(unique_pair_low_factor):

        args_x = {"file": args.file, "bond_broken": [
            x[0], x[1]], "print": False, "debug": False}
        Sts_topo: Topo = Topo(args_x["file"])
        idx_STD_L: list[int] = Sts_topo.method_broken_bond(
            argparse.Namespace(**args_x))
        args_x = {"file": args.file, "bond_broken": [
            x[1], x[0]], "print": False, "debug": False}
        idx_STD_R: list[int] = Sts_topo.method_broken_bond(
            argparse.Namespace(**args_x))

        total_STD_L: float = float(0)
        total_STD_R: float = float(0)

        if len(idx_STD_L) < 1 or len(idx_STD_R) < 1:
            print("something wrong in your List_STD ")
            ic()
            exit(1)

        elif len(idx_STD_L) < (nCONFS-2) and len(idx_STD_R) < (nCONFS-2):

            print(f" Index of atoms :      {x[0]:4d}   vs {x[1]:4d}")
            print(f" Sizes of deviation :  {int(len(idx_STD_L)): 4d}   vs {int(len(idx_STD_R)): 4d}")  # nopep8

            for y in idx_STD_L:
                total_STD_L += float(idx_STD[y])
            for y in idx_STD_R:
                total_STD_R += float(idx_STD[y])

            print(f" STD :           {total_STD_L:10.5f}   vs {total_STD_R: 10.5f}")  # nopep8
            print(f" STD/STD =    {(total_STD_L/total_STD_R):10.7f}")
            if total_STD_L/total_STD_R < 1:
                print(f" RATIO   =    {(total_STD_L/total_STD_R):10.7f}")
                idx_ratio.append([x[1], x[0]])
                list_ratio.append(total_STD_L/total_STD_R)
            else:
                print(f" RATIO   =    {(total_STD_R/total_STD_L):10.7f}")
                idx_ratio.append([x[0], x[1]])
                list_ratio.append(total_STD_R/total_STD_L)
            print("")

    print("")

    print(f" The Optimized Broken-bond location :  {idx_ratio[list_ratio.index(max(list_ratio))]}")  # nopep8
    print(f" The max ratio location :  {max(list_ratio)}")
    if max(list_ratio) <= 0.60:
        print(" Ratio is below 0.60")
        print(" It is not a good choice as broken-bond location.")
        print(" Not recommended to use it.")
        return False
    print(" ========== Finished ========== ")
    return True, idx_ratio[list_ratio.index(max(list_ratio))], max(list_ratio)
