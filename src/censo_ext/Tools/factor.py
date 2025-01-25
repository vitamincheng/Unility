#!/usr/bin/env python3
from censo_ext.Tools.xyzfile import ClassGeometryXYZs
import argparse
import numpy as np
import os
from icecream import ic


def FactorAnalysis(args) -> tuple[list[int], dict]:

    xyzfile: ClassGeometryXYZs = ClassGeometryXYZs(args.file)
    xyzfile.method_read_xyz()
    import censo_ext.Tools.calculate_rmsd as calculate_rmsd
    x: dict = {"remove_idx": None, "add_idx": None,
               "bond_broken": None, "ignore_Hydrogen": True, "debug": False, }
    coord: list = []
    CONF: list = []
    for idx in range(len(xyzfile)):
        coord_square, result_rmsd = calculate_rmsd.main_xyz(
            xyzfile, 1, idx+1, args=argparse.Namespace(**x))
        var = list(coord_square.values())
        if idx == 0:
            CONF = list(coord_square.keys())
        coord.append(var)

    np_R: np.ndarray = np.array(coord).T
    S_std: np.ndarray = np.std(np_R, axis=1)

    if len(CONF) != 0:
        idx_S = CONF
    else:
        print("CONF is not exist and so quit this program !!!")
        ic()
        os._exit(0)

    Table_S: dict = dict(zip(idx_S, S_std))

    S_average_std: np.float64 = np.average(np.array(list(Table_S.values())))
    print(" ========== Factor Analysis Processing ========== ")

    print("\n Average of STD      : ", end="")
    # print("{:>12.8}".format(S_average_std))
    print(f"{S_average_std:>12.8f}")
    print(f" Threshold {
          args.factor:3.2f} *STD : {S_average_std*args.factor:>12.8}", "\n")

    idx1_major_factor: list[int] = []
    idx1_minor_factor: list[int] = []
    print(f" Atom        STD     Major (>STD)  or    Low (<({args.factor}STD)")

    for idx, x in Table_S.items():
        if (x >= S_average_std):   # type: ignore
            print(f"{int(idx):>5d} {x:>10.5f}     Major factor")
            idx1_major_factor.append(idx)
        elif (x <= S_average_std*args.factor):
            print(f"{int(idx):>5d} {x:>10.5f}", " "*23, "Low factor")
            idx1_minor_factor.append(idx)
        else:
            print(f"{idx:>5d} {x:>10.5f}")

    print("\n Major Factor List: ", idx1_major_factor)
    print(" ========== Finished ==========")
    return idx1_minor_factor, Table_S


def FactorOpt(args, np_low_factor: list[int], Table_S: dict):
    """
    Optimize the broken-bond location based on calculated STD values.
    return: A tuple indicating whether the optimization is successful, the optimized broken-bond location, and the ratio.
    """

    print(" ")
    print(" ========== Optimized Broken-Bond Location Process ==========")
    # ic(args,np_low_factor)
    # np_low_factor = minor_list
    Bonding_low_factor: list = []
    for x in (np_low_factor):
        from censo_ext.Tools.topo import topo
        args_x: dict = {"file": args.file, "bonding": x,
                        "print": False, "debug": False}
        Sts_topo: topo = topo(args_x["file"])
        Bonding_low_factor.append(
            np.array(Sts_topo.Bonding(argparse.Namespace(**args_x))))

    Pair_low_factor: list = []

    for idx, x in enumerate(np_low_factor):
        for idy, y in enumerate(Bonding_low_factor[idx]):
            Pair_low_factor.append([x, int(Bonding_low_factor[idx][idy])])

    for idx, x in enumerate(Pair_low_factor):
        if x[0] > x[1]:
            x[0], x[1] = x[1], x[0]

    temp = set(tuple(element) for element in Pair_low_factor)
    unique_pair_low_factor: list[list] = [list(t) for t in set(temp)]

    nCONFS: int = len(np.array(list(Table_S.keys())))
    idx_STD: dict = Table_S

    idx_list_ratio: list = []
    list_ratio: list = []
    for idx, x in enumerate(unique_pair_low_factor):
        from censo_ext.Tools.topo import topo

        args_x = {"file": args.file, "bond_broken": [
            x[0], x[1]], "print": False, "debug": False}
        Sts_topo: topo = topo(args_x["file"])
        STD_L: list[int] = Sts_topo.Broken_bond(argparse.Namespace(**args_x))
        args_x = {"file": args.file, "bond_broken": [
            x[1], x[0]], "print": False, "debug": False}
        STD_R: list[int] = Sts_topo.Broken_bond(argparse.Namespace(**args_x))

        total_std_L, total_std_R = 0, 0

        if len(STD_L) < 1 or len(STD_R) < 1:
            print("something wrong in your List_STD ")
            ic()
            os._exit(0)

        elif len(STD_L) < (nCONFS-2) and len(STD_R) < (nCONFS-2):

            print(f" Index of atoms :      {x[0]:4d}   vs {x[1]:4d}")
            print(f" Sizes of deviation :  {int(len(STD_L)): 4d}   vs {int(len(STD_R)): 4d}")  # nopep8

            for y in (STD_L):
                total_std_L += idx_STD[y]
            for y in (STD_R):
                total_std_R += idx_STD[y]

            print(f" STD :           {float(total_std_L):10.5f}   vs {float(total_std_R): 10.5f}")  # nopep8
            print(f" STD/STD =    {(total_std_L/total_std_R):10.7f}")
            if total_std_L/total_std_R < 1:
                print(f" RATIO   =    {(total_std_L/total_std_R):10.7f}")
                idx_list_ratio.append([x[1], x[0]])
                list_ratio.append(total_std_L/total_std_R)
            else:
                print(f" RATIO   =    {(total_std_R/total_std_L):10.7f}")
                idx_list_ratio.append([x[0], x[1]])
                list_ratio.append(total_std_R/total_std_L)
            print("")

    print("")

    print(f" The Optimized Broken-bond location :  {idx_list_ratio[list_ratio.index(max(list_ratio))]}")  # nopep8
    print(f" The max ratio location :  {max(list_ratio)}")
    if max(list_ratio) <= 0.60:
        print(" Ratio is below 0.60")
        print(" It is not a good choice as broken-bond location.")
        print(" Not recommended to use it.")
        return False
    print(" ========== Finished ========== ")
    return True, idx_list_ratio[list_ratio.index(max(list_ratio))], max(list_ratio)
