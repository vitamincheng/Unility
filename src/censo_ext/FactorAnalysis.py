#!/usr/bin/env python3
import argparse
import os
import numpy as np
from sys import argv as sysargv
from icecream import ic
from Tools.xyzfile import ClassGeometryXYZs
from pathlib import Path
descr = """
________________________________________________________________________________
|                                          [10.05.2024] vitamin.cheng@gmail.com
| Usages : FactorAnalysis.py <geometry> [options]
| 1)
| Analysis  : -A Find the Factor for atom's index of broken location
| Input     : -i input file [default traj.xyz]
| Factor    : -f factor of threshold limit of std [default 0.50]
| Opt       : -o or --opt To optimize broken-bond location [default False]
| 2)
| Filter    : -F Filter the std for atom's index of broken location
| Input     : -i input file [default traj.xyz]
| Output    : Under Final_Result folder 
| Atom      : -bb Location of Bond Broken: 1(included) 2(not included)
| ignore H  : -nh Ignore hydrogen in calculated standard deviation
| Add       : --add-idx To add the atom's index
| Remove    : --remove-idx To remove the atom's index
| Factor    : -ff factor of threshold limit of std [default 0.20]
| Threshold : -t the minimum threshold number [default 2]
| Packages  : Tools 
| Module    : xyzfile.py calculate_rmsd.py factor.py
|______________________________________________________________________________
"""


def cml(descr) -> argparse.Namespace:
    """ Get args object from commandline interface.
        Needs argparse module."""
    parser = argparse.ArgumentParser(
        description="",
        #        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage=argparse.SUPPRESS,
        add_help=False
    )  # argparse.RawDescriptionHelpFormatter) #,

    parser.add_argument(
        "-i",
        "--input",
        dest="file",
        action="store",
        required=False,
        type=str,
        default="traj.xyz",
        help="Provide input_file name [default traj.xyz] ",
    )

    parser.add_argument(
        "-bb",
        "--bond-broken",
        dest="bond_broken",
        action="store",
        nargs=2,
        required=False,
        type=int,
        help="Location of Bond Broken: 1(included) 2(not included) ",
    )

    parser.add_argument(
        "--remove-idx",
        nargs="+",
        dest="remove_idx",
        action="store",
        type=int,
        required=False,
        help="To remove atom's index",
    )

    parser.add_argument(
        "--add-idx",
        nargs="+",
        dest="add_idx",
        action="store",
        type=int,
        required=False,
        help="To add atom's index",
    )

    parser.add_argument(
        "-f",
        "--factor",
        dest="factor",
        action="store",
        default=None,
        required=False,
        type=float,
        help="Factor of threshold limit of std [default 0.50]",
    )

    parser.add_argument(
        "-o",
        "--opt",
        dest="opt",
        action="store_true",
        default=False,
        help="Optimized the broken-bond location [default False]",
    )

    parser.add_argument(
        "-t",
        "--thr",
        "--threshold",
        dest="thr",
        action="store",
        required=False,
        type=int,
        help="the minimum threshold number of residue conformers [default 2]",
    )
    parser.add_argument(
        "-nh",
        "--ignore-Hydrogen",
        "--no-Hydrogen",
        dest="ignore_Hydrogen",
        action="store_true",
        default=False,
        help="Ignore Hydrogen in calculated standard deviation [default: False]",
    )

    index_group = parser.add_mutually_exclusive_group()

    index_group.add_argument(
        "-A",
        "--Analysis",
        dest="Analysis",
        action="store_true",
        help="Factor Analysis",
    )

    index_group.add_argument(
        "-F",
        "--Filter",
        dest="Filter",
        action="store_true",
        help="Factor Filter",
    )

    args: argparse.Namespace = parser.parse_args()
    return args


def cal_rmsd_coord(args, xyzfile, idx_cal) -> np.ndarray:
    idx0_cal = [x-1 for x in idx_cal]              # start from 0 to num-1
    import Tools.calculate_rmsd
    x = {"remove_idx": args.remove_idx, "add_idx": args.add_idx,
         "bond_broken": args.bond_broken, "ignore_Hydrogen": args.ignore_Hydrogen, "debug": False}
    coord: list = []
    for idx0 in (idx0_cal):
        coord_square, result_rmsd = Tools.calculate_rmsd.main_xyz(
            xyzfile, idx0_cal[0]+1, idx0+1, args=argparse.Namespace(**x))
        A = list(coord_square.values())
        coord.append(A)
    return np.array(coord)


def FactorFilter(args) -> None:
    xyzfile: ClassGeometryXYZs = ClassGeometryXYZs(args.file)
    xyzfile.method_read_xyz()
    # start from 1 to num
    idx1_cal: list[int] = [x+1 for x in [*range(len(xyzfile))]]

    conformers_major: int = len(xyzfile)

    if not args.thr:
        args.thr = 2
        if int(conformers_major / 10) > args.thr:
            args.thr = int(conformers_major / 10)

    print("")
    print(" Total conformers     : ", conformers_major)
    print(" Threshold conformers : ", args.thr)
    print("")

    idx1_sperate: int = 1
    idx1_minor_list: list = []
    major_idx, minor_idx = [], []

    while (conformers_major > args.thr):
        print(" ========== Processing ", idx1_sperate, " ==========")
        np_S: np.ndarray = cal_rmsd_coord(args, xyzfile, idx1_cal).T
        S_std: np.ndarray = np.std(np_S, axis=0)
        S_average_std: np.float64 = np.average(S_std)
        print(" Average of STD       : ", f"{S_average_std:10.5f}")
        print(" Factor of STD ranges : ", f"{args.factor:10.5f}")
        print(" Limits of STD        : ", f"{S_average_std*args.factor: 10.5f}", "\n")  # nopep8

        S_std_array: np.ndarray = S_std
        counter_major, counter_minor = 0, 0
        major_idx, minor_idx = [], []
        print(" CONF        STD  in major.xyz      idx     in minor.xyz    idx")
        import copy
        idx_cal_CONF: list = copy.deepcopy(idx1_cal)

        for idx in range(len(S_std_array)):
            if (S_std_array[idx] >= S_average_std*args.factor):
                print(f"{idx_cal_CONF[idx]:>5d} {
                      S_std_array[idx]: > 10.5f}  major factor    {(counter_major+1): > 5d}")
                major_idx.append(idx_cal_CONF[idx])
                counter_major += 1
            else:
                print(f"{idx_cal_CONF[idx]:>5d} {
                      S_std_array[idx]: > 10.5f}", " "*26, f"minor factor  {(counter_minor+1): > 5d}")
                minor_idx.append(idx_cal_CONF[idx])
                idx1_cal.remove(idx_cal_CONF[idx])
                counter_minor += 1

        print("")
        print(" Major idx: \n", major_idx)
        print(" Minor idx: \n", minor_idx)
        idx1_minor_list.append(minor_idx)
        print(""*2)
        idx1_sperate += 1
        conformers_major = len(idx1_cal)

    path: Path = Path("Final_Result")
    from Tools.utility import IsExistReturnBool
    isExist: bool = os.path.exists(path)
    if not isExist:
        path.mkdir()

    print(" ========== Finally Data ==========")

    nMinor: list = []
    for idx, x in enumerate(idx1_minor_list):
        xyzfile.set_filename(path / Path("minor"+str(idx+1)+".xyz"))
        xyzfile.method_save_xyz(x)
        print(f" minor{str(idx+1)}.xyz  : {x}")
        nMinor.append(len(x))

    residue_fileName: Path = Path("residue.xyz")
    xyzfile.set_filename(path / residue_fileName)
    xyzfile.method_save_xyz(major_idx)
    print(f" {residue_fileName} : {major_idx}")

    np_nMinor: np.ndarray = np.array(nMinor)
    print(" Coefficient of variation : ", np_nMinor.std()/np_nMinor.mean())

    import subprocess
    subprocess.call("rm -rf *_tmp", shell=True)
    subprocess.call("rm -rf CONF", shell=True)


def main(args: argparse.Namespace = argparse.Namespace()) -> None:

    if args == argparse.Namespace():
        args = cml(descr)
    print(descr)  # Program description
    print("    provided arguments: {}".format(" ".join(sysargv)))

    if args.Analysis == True:
        if args.factor == None:
            args.factor = 0.50
        minor_list: list[int]
        Table_S: dict
        from Tools.factor import FactorAnalysis, FactorOpt
        minor_list, Table_S = FactorAnalysis(args)
        if args.opt:
            FactorOpt(args=args, np_low_factor=minor_list, Table_S=Table_S)

    if args.Filter == True:
        if args.factor == None:
            args.factor = 0.20
        if args.thr == None:
            args.thr = 2
        FactorFilter(args)


if __name__ == "__main__":
    main()

#   test
#   python3 FactorAnalysis.py -A -i tests/crest_conformers.xyz --opt
