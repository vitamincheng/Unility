#!/usr/bin/env python
import argparse
import os
import numpy as np
import numpy.typing as npt
from sys import argv as sysargv
# from icecream import ic
from censo_ext.Tools.xyzfile import GeometryXYZs
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
        help="Factor of threshold limit of std [default 0.50 for analysis 0.20 for Filter]",
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


def cal_RMSD_coord(args, xyzfile: GeometryXYZs, idx1_cal: list[int]) -> npt.NDArray[np.float64]:
    # start from 0 to num-1
    idx0_cal: list[int] = [x-1 for x in idx1_cal]
    from censo_ext.Tools.calculate_rmsd import cal_RMSD_xyz
    x = {"remove_idx": args.remove_idx, "add_idx": args.add_idx,
         "bond_broken": args.bond_broken, "ignore_Hydrogen": args.ignore_Hydrogen, "debug": False}
    coord: list[list[float]] = []
    for idx0 in (idx0_cal):
        coord_square, _ = cal_RMSD_xyz(
            xyzfile, idx0_cal[0]+1, idx0+1, args=argparse.Namespace(**x))
        A: list[float] = list(coord_square.values())
        coord.append(A)
    return np.array(coord, dtype=np.float64)


def FactorFilter(args) -> None:
    xyzfile: GeometryXYZs = GeometryXYZs(args.file)
    xyzfile.method_read_xyz()
    # start from 1 to num
    idx1_cal: list[int] = [x+1 for x in [*range(len(xyzfile))]]

    nConfs: int = len(xyzfile)

    if not args.thr:
        args.thr = 2
        if int(nConfs / 10) > args.thr:
            args.thr = int(nConfs / 10)

    print("")
    print(" Total conformers     : ", nConfs)
    print(" Threshold conformers : ", args.thr)
    print("")

    idx1_separate: int = 1
    idx1_minor: list[list[int]] = []
    major_idx: list[int] = []
    minor_idx: list[int] = []

    while (nConfs > args.thr):
        print(" ========== Processing ", idx1_separate, " ==========")
        np_S: npt.NDArray[np.float64] = cal_RMSD_coord(
            args, xyzfile, idx1_cal).T
        S_STD: npt.NDArray[np.float64] = np.std(np_S, axis=0)
        S_avg_STD: np.float64 = np.float64(np.average(S_STD))
        print(" Average of STD       : ", f"{S_avg_STD:10.5f}")
        print(" Factor of STD ranges : ", f"{args.factor:10.5f}")
        print(" Limits of STD        : ", f"{S_avg_STD*args.factor: 10.5f}", "\n")  # nopep8

        counter_major, counter_minor = 0, 0
        major_idx, minor_idx = [], []
        print(" CONF        STD  in major.xyz      idx     in minor.xyz    idx")
        import copy
        idx1_cal_CONF: list[int] = copy.deepcopy(idx1_cal)

        for idx in range(len(S_STD)):
            if (S_STD[idx] >= S_avg_STD*args.factor):
                print(f"{idx1_cal_CONF[idx]:>5d} {S_STD[idx]: > 10.5f}  major factor    {(counter_major+1): > 5d}")  # nopep8
                major_idx.append(idx1_cal_CONF[idx])
                counter_major += 1
            else:
                print(f"{idx1_cal_CONF[idx]:>5d} {S_STD[idx]: > 10.5f}", " "*26, f"minor factor  {(counter_minor+1): > 5d}")  # nopep8
                minor_idx.append(idx1_cal_CONF[idx])
                idx1_cal.remove(idx1_cal_CONF[idx])
                counter_minor += 1

        print("")
        print(f" Major idx: \n {major_idx}")
        print(f" Minor idx: \n {minor_idx}")
        print(""*2)
        idx1_minor.append(minor_idx)
        idx1_separate += 1
        nConfs = len(idx1_cal)

    Dir_Res: Path = Path("Final_Result")
    isExist: bool = os.path.exists(Dir_Res)
    if not isExist:
        Dir_Res.mkdir()

    print(" ========== Finally Data ==========")

    nMinor: list[int] = []
    for idx, x in enumerate(idx1_minor):
        xyzfile.set_filename(Dir_Res / Path(f"minor{str(idx+1)}.xyz"))
        xyzfile.method_save_xyz(x)
        print(f" minor{str(idx+1)}.xyz  : {x}")
        nMinor.append(len(x))

    residue_file: Path = Path("residue.xyz")
    xyzfile.set_filename(Dir_Res / residue_file)
    xyzfile.method_save_xyz(major_idx)
    print(f" {residue_file} : {major_idx}")

    np_nMinor: npt.NDArray[np.float64] = np.array(nMinor)
    print(" Coefficient of variation : ", np_nMinor.std()/np_nMinor.mean())

    import subprocess
    subprocess.call("rm -rf *_tmp", shell=True)
    subprocess.call("rm -rf CONF", shell=True)


def main(args: argparse.Namespace = argparse.Namespace()) -> None:

    if args == argparse.Namespace():
        args = cml(descr)
    print(descr)  # Program description
    print("    provided arguments: {}".format(" ".join(sysargv)))

    if args.Analysis:
        if not args.factor:
            args.factor = 0.50
        minor_list: list[int]
        Table_S: dict[int, float]
        from censo_ext.Tools.factor import method_factor_analysis, method_factor_opt
        minor_list, Table_S = method_factor_analysis(args)
        if args.opt:
            method_factor_opt(
                args=args, low_factor=minor_list, Table_S=Table_S)

    if args.Filter:
        if not args.factor:
            args.factor = 0.20
        if not args.thr:
            args.thr = 2
        FactorFilter(args)


if __name__ == "__main__":
    main()

#   test
#   python3 FactorAnalysis.py -A -i tests/data/crest_conformers.xyz --opt
#   python3 FactorAnalysis.py -F -i tests/data/crest_conformers.xyz --bb 40 44 -nh
