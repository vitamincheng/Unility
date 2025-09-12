#!/usr/bin/env python
import argparse
import numpy as np
import numpy.typing as npt
from sys import argv as sysargv
# from icecream import ic
from censo_ext.Tools.xyzfile import GeometryXYZs
from pathlib import Path
descr = """
________________________________________________________________________________
| For Analysis Factor of molecule xyz files 
| Usages         : FactorAnalysis.py <geometry> [options]
| [options]
| 1) Analysis    : -A Find the Factor for atom's index of broken location
|      Input     : -i input file [default traj.xyz]
|      Factor    : -f factor of threshold limit of std [default 0.50]
|      Opt       : -o or --opt To optimize broken-bond location [default False]
| 2) Filter      : -F Filter the std for atom's index of broken location
|      Input     : -i input file [default traj.xyz]
|      Output    : Under Final_Result folder 
|      Atom      : -bb Location of Bond Broken: 1(included) 2(not included)
|      ignore H  : -nh Ignore hydrogen in calculated standard deviation [default False]
|      Add       : --add-idx To add the atom's index
|      Remove    : --remove-idx To remove the atom's index
|      Factor    : -ff factor of threshold limit of std [default 0.20]
|      Threshold : -t the minimum threshold number [default 2]
|______________________________________________________________________________
"""


def cml() -> argparse.Namespace:
    """ Get args object from commandline interface.
        Needs argparse module."""
    parser = argparse.ArgumentParser(
        description="",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage=argparse.SUPPRESS,
        add_help=False
    )

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
        help="Ignore Hydrogen in calculated standard deviation [default False]",
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


def cal_RMSD_coord(args, xyzFile: GeometryXYZs, idx1_cal: list[int]) -> npt.NDArray[np.float64]:
    # start from 0 to num-1
    idx0_cal: list[int] = [x-1 for x in idx1_cal]
    from censo_ext.Tools.calculate_rmsd import cal_RMSD_xyz
    x = {"remove_idx": args.remove_idx, "add_idx": args.add_idx,
         "bond_broken": args.bond_broken, "ignore_Hydrogen": args.ignore_Hydrogen, "debug": False}
    coordSquare: list[list[float]] = []
    for idx0 in (idx0_cal):
        idx_coordSquare, _ = cal_RMSD_xyz(
            xyzFile, idx0_cal[0]+1, idx0+1, args=argparse.Namespace(**x))
        A: list[float] = list(idx_coordSquare.values())
        coordSquare.append(A)
    return np.array(coordSquare, dtype=np.float64)


def FactorFilter(args) -> None:
    xyzFile: GeometryXYZs = GeometryXYZs(args.file)
    xyzFile.method_read_xyz()
    # start from 1 to num
    idx1_xyz: list[int] = [x+1 for x in [*range(len(xyzFile))]]

    nConfs: int = len(xyzFile)

    if not args.thr:
        args.thr = 2
        if int(nConfs / 10) > args.thr:
            args.thr = int(nConfs / 10)

    print("")
    print(f" Total conformers     : {nConfs}")
    print(f" Threshold conformers : {args.thr}")
    print("")

    idx1_separate: int = 1
    idx1_minor: list[list[int]] = []
    major_idx: list[int] = []
    minor_idx: list[int] = []

    while (nConfs > args.thr):
        print(f" ========== Processing {idx1_separate} ==========")
        coord_STD: npt.NDArray[np.float64] = cal_RMSD_coord(
            args, xyzFile, idx1_xyz).T
        Column_STD: npt.NDArray[np.float64] = np.std(coord_STD, axis=0)
        Average_STD: np.float64 = np.float64(np.average(Column_STD))
        print(f" Average of STD       : {Average_STD:10.5f}")
        print(f" Factor of STD ranges : {args.factor:10.5f}")
        print(f" Limits of STD        : {Average_STD*args.factor: 10.5f} \n")  # nopep8

        counter_major, counter_minor = 0, 0
        major_idx, minor_idx = [], []
        print(" CONF        STD  in major.xyz      idx     in minor.xyz    idx")
        import copy
        idx1_calc: list[int] = copy.deepcopy(idx1_xyz)

        for idx in range(len(Column_STD)):
            if (Column_STD[idx] >= Average_STD*args.factor):
                print(f"{idx1_calc[idx]:>5d} {Column_STD[idx]: > 10.5f}  major factor    {(counter_major+1): > 5d}")  # nopep8
                major_idx.append(idx1_calc[idx])
                counter_major += 1
            else:
                print(f"{idx1_calc[idx]:>5d} {Column_STD[idx]: > 10.5f}", " "*26, f"minor factor  {(counter_minor+1): > 5d}")  # nopep8
                minor_idx.append(idx1_calc[idx])
                idx1_xyz.remove(idx1_calc[idx])
                counter_minor += 1

        print("")
        print(f" Major idx: \n {major_idx}")
        print(f" Minor idx: \n {minor_idx}")
        print(""*2)
        idx1_minor.append(minor_idx)
        idx1_separate += 1
        nConfs = len(idx1_xyz)

    reDir: Path = Path("Final_Result")
    if not reDir.is_dir():
        reDir.mkdir()

    print(" ========== Finally Data ==========")

    nMinor: list[int] = []
    for idx1, x in enumerate(idx1_minor, 1):
        xyzFile.set_filename(reDir / Path(f"minor{idx1}.xyz"))
        xyzFile.method_save_xyz(x)
        print(f" minor{idx1}.xyz  : {x}")
        nMinor.append(len(x))
    np_nMinor: npt.NDArray[np.float64] = np.array(nMinor)

    residueFile: Path = Path("residue.xyz")
    xyzFile.set_filename(reDir / residueFile)
    xyzFile.method_save_xyz(major_idx)
    print(f" {residueFile} : {major_idx}")

    print(f" Coefficient of variation : {np_nMinor.std()/np_nMinor.mean()}")

    import subprocess
    subprocess.call("rm -rf *_tmp", shell=True)
    subprocess.call("rm -rf CONF", shell=True)


def main(args: argparse.Namespace = argparse.Namespace()) -> None:

    if args == argparse.Namespace():
        args = cml()
    print(descr)  # Program description
    print(f"    provided arguments: {" ".join(sysargv)}")

    if args.Analysis:
        if not args.factor:
            args.factor = 0.50
        idx_minor: list[int]
        idx_dev: dict[int, float]
        from censo_ext.Tools.factor import method_factor_analysis, method_factor_opt
        idx_minor, idx_dev = method_factor_analysis(args)
        if args.opt:
            method_factor_opt(
                args=args, low_factor=idx_minor, Table_S=idx_dev)

    if args.Filter:
        if not args.factor:
            args.factor = 0.20
        if not args.thr:
            args.thr = 2
        FactorFilter(args)


if __name__ == "__main__":
    main()

#   test
#   FactorAnalysis.py -A -i tests/data/crest_conformers.xyz --opt
#   FactorAnalysis.py -F -i tests/data/crest_conformers.xyz --bb 40 44 -nh
