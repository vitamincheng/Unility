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
| Usages : FactorCompare.py <geometry> [options]
| Structure_Integrity_Compare 
| Input     : -i first.xyz second.xyz 
| Packages  : Tools 
| Module    : xyzfile.py / calculate_rmsd.py
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
        nargs=2,
        help="Provide two input_file name ",
    )

    args: argparse.Namespace = parser.parse_args()
    return args


def cal_RMSD(xyzfile, idx_p, idx_q) -> float:
    from censo_ext.Tools.calculate_rmsd import cal_RMSD_xyz
    x: dict = {"remove_idx": None, "add_idx": None, "bond_broken": None,
               "ignore_Hydrogen": True, "quiet": True, "debug": False}
    _, RMSD = cal_RMSD_xyz(
        xyzfile, idx_p, idx_q, args=argparse.Namespace(**x))
    return RMSD


def Factor_xyzCompare(args) -> None:
    import subprocess
    import re
    merge_FileName: Path = Path("temp_save.xyz")
    subprocess.call(
        "cat "+str(args.file[0])+" "+str(args.file[1]) + "> " + str(merge_FileName), shell=True)
    xyzfile_P: GeometryXYZs = GeometryXYZs(args.file[0])
    xyzfile_P.method_read_xyz()
    xyzfile_Q: GeometryXYZs = GeometryXYZs(args.file[1])
    xyzfile_Q.method_read_xyz()

    xyzfile_Merge: GeometryXYZs = GeometryXYZs(merge_FileName)
    xyzfile_Merge.method_read_xyz()
    nSts_P: int = len(xyzfile_P)
    nSts_Q: int = len(xyzfile_Q)

    result: list | npt.NDArray[np.float64] = []
    for idx_P in range(nSts_P):
        for idx_Q in range(nSts_Q):
            result.append(cal_RMSD(xyzfile_Merge, idx_P+1, idx_Q+nSts_P+1))

    result = np.array(result)
    result = result.reshape(nSts_P, nSts_Q).T

    import shutil
    prog: str = "crest"
    if not shutil.which(prog):
        print(f"  Need the {prog} program !!!")
        print("  Exit and Close the program !!!")
        exit(1)

    Dir_str: Path = Path("CREST_P")
    workDir: Path = Path(os.getcwd())
    Compare_Res: Path = Path("weight_P")
    New_cwd: Path = workDir / Dir_str
    from os.path import exists
    if not exists(New_cwd):
        os.makedirs(New_cwd)
    from censo_ext.Tools.utility import IsExists_DirFileName
    fileName_Dir, FileName_str = IsExists_DirFileName(Path(args.file[0]))
    FileName: Path = Path(FileName_str)
    shutil.copyfile(workDir / fileName_Dir / FileName, New_cwd / FileName)

    os.chdir(New_cwd)
    subprocess.call(prog + " " + str(FileName) + " --cregen " + str(FileName) +
                    " --rthr 0.0175 --bthr 0.003 --ethr 0.015 --ewin 40.0 > " + str(Compare_Res), shell=True)
    os.chdir(workDir)
    shutil.copyfile(New_cwd / Compare_Res, workDir / Compare_Res)

    Path_weight_P: Path = New_cwd / Compare_Res
    from censo_ext.Tools.utility import IsExist
    IsExist(Path_weight_P)

    print(" Reading the ", Path_weight_P, " file ")

    lines: list[str] = open(Path_weight_P, "r").readlines()
    start_idx0: int = 0
    end_idx0: int = 0
    for idx0, line in enumerate(lines):
        if re.search(r"Erel/kcal", line):
            start_idx0 = idx0 + 1
        if re.search(r"ensemble average energy", line):
            end_idx0 = idx0 - 3
    St_crest: npt.NDArray[np.float64] = np.array([])
    for idx0, line in enumerate(lines):
        if idx0 >= start_idx0 and idx0 <= end_idx0:
            St_crest = np.append(
                St_crest, [float(line.split()[1]), float(line.split()[3])])
    St_crest = np.reshape(St_crest, (int(len(St_crest)/2), 2))

    np.set_printoptions(suppress=True)

    # Output results
    print("")
    print(" ========== Structure_Integrity_Compare ========== ")

    np_Res: npt.NDArray[np.float64] = result
    min_Res: npt.NDArray[np.float64] = np_Res.min(0)
    idx_Res: npt.NDArray[np.float64] = np.array([], dtype=int)

    for idx in range(len(np_Res[0])):
        idx_Res = np.append(idx_Res, np.where(
            np_Res.T[idx] == np_Res.min(0)[idx])[0][0]+1)

    sort_Res: npt.NDArray[np.float64] = np.copy(min_Res)
    sort_Res.sort()
    if sort_Res[-1] <= 1e-14:
        print("")
        print("  In your input two files are the same")
        print("  Exit and Close the program !!!")
        exit(0)

    diff2_Res: npt.NDArray[np.float64] = np.diff(np.diff(sort_Res))
    STD_diff2_Res: float = float(diff2_Res.std())

    idx_max_diff2_R: npt.NDArray[np.float64] = np.array([], dtype=int)
    for idx, num in enumerate(diff2_Res):
        if num > STD_diff2_Res:
            idx_max_diff2_R = np.append(idx_max_diff2_R, idx)

    thr: float = float(sort_Res[idx_max_diff2_R[0]+2])

    print(f" threhsold(thr) is : {thr}")
    print("")
    print("   P_idx Erel/kcal weight/tot     STD<thr    Q_idx         STD>thr    Q_idx")

    weight_total = 0
    for idx0, x in enumerate(min_Res):

        print(f"{(idx0+1):>8d}", end="")
        print(f"{(St_crest.T[0][idx0]):10.3f} {(St_crest.T[1][idx0]):10.5f}", end="")  # nopep8

        if x < thr:
            print(f"{x:>12.5f} {(idx_Res[idx0]):>8d}")
            weight_total = weight_total + St_crest.T[1][idx0]
        else:
            print(" "*25, end="")
            print(f"{x:>12.5f} {idx_Res[idx0]:>8d}")

    print("")
    print(f"Weight_total  :  {weight_total:>12.5f}")
    print("")
    print(" ========== Finished ==========")
    print("")
    subprocess.call("rm -rf "+str(Dir_str)+" " + str(Compare_Res) +
                    " " + str(merge_FileName), shell=True)
    # subprocess.call("rm -rf tmp_save.xyz weight_P ", shell=True)
    print(" Removed the temp file ")


def main(args: argparse.Namespace = argparse.Namespace()) -> None:
    if args == argparse.Namespace():
        args = cml(descr)
    print(descr)  # Program description
    print("    provided arguments: {}".format(" ".join(sysargv)))

    if args.file is None or len(args.file) != 2:
        print("  Your input files are wrong (two input file) !!!")
        print("  Exit and Close the program !!!")
        exit(1)

    Factor_xyzCompare(args)


if __name__ == "__main__":
    main()

#
#   test
#   python FactorCompare.py -i tests/data/crest_conformers.xyz tests/data/crest_conformers.xyz
#
#
