#!/usr/bin/env python3
import argparse
import os
import numpy as np
from sys import argv as sysargv
from icecream import ic
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


def cal_rmsd(xyzfile, idx_p, idx_q) -> float:
    from censo_ext.Tools.calculate_rmsd import cal_rmsd_xyz
    x: dict = {"remove_idx": None, "add_idx": None, "bond_broken": None,
               "ignore_Hydrogen": True, "quiet": True, "debug": False}
    coord_square, result_rmsd = cal_rmsd_xyz(
        xyzfile, idx_p, idx_q, args=argparse.Namespace(**x))
    return result_rmsd


def Factor_xyzCompare(args) -> None:
    import subprocess
    import re
    import sys
    subprocess.call(
        "cat "+str(args.file[0])+" "+str(args.file[1])+"> tmp_save.xyz", shell=True)
    xyzfile_P: GeometryXYZs = GeometryXYZs(args.file[0])
    xyzfile_P.method_read_xyz()
    xyzfile_Q: GeometryXYZs = GeometryXYZs(args.file[1])
    xyzfile_Q.method_read_xyz()

    xyzfile_One: GeometryXYZs = GeometryXYZs(Path("tmp_save.xyz"))
    xyzfile_One.method_read_xyz()
    nSt_P: int = len(xyzfile_P)
    nSt_Q: int = len(xyzfile_Q)

    result: list | np.ndarray = []
    for idp in range(nSt_P):
        for idq in range(nSt_Q):
            result.append(cal_rmsd(xyzfile_One, idp+1, idq+nSt_P+1))
    result = np.array(result)
    result = (result.reshape(nSt_P, nSt_Q).T)

    import shutil
    if not shutil.which("crest"):
        print(" Need the crest program !!!!")
        print(" Exit the program ")
        exit(1)

    Dir_str = "CREST_P"
    original_cwd: str = os.getcwd()
    new_cwd: str = os.getcwd()+"/"+Dir_str
    from os.path import exists
    if not exists(Dir_str):
        os.makedirs(Dir_str)

    shutil.copyfile(original_cwd+"/"+args.file[0], new_cwd+"/"+args.file[0])
    os.chdir(new_cwd)
    subprocess.call("crest "+args.file[0]+" --cregen "+args.file[0] +
                    " --rthr 0.0175 --bthr 0.003 --ethr 0.015 --ewin 40.0 > weight_P", shell=True)
    os.chdir(original_cwd)
    shutil.copyfile(new_cwd+"/weight_P", original_cwd+"/weight_P")

    file_weight: Path = Path("weight_P")
    from censo_ext.Tools.utility import is_exist
    is_exist(file_weight)

    print(" Reading the ", file_weight, " file ")

    lines: list[str] = open(file_weight, "r").readlines()

    start_idx0: int = 0
    end_idx0: int = 0
    for idx0, line in enumerate(lines):
        if re.search(r"Erel/kcal", line):
            start_idx0 = idx0+1
        if re.search(r"ensemble average energy", line):
            end_idx0 = idx0-3
    struct_crest: np.ndarray = np.array([])
    for idx0, line in enumerate(lines):
        if idx0 >= start_idx0 and idx0 <= end_idx0:
            struct_crest: np.ndarray = np.append(
                struct_crest, [float(line.split()[1]), float(line.split()[3])])
    struct_crest = np.reshape(struct_crest, (int(len(struct_crest)/2), 2))

    np.set_printoptions(suppress=True)

    # Output results
    print("")
    print(" ========== Structure_Integrity_Compare ========== ")

    np_Result: np.ndarray = result
    min_Result = np_Result.min(0)
    idx_Result: np.ndarray = np.array([], dtype=int)

    for idx in range(len(np_Result[0])):
        idx_Result = np.append(idx_Result, np.where(
            np_Result.T[idx] == np_Result.min(0)[idx])[0][0]+1)

    sort_Result: np.ndarray = np.copy(min_Result)
    sort_Result.sort()

    diff2_Result: np.ndarray = np.diff(np.diff(sort_Result))
    STD_diff2_R: float = diff2_Result.std()

    idx_max_diff2_R: np.ndarray = np.array([], dtype=int)
    for idx, num in enumerate(diff2_Result):
        if num > STD_diff2_R:
            idx_max_diff2_R: np.ndarray = np.append(idx_max_diff2_R, idx)

    thr = sort_Result[idx_max_diff2_R[0]+2]

    print(f" threhsold(thr) is : {thr}")
    print("")
    print("   P_idx Erel/kcal weight/tot     STD<thr    Q_idx         STD>thr    Q_idx")

    weight_total = 0
    for idx0, x in enumerate(min_Result):

        print(f"{(idx0+1):>8d}", end="")
        print(f"{(struct_crest.T[0][idx0]):10.3f} {(struct_crest.T[1][idx0]):10.5f}", end="")  # nopep8

        if x < thr:
            print(f"{x:>12.5f} {(idx_Result[idx0]):>8d}")
            weight_total = weight_total + struct_crest.T[1][idx0]
        else:
            print(" "*25, end="")
            print(f"{x:>12.5f} {idx_Result[idx0]:>8d}")

    print("")
    print(f"Weight_total  :  {weight_total:>12.5f}")
    print("")
    print(" ========== Finished ==========")
    print("")
    subprocess.call("rm -rf tmp_save.xyz CREST_P weight_P ", shell=True)
    print(" Removed the temp file ")


def main(args: argparse.Namespace = argparse.Namespace()) -> None:
    if args == argparse.Namespace():
        args = cml(descr)
    print(descr)  # Program description
    print("    provided arguments: {}".format(" ".join(sysargv)))

    if args.file is None or len(args.file) != 2:
        print(" Your input files are wrong (two input file) !!!")
        print(" Exit to quit the program ")
        exit(1)

    Factor_xyzCompare(args)


if __name__ == "__main__":
    main()

#
#   test
#   python3 FactorCompare.py -i tests/data/crest_conformers.xyz tests/data/crest_conformers.xyz
#
#
