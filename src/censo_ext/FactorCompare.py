#!/usr/bin/env python
import argparse
import numpy as np
import numpy.typing as npt
from pathlib import Path
from censo_ext.Tools.Parameter import Eh
from censo_ext.Tools.utility import delete_all_files, print_arguments
from censo_ext.Tools.xyzfile import GeometryXYZs
from censo_ext.Tools.spectra import Boltzmann_Weighting
descr = """
________________________________________________________________________________
| Compare two different xyz files by using factor analysis (Structure Integrity)
| Usages  : FactorCompare.py <geometry> 
| Input   : -i first.xyz second.xyz 
|______________________________________________________________________________
"""


def cml() -> argparse.Namespace:
    """ Get args object from commandline interface.
        Needs argparse module."""
    parser = argparse.ArgumentParser(
        description=f"{descr}",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage=argparse.SUPPRESS,
        add_help=True
    )

    parser.add_argument(
        "-i",
        "--input",
        dest="file",
        action="store",
        required=True,
        type=str,
        nargs=2,
        help="Provide two input_file name ",
    )

    args: argparse.Namespace = parser.parse_args()
    return args


def cal_RMSD(xyzFile: GeometryXYZs, idx1_p: int, idx1_q: int) -> float:
    from censo_ext.Tools.calculate_rmsd import cal_RMSD_xyz
    _, RMSD = cal_RMSD_xyz(
        xyzFile=xyzFile, idx1_p=idx1_p, idx1_q=idx1_q,
        _remove_idx=None, _add_idx=None, _bond_broken=None, _ignore_Hydrogen=True)
    return RMSD


def Factor_xyzCompare(args) -> None:
    import subprocess
    merge_FileName: Path = Path("temp_save.xyz")
    subprocess.call(
        f"cat {args.file[0]} {args.file[1]} > {merge_FileName}", shell=True)

    xyzFile_P: GeometryXYZs = GeometryXYZs(args.file[0])
    xyzFile_P.method_read_xyz()
    xyzFile_Q: GeometryXYZs = GeometryXYZs(args.file[1])
    xyzFile_Q.method_read_xyz()
    xyzFile_Merge: GeometryXYZs = GeometryXYZs(merge_FileName)
    xyzFile_Merge.method_read_xyz()

    nSts_p: int = len(xyzFile_P)
    nSts_q: int = len(xyzFile_Q)

    result: list[float] | npt.NDArray[np.float64] = []
    for idx1_p in range(1, nSts_p+1):
        for idx1_q in range(1, nSts_q+1):
            result.append(cal_RMSD(xyzFile_Merge, idx1_p, idx1_q+nSts_p))

    result = np.array(result).reshape(nSts_p, nSts_q).T

    Energy: list[float] = [
        St._comment_energy for St in np.array(xyzFile_P.Sts)]

    np_Energy = np.array(Energy)*Eh
    np_Energy = np_Energy - np_Energy.min()
    intp_Energy: npt.NDArray[np.intp] = np.argsort(np_Energy)
    BW: npt.NDArray[np.float64] = Boltzmann_Weighting(
        np_Energy, TEMP=298.15)
    # print(np_Energy[intp_Energy])
    # print(BW[intp_Energy])

    np.set_printoptions(suppress=True)

    # Output results
    print("")
    print(" ========== Structure_Integrity_Compare ========== ")

    np_Res: npt.NDArray[np.float64] = result
    min_Res: npt.NDArray[np.float64] = np.min(np_Res, axis=0)
    idx1_Res: npt.NDArray[np.int64] = np.array([], dtype=int)

    for idx0 in range(len(np_Res[0])):
        idx1_Res = np.append(idx1_Res, np.where(
            np_Res.T[idx0] == np_Res.min(0)[idx0])[0][0]+1)

    sort_Res: npt.NDArray[np.float64] = np.copy(min_Res)
    sort_Res.sort()
    if sort_Res[-1] <= 1e-14:
        print("")
        print("  In your input two files are the same")
        print("  Exit and Close the program !!!")
        exit(0)

    second_diff_Res: npt.NDArray[np.float64] = np.diff(sort_Res, 2)

    idx0_max_second_diff_Res: npt.NDArray[np.int64] = np.array(
        [], dtype=np.int64)
    for idx0, x in enumerate(second_diff_Res):
        if x > float(second_diff_Res.std()):
            idx0_max_second_diff_Res = np.append(
                idx0_max_second_diff_Res, idx0)

    thr: float = float(sort_Res[idx0_max_second_diff_Res[0]+2])

    print(f" threhsold(thr) is : {thr}")
    print("")
    print("   P_idx Erel/kcal weight/tot     STD<thr    Q_idx         STD>thr    Q_idx")

    total_BW: float = 0
    for idx0, x in enumerate(min_Res):
        print(f"{(idx0+1):>8d}", end="")
        print(f"{(np_Energy[intp_Energy][idx0]):10.3f} {(BW[intp_Energy][idx0]):10.5f}", end="")  # nopep8

        if x < thr:
            print(f"{x:>12.5f} {(idx1_Res[idx0]):>8d}")
            total_BW = total_BW + BW[intp_Energy][idx0]
        else:
            print(" "*25, end="")
            print(f"{x:>12.5f} {idx1_Res[idx0]:>8d}")

    delete_all_files(merge_FileName)
    print("")
    print(f"Weight_total  :  {total_BW:>12.5f}")
    print("")
    print(" ========== Finished ==========")
    print("")


def main(args: argparse.Namespace = argparse.Namespace()) -> None:
    if args == argparse.Namespace():
        args = cml()
    print_arguments()

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
