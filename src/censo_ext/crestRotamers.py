#!/usr/bin/env python
from pathlib import Path
from censo_ext.Tools.calculate_rmsd import cal_RMSD_xyz
from censo_ext.Tools.utility import IsExist_bool, print_arguments
import argparse
import numpy as np
import numpy.typing as npt

from censo_ext.Tools.xyzfile import GeometryXYZs
from censo_ext.xyzDuplicate import Duplicate_process
from censo_ext.xyzMirror import Mirror_process
descr = """
________________________________________________________________________________
| For Generation of xyz molecule
| Usages   : xyz.py <geometry> [options]
| [options]
|______________________________________________________________________________
"""


def cml() -> argparse.Namespace:
    """ Get args object from commandline interface. Needs argparse module."""
    parser = argparse.ArgumentParser(
        description="descr",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage=argparse.SUPPRESS)
    parser.add_argument(
        "-i",
        "--input",
        dest="file",
        action="store",
        required=False,
        default="crest_rotamers.xyz",
        help="Provide one input xyz file [default crest_rotamers.xyz]",
    )
    parser.add_argument(
        "-rthr",
        "--rthr",
        dest="rthr",
        action="store",
        required=False,
        type=float,
        default=0.125,
        help="the threshold of RMSD [default 0.125]",
    )
    parser.add_argument(
        "--auto",
        dest="auto",
        action="store_true",
        help="Auto mode of saved files by use xyzReturnOandZ [default False]",
    )

    parser.add_argument(
        "--add-idx",
        nargs="+",
        dest="add_idx",
        action="store",
        type=int,
        required=False,
        help="Add atom's index (for -SH -OH -NH)",
    )

    args: argparse.Namespace = parser.parse_args()
    return args


def main(args: argparse.Namespace = argparse.Namespace()) -> None:
    if args == argparse.Namespace():
        args = cml()
    print_arguments()

    inFile: Path = Path(args.file)
    xyzFile: GeometryXYZs = GeometryXYZs(inFile)
    xyzFile.method_read_xyz()
    xyzFile.method_comment_new()
    fileName = Path("cre_members")
    Result = GeometryXYZs()

    print("")
    print("  ===== Parameter ======")
    print(f"  The threshold of RMSD                 : {args.rthr}")
    print(f"  Add atom's index (for -SH -OH -NH)    : {args.add_idx}")
    print("")

    if IsExist_bool(fileName):
        _inData: npt.NDArray[np.int64] = np.genfromtxt(
            fileName, skip_header=1, dtype=int)
        counter = 0
        for _, idx1_start, idx1_end in _inData:
            counter = counter + 1
            import copy
            outFile: GeometryXYZs = copy.deepcopy(xyzFile)
            outFile.method_xyzExtract([*range(idx1_start - 1, idx1_end)])
            outFile.set_filename(f"{idx1_start}_{idx1_end}.xyz")
            if args.auto:
                outFile.method_xyzReturnOandZ_auto()
            outFile.method_save_xyz([])
            print(f"  ===== {counter} =====")
            print(f"  Data saved to : {idx1_start}_{idx1_end}.xyz")

            Duplicate_process(_rthr=args.rthr, _xyzFile=outFile,
                              _add_idx=args.add_idx)
            outFile.set_filename(f"{idx1_start}_{idx1_end}_ext.xyz")
            print(f"  Data saved to : {idx1_start}_{idx1_end}_ext.xyz")
            print("")
            outFile.method_save_xyz([])

            StFile: GeometryXYZs = copy.deepcopy(outFile)
            StFile.method_xyzExtract([0])
            import sys
            import os
            sys.stdout = open(os.devnull, 'w')
            Mirror_process(StFile, _atom=None)
            sys.stdout = sys.__stdout__

            if len(StFile) == 0:
                continue
            else:
                outFile.method_Sts_append(StFile)

            numbers: list[int] = [*range(2, len(outFile))]
            _idx1: list[int] = [
                outFile.Sts[x-1].comment_Cluster for x in [*range(1, len(outFile))]]
            if len(_idx1) == 1:
                print("  Only one confomer in this file, and pass this")
            else:
                print(_idx1)

            if len(numbers) > 0:
                print("  idx1_p       #    |  idx1_q       #                    RMSD")
                for idx1_q in numbers:
                    idx1_p = len(outFile)
                    _, result_rmsd = cal_RMSD_xyz(xyzFile=outFile, idx1_p=idx1_p, idx1_q=idx1_q, _add_idx=None, _remove_idx=None,
                                                  _bond_broken=None, _ignore_Hydrogen=True)
                    print(
                        f"      1'   {outFile.Sts[1-1].comment_Cluster:5d}    |     {idx1_q:3d}   {outFile.Sts[idx1_q-1].comment_Cluster:5d}                {result_rmsd:12.6f}", end="")
                    if result_rmsd <= args.rthr:
                        # print(" ", len(Result.Sts), end="")
                        Result.Sts.append(outFile.Sts[idx1_q-1])
                        # print(" ", len(Result.Sts), end="")
                        # for x in Result.Sts:
                        #    print(x.comment_Cluster, end=" ")

                        print(
                            f"    Added idx1_q {idx1_q} #{outFile.Sts[idx1_q-1].comment_Cluster} to append.xyz")
                    else:
                        print("")
            print("")
            print("")

        Result.set_filename(Path("append.xyz"))
        Result.method_save_xyz([])

    else:
        print(f"  Your file {fileName} is not Exist !!! ")
        print("  Close and Exit the program !!!")
        exit(0)


if __name__ == "__main__":
    main()
