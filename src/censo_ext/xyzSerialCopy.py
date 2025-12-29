#!/usr/bin/env python
from pathlib import Path
from censo_ext.Tools.calculate_rmsd import cal_RMSD_xyz
from censo_ext.Tools.utility import print_arguments
import argparse
import numpy as np

from censo_ext.Tools.xyzfile import GeometryXYZs
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
        default="crest_conformers.xyz",
        help="Provide one input xyz file [default crest_conformers.xyz]",
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="out",
        action="store",
        required=False,
        default="output.xyz",
        help="Provide one output xyz file [default output.xyz]",
    )
    args: argparse.Namespace = parser.parse_args()
    return args


def main(args: argparse.Namespace = argparse.Namespace()) -> None:
    if args == argparse.Namespace():
        args = cml()
    print_arguments()

    _inFile = Path(args.file)
    _outFile = Path(args.out)
    inFile: GeometryXYZs = GeometryXYZs(_inFile)
    outFile: GeometryXYZs = GeometryXYZs(_outFile)
    inFile.method_read_xyz()

    import sys
    import os
    sys.stdout = open(os.devnull, 'w')
    outFile.method_read_xyz()
    sys.stdout = sys.__stdout__
    nInFile_Sts = len(inFile.Sts)
    nOutFile_Sts = len(outFile.Sts)

    outFile.method_Sts_append(inFile)

    print("  ===== Serial Number Copy =====")

    print(
        f"  The numbers of the structures idx1_p in {args.file} : {nInFile_Sts}")
    print(
        f"  The numbers of the structures idx1_q in {args.out} : {nOutFile_Sts}")
    print("  idx1_p is reference and copy the nClusters of idx1_p to idx1_q")
    print("")
    # print(len(outFile))
    print("  idx1_q idx1_p new_nClusters")

    for idx0_q in range(nOutFile_Sts):
        list_rmsd: list[float] = []
        for idx0_p in range(nInFile_Sts):
            _, RMSD = cal_RMSD_xyz(xyzFile=outFile, idx1_p=idx0_q+1,
                                   idx1_q=nOutFile_Sts+idx0_p+1, _remove_idx=None,
                                   _add_idx=None, _bond_broken=None, _ignore_Hydrogen=True)

            list_rmsd.append(RMSD)
        intp = np.argmin(np.array(list_rmsd))
        outFile.Sts[idx0_q].comment_Cluster = outFile.Sts[nOutFile_Sts +
                                                          intp].comment_Cluster
        print(
            f"   {idx0_q+1:5d}  {int(intp+1):5d}    {outFile.Sts[idx0_q].comment_Cluster:5d}")
    outFile.method_xyzExtract([*range(nOutFile_Sts)])
    outFile.method_rewrite_comment()
    outFile.method_save_xyz([])


if __name__ == "__main__":
    main()
