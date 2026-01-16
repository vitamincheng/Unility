#!/usr/bin/env python
from pathlib import Path
# from censo_ext import xyzReturnOandZ
from censo_ext.Tools.calculate_rmsd import cal_RMSD_xyz
from censo_ext.Tools.utility import IsExists_DirFileName, print_arguments
import argparse
from censo_ext.Tools.xyzfile import GeometryXYZs
descr = """
    ________________________________________________________________________________
    | For remove Duplicate confomrers by RMSD 
    | Usages   : xyzDuplicate.py <geometry> [options]
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
        default="isomers.xyz",
        help="Provide one input xyz file [default isomers.xyz]",
    )
    parser.add_argument(
        "-o",
        "--out",
        dest="out",
        action="store",
        required=False,
        help="Provide one input xyz file [default input_ext.xyz]",
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


def Duplicate_process(_rthr: float, _xyzFile: GeometryXYZs, _add_idx: list[int] | None = None):

    idx1_Sts: list[int] = [x+1 for x in [*range(len(_xyzFile.Sts))]]

    print("")
    print("  The following list has the same structure as in RMSD and")
    print("  the list of idx1_q will be removed in the next step.")
    print("")
    print("  idx1_p       #    |  idx1_q       #                    RMSD")
    for idx1_p in range(1, len(_xyzFile.Sts)+1):
        if len(idx1_Sts) >= 1:
            for idx1_q in idx1_Sts.copy():
                if idx1_p != idx1_q and idx1_p < idx1_q:
                    _, result_RMSD = cal_RMSD_xyz(xyzFile=_xyzFile, idx1_p=idx1_p,
                                                  idx1_q=idx1_q, _add_idx=_add_idx, _remove_idx=None, _bond_broken=None, _ignore_Hydrogen=True)
                    if result_RMSD <= _rthr:
                        print(
                            f"   {idx1_p:5d}   {_xyzFile.Sts[idx1_p-1].comment_Cluster:5d}    |   {idx1_q:5d}   {_xyzFile.Sts[idx1_q-1].comment_Cluster:5d}", end="")
                        print(f"                    {result_RMSD:5.3f}")
                        idx1_Sts.remove(idx1_q)
                        # print(f" {idx1_p:5d}    | {idx1_q:5d}")
        else:
            pass

    _xyzFile.method_xyzExtract([x-1 for x in idx1_Sts])


def main(args: argparse.Namespace = argparse.Namespace()) -> None:
    if args == argparse.Namespace():
        args = cml()
    print_arguments()

    xyzFile: GeometryXYZs = GeometryXYZs(args.file)
    xyzFile.method_read_xyz()

    print(
        f"\n  The numbers of the structures in xyz file is {len(xyzFile.Sts)}")

    Duplicate_process(_rthr=args.rthr, _xyzFile=xyzFile,
                      _add_idx=args.add_idx)

    if args.out:
        xyzFile.set_filename(args.out)
    else:
        Dir, file = IsExists_DirFileName(Path(args.file))
        fileName: str = file.split(".")[0]
        xyzFile.set_filename(fileName + "_ext.xyz")

    xyzFile.method_save_xyz([])
    print("\n  After removing duplicate structures,")
    print(
        f"  The numbers of the structures in xyz file is {len(xyzFile.Sts)}")
    print(f"\n  Saved reduced file : {xyzFile.get_fileName()}")


if __name__ == "__main__":
    main()
