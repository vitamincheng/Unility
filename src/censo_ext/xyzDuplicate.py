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
    args: argparse.Namespace = parser.parse_args()
    return args


def Duplicate_process(_rthr: float, _xyzFile: GeometryXYZs, _add_idx: list[int] | None = None):
    print(
        f"\n  The numbers of the structures in xyz files : {len(_xyzFile.Sts)}")

    idx1_Sts: list[int] = [x+1 for x in [*range(len(_xyzFile.Sts))]]

    for idx1_p in range(1, len(_xyzFile.Sts)+1):
        if len(idx1_Sts) >= 1:
            for idx1_q in idx1_Sts.copy():
                if idx1_p != idx1_q and idx1_p < idx1_q:
                    _, result_RMSD = cal_RMSD_xyz(xyzFile=_xyzFile, idx1_p=idx1_p,
                                                  idx1_q=idx1_q, _add_idx=_add_idx, _remove_idx=None, _bond_broken=None, _ignore_Hydrogen=True)
                    if result_RMSD <= _rthr:
                        idx1_Sts.remove(idx1_q)
        else:
            pass

    _xyzFile.method_xyzExtract([x-1 for x in idx1_Sts])


def main(args: argparse.Namespace = argparse.Namespace()) -> None:
    if args == argparse.Namespace():
        args = cml()
    print_arguments()

    # moment: dict = {"file": args.file, "auto": True,
    #                "atom": None, "print": False, "replace": True, "out": None}
    # x_args = argparse.Namespace(**moment)
    # xyzReturnOandZ.main(x_args)

    xyzFile: GeometryXYZs = GeometryXYZs(args.file)
    xyzFile.method_read_xyz()

    Duplicate_process(_rthr=args.rthr, _xyzFile=xyzFile)

    Dir, file = IsExists_DirFileName(Path(args.file))
    file_split: list[str] = file.split(".")
    # file_ext: str = file_split[-1]
    fileName: str = file_split[0]

    if args.out:
        xyzFile.set_filename(args.out)
    else:
        xyzFile.set_filename(fileName + "_ext.xyz")

    xyzFile.method_save_xyz([])
    print("\n  After removed duplicated the structures,")
    print(f"  The numbers of the structures in xyz files : {len(xyzFile.Sts)}")
    print(f"\n  Saved the extracted file : {xyzFile.get_fileName()}")


if __name__ == "__main__":
    main()
