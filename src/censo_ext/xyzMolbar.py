#!/usr/bin/env python
from copy import deepcopy
from pathlib import Path
from censo_ext.Tools.utility import print_arguments
import argparse
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
        # description=f"{descr}",
        # formatter_class=argparse.RawDescriptionHelpFormatter,
        usage=argparse.SUPPRESS)
    parser.add_argument(
        "-i",
        "--input",
        dest="file",
        action="store",
        required=False,
        default="traj.xyz",
        help="Provide one input xyz file [default traj.xyz]",
    )

    parser.add_argument(
        "-o",
        "--output",
        dest="out",
        action="store",
        required=False,
        default="output.xyz",
        help="Provide one input xyz file [default output.xyz]",
    )

    return parser.parse_args()


def main(args: argparse.Namespace = argparse.Namespace()) -> None:
    if args == argparse.Namespace():
        args = cml()
    print_arguments()

    # default to read the file
    inFile = Path(args.file)
    outFile = Path(args.out)
    xyzFile: GeometryXYZs = GeometryXYZs(inFile)
    xyzFile.method_read_xyz()
    res_File = deepcopy(xyzFile)

    result, idx1_molbars_false = xyzFile.method_molbar(_verbose=True)
    # print(idx1_molbars_false)
    if result is False:
        list_Sts = set([*range(len(xyzFile.Sts))])
        if idx1_molbars_false is not None:
            diff_idx1 = list_Sts.difference(
                set([x-1 for x in idx1_molbars_false]))
            xyzFile.method_xyzExtract(list(diff_idx1))
            xyzFile.set_filename(outFile)
            xyzFile.method_save_xyz([])
            print(f"  Saved the file in {outFile}")
            res_File.method_xyzExtract(
                list(set([x-1 for x in idx1_molbars_false])))
            res_File.set_filename("residue.xyz")
            res_File.method_save_xyz([])
            print("  Saved the deleted file in residue.xyz")


if __name__ == "__main__":
    main()
