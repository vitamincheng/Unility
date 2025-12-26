#!/usr/bin/env python
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
        description="descr",
        formatter_class=argparse.RawDescriptionHelpFormatter,
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
    args: argparse.Namespace = parser.parse_args()
    return args


def main(args: argparse.Namespace = argparse.Namespace()) -> None:
    if args == argparse.Namespace():
        args = cml()
    print_arguments()

    # default to read the file
    inFile = Path(args.file)
    xyzFile: GeometryXYZs = GeometryXYZs(inFile)
    xyzFile.method_read_xyz()

    result, idx1_molbars_false = xyzFile.method_molbar(_verbose=True)


if __name__ == "__main__":
    main()
