#!/usr/bin/env python3
import argparse
import os
import numpy as np
from pathlib import Path
from icecream import ic
from censo_ext.Tools.xyzfile import GeometryXYZs
from sys import argv as sysargv

descr = """
________________________________________________________________________________
|                                          [10.05.2024] vitamin.cheng@gmail.com
| xyzTranslate.py  
| Usages   : xyzTranslate.py <geometry> [options]
| [options]
| Input    : -i input xyz file to be translate and cut if you assign -c 
| Move     : -m use a vector [x,y,z]
| Cut      : -c cut of the line 
| Output   : -o output xyz file [default output.xyz]
| Packages : Tools 
| Module   : xyzfile.py
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
    )  # argparse.RawDescriptionHelpFormatter) #,

    parser.add_argument(
        "-c",
        "--cut",
        dest="cut",
        action="store",
        required=False,
        type=int,
        help="Provide the numbers of cut (excluding origin point)",
    )

    parser.add_argument(
        "-m",
        "--move",
        dest="move",
        action="store",
        required=False,
        type=float,
        nargs=3,
        help="Provide a translation vector [x,y,z]. Required if -c is used.",
    )

    parser.add_argument(
        "-i",
        "--input",
        dest="file",
        action="store",
        required=False,
        default="traj.xyz",
        help="Input xyz file [default traj.xyz]",
    )

    parser.add_argument(
        "-o",
        "--output",
        dest="out",
        action="store",
        required=False,
        default="output.xyz",
        help="Output xyz file [default output.xyz]",
    )
    args: argparse.Namespace = parser.parse_args()
    return args


def read_xyz_file(file: str) -> GeometryXYZs:
    try:
        infile = GeometryXYZs()
        infile.set_filename(Path(file))
        infile.method_read_xyz()
        return infile

    except Exception as e:
        print(f"Failed to read file {file}: {e}")
        raise


def write_xyz_file(outfile: GeometryXYZs, file: str) -> None:
    """Write XYZ data to a file."""
    try:
        outfile.set_filename(Path(file))
        outfile.method_save_xyz([])
    except Exception as e:
        print(f"Failed to write file {file}: {e}")
        raise


def main(args: argparse.Namespace = argparse.Namespace()) -> None:

    if args == argparse.Namespace():
        args = cml(descr)
    print(descr)  # Program description
    print("    provided arguments: {}".format(" ".join(sysargv)))

    try:
        infile: GeometryXYZs = read_xyz_file(args.file)

        outfile = GeometryXYZs()
        if args.cut == None:
            outfile: GeometryXYZs = infile.method_translate_xyzs(
                np.array(args.move))
        elif args.cut:
            outfile: GeometryXYZs = infile.method_translate_cut_xyzs(
                delta=np.array(args.move), cut=args.cut+1)

        write_xyz_file(outfile, args.out)

    except Exception as e:
        print(f"An error occurred: {e}")
        exit(1)


if __name__ == "__main__":
    main()

    # For test
    # python3 xyzTranslate.py -i ../tests/crest_conformers.xyz -m 5 0 0
    # python3 xyzTranslate.py -i ../tests/crest_conformers.xyz -m 5 0 0 -c 10
    # python3 xyzTranslate.py -i output.xyz -m 0 0 5 -c 10 -o output2.xyz
