#!/usr/bin/env python
import argparse
# import os
import numpy as np
from pathlib import Path
# from icecream import ic
from censo_ext.Tools.xyzfile import GeometryXYZs
from sys import argv as sysargv

descr = """
________________________________________________________________________________
| xyzTranslate.py  
| Usages   : xyzTranslate.py <geometry> [options]
| Input    : -i input xyz file to be translate and cut if you assign -c [default traj.xyz]
| Output   : -o output xyz file [default output.xyz]
| [options]
| Move     : -m use a vector [x,y,z]
| Cut      : -c cut of the line 
|______________________________________________________________________________
"""


def cml(descr) -> argparse.Namespace:
    """ Get args object from commandline interface.
        Needs argparse module."""
    parser = argparse.ArgumentParser(
        description="",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage=argparse.SUPPRESS,
    )

    parser.add_argument(
        "-c",
        "--nCuts",
        dest="cuts",
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


def read_xyz_file(file: str | Path) -> GeometryXYZs:
    try:
        reFile = GeometryXYZs()
        reFile.set_filename(Path(file))
        reFile.method_read_xyz()
        return reFile

    except Exception as e:
        print(f"Failed to read file {file}: {e}")
        raise FileNotFoundError(f"{file}")


def write_xyz_file(outfile: GeometryXYZs, file: str | Path) -> None:
    """Write XYZ data to a file."""
    try:
        outfile.set_filename(Path(file))
        outfile.method_save_xyz([])
    except Exception as e:
        print(f"Failed to write file {file}: {e}")
        raise FileNotFoundError(f"{file}")


def main(args: argparse.Namespace = argparse.Namespace()) -> None:

    if args == argparse.Namespace():
        args = cml(descr)
    print(descr)  # Program description
    print(f"    provided arguments: {" ".join(sysargv)}")

    inFile = Path(args.file)

    try:
        xyzFile: GeometryXYZs = read_xyz_file(inFile)

        outFile = GeometryXYZs()
        if args.cuts:
            outFile: GeometryXYZs = xyzFile.method_translate_cut_xyzs(
                delta=np.array(args.move), cut=args.cuts+1)
        else:
            outFile: GeometryXYZs = xyzFile.method_translate_xyzs(
                np.array(args.move))

        write_xyz_file(outFile, Path(args.out))

    except Exception as e:
        print(f"An error occurred: {e}")
        print("  Exit and Close the program !!!")
        exit(1)


if __name__ == "__main__":
    main()

    # For test
    # xyzTranslate.py -i tests/data/crest_conformers.xyz -m 5 0 0
    # xyzTranslate.py -i tests/data/crest_conformers.xyz -m 5 0 0 -c 10
    # xyzTranslate.py -i output.xyz -m 0 0 5 -c 10 -o output2.xyz
