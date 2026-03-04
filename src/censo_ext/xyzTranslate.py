#!/usr/bin/env python
import argparse
import numpy as np
from pathlib import Path
from censo_ext.Tools.utility import print_arguments
from censo_ext.Tools.xyzfile import GeometryXYZs

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


def cml() -> argparse.Namespace:
    """ Get args object from commandline interface.
        Needs argparse module."""
    parser = argparse.ArgumentParser(
        # description=f"{descr}",
        # # formatter_class=argparse.RawDescriptionHelpFormatter,
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
    return parser.parse_args()


def read_xyz_file(fileName: str | Path) -> GeometryXYZs:
    try:
        geometryXYZs = GeometryXYZs()
        geometryXYZs.set_filename(fileName)
        geometryXYZs.method_read_xyz()
        return geometryXYZs

    except Exception as e:
        print(f"Failed to read file {fileName}: {e}")
        raise FileNotFoundError(f"{fileName}")


def write_xyz_file(outFile: GeometryXYZs, fileName: str | Path) -> None:
    """Write XYZ data to a file."""
    fileName = Path(fileName)
    try:
        outFile.set_filename(fileName)
        outFile.method_save_xyz([])
    except Exception as e:
        print(f"Failed to write file {fileName}: {e}")
        raise FileNotFoundError(f"{fileName}")


def main(args: argparse.Namespace = argparse.Namespace()) -> None:

    if args == argparse.Namespace():
        args = cml()
    print_arguments()

    inFile = Path(args.file)

    try:
        xyzFile: GeometryXYZs = read_xyz_file(fileName=inFile)

        outFile = GeometryXYZs()
        if args.cuts:
            outFile: GeometryXYZs = xyzFile.method_translate_cut_xyzs(
                delta=np.array(args.move), cut=args.cuts+1)
        else:
            outFile: GeometryXYZs = xyzFile.method_translate_xyzs(
                np.array(args.move))

        write_xyz_file(outFile=outFile, fileName=args.out)

    except FileNotFoundError as e:
        print(f"Input file not found: {e}")
        exit(1)

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
