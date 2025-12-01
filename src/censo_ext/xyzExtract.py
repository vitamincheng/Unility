#!/usr/bin/env python
from censo_ext.Tools.xyzfile import GeometryXYZs
from pathlib import Path
import argparse
from censo_ext.Tools.utility import print_arguments
descr = """
________________________________________________________________________________
| Extract the index numbers in xyz file
| Usage    : xyzExtract.py <geometry> [options]
| Input    : -i one xyz file [default isomers.xyz]
| Output   : -o one xyz file [default output.xyz]
| [Options]
| Index    : -d To index numbers of xyz file and index numbers from No. 1
|______________________________________________________________________________
"""


def cml() -> argparse.Namespace:
    """ Get args object from commandline interface.
        Needs argparse module."""
    parser = argparse.ArgumentParser(
        description=f"{descr}",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage=argparse.SUPPRESS,
    )

    parser.add_argument(
        "-i",
        "--input",
        dest="file",
        action="store",
        required=False,
        default="isomers.xyz",
        help="Input xyz file [default isomers.xyz]",
    )

    parser.add_argument(
        "-o",
        "--output",
        dest="out",
        action="store",
        required=False,
        default="output.xyz",
        help="Output xyz file [dafault output.xyz]",
    )

    parser.add_argument(
        "-d",
        "--index",
        dest="index",
        action="store",
        type=int,
        nargs="+",
        required=True,
        help="Index numbers of xyz file [required]",
    )

    args: argparse.Namespace = parser.parse_args()
    return args


def main(args: argparse.Namespace = argparse.Namespace()) -> None:

    if args == argparse.Namespace():
        args = cml()
    print_arguments()

    inFile: Path = Path(args.file)
    outFile: Path = Path(args.out)
    xyzFile: GeometryXYZs = GeometryXYZs(inFile)
    xyzFile.method_read_xyz()
    xyzFile.set_filename(outFile)
    xyzFile.method_save_xyz(args.index)
    print(f"Data saved to : {outFile}")


if __name__ == "__main__":
    main()
