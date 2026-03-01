#!/usr/bin/env python
from censo_ext.Tools.xyzfile import GeometryXYZs
from pathlib import Path
import argparse
from censo_ext.Tools.utility import IsExists_DirFileName, print_arguments
descr = """
________________________________________________________________________________
| Extract the index number in xyz file
| Usage    : xyzExtract.py <geometry> [options]
| Input    : -i one xyz file [default isomers.xyz]
| Output   : -o one xyz file [default output.xyz]
| [Options]
| Index    : -d To extract the index number of xyz file and index numbers from No. 1
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
        help="Output xyz file [dafault isomers_num.xyz]",
    )

    parser.add_argument(
        "-d",
        "--index",
        dest="index",
        action="store",
        type=int,
        nargs="+",
        required=True,
        help="To extract the index number of xyz file [required]",
    )

    return parser.parse_args()


def main(args: argparse.Namespace = argparse.Namespace()) -> None:

    if args == argparse.Namespace():
        args = cml()
    print_arguments()

    inFile: Path = Path(args.file)
    xyzFile: GeometryXYZs = GeometryXYZs(inFile)
    xyzFile.method_read_xyz()

    if args.out:
        xyzFile.set_filename(args.out)
    else:
        Dir, file = IsExists_DirFileName(Path(args.file))
        fileName: str = file.split(".")[0]
        args.out = fileName + "_" + \
            "_".join([str(x) for x in args.index]) + ".xyz"
        xyzFile.set_filename(args.out)

    xyzFile.method_save_xyz(args.index)
    print(f"Data saved to : {args.out}")


if __name__ == "__main__":
    main()
