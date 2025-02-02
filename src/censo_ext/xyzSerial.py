#!/usr/bin/env python3
from censo_ext.Tools.xyzfile import GeometryXYZs
from pathlib import Path
import argparse
import os
from sys import argv as sysargv
from icecream import ic
descr = """
________________________________________________________________________________
|                                          [07.07.2023] vitamin.cheng@gmail.com
| Reorder the Serial No. in xyz file
| Usage : Serial.py <geometry> [options]
| [Options]
| Input    : -i one xyz file [default isomers.xyz]
| Output   : -o one xyz file [default output.xyz]
| New      : -n or --new  To reorder the Serial No. from No. 1
| Keep     : -k or --keep To Keep the Serial No. in #Clusters
| Print    : -p Print the final data on screen 
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
        "-p",
        "--print",
        dest="print",
        action="store_true",
        default=False,
        help="Print the final data on screen (stdout)",
    )

    index_group = parser.add_mutually_exclusive_group()

    index_group.add_argument(
        "-k",
        "--keep",
        dest="keep",
        action="store_true",
        required=False,
        help="To keep the Serial No. in #Clusters",
    )

    index_group.add_argument(
        "-n",
        "--new",
        dest="new",
        action="store_true",
        required=False,
        help="To reorganize the Serial No. from No. 1",
    )

    args: argparse.Namespace = parser.parse_args()
    return args


def main(args: argparse.Namespace = argparse.Namespace()) -> None:

    if args == argparse.Namespace():
        args = cml(descr)
    if not args.print:
        print(descr)  # Program description
        print("    provided arguments: {}".format(" ".join(sysargv)))

    xyzfile: GeometryXYZs = GeometryXYZs(Path(args.file))
    xyzfile.method_read_xyz()

    if args.keep:
        if not args.print:
            print("keep serial numbers")
        xyzfile.method_comment_keep()

    elif args.new:
        if not args.print:
            print("Reordering serials numbers from 1")
        xyzfile.method_comment_new()

    if args.print:
        xyzfile.method_print([])
    else:
        xyzfile.set_filename(Path(args.out))
        xyzfile.method_save_xyz([])
        print(f"Data saved to : {args.out}")


if __name__ == "__main__":
    main()

    # For example
    # python3 xyzSerial.py -i ../tests/crest_conformers.xyz -n
    # python3 xyzSerial.py -i ../tests/crest_conformers.xyz -k
    # python3 xyzSerial.py -i ../tests/crest_conformers.xyz -n -p
    #
