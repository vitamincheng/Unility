#!/usr/bin/env python3
from Tools.xyzfile import ClassGeometryXYZs
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
        help="Provide one xyz file to reorder the Serial No. [default isomers.xyz]",
    )

    parser.add_argument(
        "-o",
        "--output",
        dest="out",
        action="store",
        required=False,
        default="output.xyz",
        help="Provide one xyz file to be saved of the data [dafault output.xyz]",
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

    xyzfile: ClassGeometryXYZs = ClassGeometryXYZs(args.file)
    xyzfile.method_read_xyz()

    if args.keep:
        if not args.print:
            print("keep")
        xyzfile.method_comment_keep()

    elif args.new:
        if not args.print:
            print("new")
        xyzfile.method_comment_new()
    else:
        print(" Either --new or --keep in your arugment.")
        print(" Exit the program !!!")
        ic()
        os._exit(0)

    if args.print:
        xyzfile.method_print([])
    else:
        print(f"    Saved to {args.out}")
        xyzfile.set_filename(args.out)
        xyzfile.method_save_xyz([])


if __name__ == "__main__":
    main()

    # For example
    # python3 xyzSerial.py -i ../tests/crest_conformers.xyz -n
    # python3 xyzSerial.py -i ../tests/crest_conformers.xyz -k
    # python3 xyzSerial.py -i ../tests/crest_conformers.xyz -n -p
    #
