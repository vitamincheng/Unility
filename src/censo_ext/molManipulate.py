#!/usr/bin/env python
import argparse
from sys import argv as sysargv
# from icecream import ic
from censo_ext.Tools.xyzfile import GeometryXYZs


descr = """
________________________________________________________________________________
|                                          [09.11.2024] vitamin.cheng@gmail.com
| molManipulate.py  
| Usages   : molManipulate.py <geometry> [options]
| [options]
| Input    : -s input xyz file to be separated to small molecule 
|               [required]
| Output   : [default in Separation folder]
| Input    : -m input xyz file to be merged to big molecule [required] 
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

    index_group = parser.add_mutually_exclusive_group()
    index_group.add_argument(
        "-s",
        "--separate",
        dest="separate",
        action="store",
        type=str,
        help="Provide input_file name [required]",
    )

    index_group.add_argument(
        "-m",
        "--merge",
        dest="merge",
        action="store",
        nargs=2,
        type=str,
        help="Provide merge two file name [required]",
    )

    parser.add_argument(
        "-o",
        "--output",
        dest="out",
        action="store",
        required=False,
        default="output.xyz",
        help="Provide one output xyz file [default output.xyz]",
    )
    args: argparse.Namespace = parser.parse_args()
    return args


def main(args: argparse.Namespace = argparse.Namespace()) -> None:

    if args == argparse.Namespace():
        args = cml(descr)
    print(descr)  # Program description
    print("    provided arguments: {}".format(" ".join(sysargv)))

    if args.separate:
        infile: GeometryXYZs = GeometryXYZs()
        infile.set_filename(args.separate)
        infile.method_read_xyz()
        if len(infile) == 1:
            infile.method_idx_molecules_xyzs(idx1=1)
        else:
            print("Only use single conformer in your xyz file")
            exit(0)

    elif args.merge:
        infile0: GeometryXYZs = GeometryXYZs()
        infile0.set_filename(args.merge[0])
        infile0.method_read_xyz()
        infile1: GeometryXYZs = GeometryXYZs()
        infile1.set_filename(args.merge[1])
        infile1.method_read_xyz()

        outfile: GeometryXYZs = GeometryXYZs()
        outfile = infile0 + infile1
        outfile.set_filename(args.out)
        outfile.method_save_xyz([])


if __name__ == "__main__":
    main()

    # molManipulate.py -s tests/data/crest_conformers3.xyz
    # molManipulate.py -m tests/data/Separation/1.xyz tests/data/Separation/2.xyz
