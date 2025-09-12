#!/usr/bin/env python
import argparse
from sys import argv as sysargv
# from icecream import ic
from censo_ext.Tools.xyzfile import GeometryXYZs


descr = """
________________________________________________________________________________
| For separate or merge of molecules 
| Usages   : molManipulate.py <geometry> [options]
| [options]
| Input    : -s input xyz file to be separated to small molecule 
|               [required]
| Output   : [default in Separation folder]
| Input    : -m input xyz file to be merged to big molecule [required] 
| Output   : -o output xyz file [default output.xyz]
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
    print(f"    provided arguments: {" ".join(sysargv)}")

    if args.separate:
        inFile: GeometryXYZs = GeometryXYZs()
        inFile.set_filename(args.separate)
        inFile.method_read_xyz()
        if len(inFile) == 1:
            inFile.method_idx_molecules_xyzs(idx1=1)
        else:
            print("  Only use single conformer in your xyz file")
            print("  Exit and Close the program !!!")
            exit(0)

    elif args.merge:
        inFile0: GeometryXYZs = GeometryXYZs()
        inFile0.set_filename(args.merge[0])
        inFile0.method_read_xyz()
        inFile1: GeometryXYZs = GeometryXYZs()
        inFile1.set_filename(args.merge[1])
        inFile1.method_read_xyz()

        outFile: GeometryXYZs = GeometryXYZs()
        outFile = inFile0 + inFile1
        outFile.set_filename(args.out)
        outFile.method_save_xyz([])


if __name__ == "__main__":
    main()

    # molManipulate.py -s tests/data/crest_conformers3.xyz
    # molManipulate.py -m tests/data/Separation/1.xyz tests/data/Separation/2.xyz
