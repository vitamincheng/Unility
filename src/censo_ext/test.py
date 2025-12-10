#!/usr/bin/env python
from censo_ext.Tools.utility import print_arguments
import argparse

from censo_ext.Tools.xyzfile import GeometryXYZs
from censo_ext.TopoAnalysis import cal_RMSD
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

    args.file = "/Users/chengwen-cheng/Desktop/Simulation/@@Structures/38.Ergocalciferol(Vitamin D2)/04.Censo(Hydrogen)/crest_conformers.xyz"
    xyzFile: GeometryXYZs = GeometryXYZs(args.file)
    xyzFile.method_read_xyz()
    print(cal_RMSD(xyzFile, 1, 4, bond_broken=(55, 57)))
    print(cal_RMSD(xyzFile, 1, 4, bond_broken=(57, 55)))


if __name__ == "__main__":
    main()
