#!/usr/bin/env python
# from censo_ext.Tools.ml4nmr import read_mol_neighbors
from censo_ext.Tools.utility import AtomID, print_arguments
import argparse

from censo_ext.Tools.xyzfile import GeometryXYZs
# from censo_ext.TopoAnalysis import cal_RMSD
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
    from icecream import ic
    # args.file = "/Users/chengwen-cheng/Desktop/Simulation/@@Structures/38.Ergocalciferol(Vitamin D2)/04.Censo(Hydrogen)/crest_conformers.xyz"
    # mol, neighbors = read_mol_neighbors(args.file, check=False)
    # ic(mol)
    # for x in mol:
    #    print(x)
    # ic(neighbors)

    xyzFile: GeometryXYZs = GeometryXYZs(args.file)
    xyzFile.method_read_xyz()
    from ase import Atom
    from ase import Atoms

    New: Atoms = Atoms()
    ic(xyzFile.Sts[0].names)
    for idx0, x in enumerate(xyzFile.Sts[0].coord):
        New.append(
            Atom(xyzFile.Sts[0].names[AtomID(idx0+1)], x))
    # for x in New:
    #    print(x)

    # print(cal_RMSD(xyzFile, 1, 4, bond_broken=(55, 57), check=False))
    # print(cal_RMSD(xyzFile, 1, 4, bond_broken=(57, 55), check=False))


if __name__ == "__main__":
    main()
