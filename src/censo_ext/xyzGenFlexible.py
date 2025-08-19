#!/usr/bin/env python
import argparse
import os
import sys
# from graph import Graph
import numpy as np
import numpy.typing as npt
from icecream import ic
from censo_ext.Tools.xyzfile import GeometryXYZs
from sys import argv as sysargv
from pathlib import Path

descr = """
________________________________________________________________________________
|                                          [11.08.2024] vitamin.cheng@gmail.com
| xyzGenFlexible.py
| Usages   : xyzGenFlexible.py <geometry> [options]
| [options]
| Input    : -i xyz file [default traj.xyz]
| Output   : -o output xyz file [default output.xyz]
| Manual   : -m Manually check out the function (xtb/orca/thermo) [defalut False]
| Packages : Tools  
| Module   : xyzfile.py / ml4nmr.py / unility.py / anmrfile.py /topo.py 
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
        default="traj.xyz",
        help="Provide one input xyz file [default traj.xyz]",
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

    parser.add_argument(
        "-m",
        "--manual",
        dest="manual",
        action="store_true",
        help="Assign the splitting position of Atoms [static Atoms, rotation Atoms] [default False]",
    )

    args: argparse.Namespace = parser.parse_args()
    return args


def read_data(args) -> tuple[dict[int, npt.NDArray], list[list[int]], list[list[np.int64]], dict[int, int], dict[int, int]]:
    from censo_ext.Tools.topo import Topo
    Sts_topo: Topo = Topo(args.file)
    _, neighbor, circleMols, residualMols = Sts_topo.topology()
    idx_atomsCN: dict[int, int] = Sts_topo.get_cn()
    # ic(neighbor, circleMols, residualMols)
    # ic(atomsCN)
    from censo_ext.Tools.ml4nmr import read_mol_neighbors_bond_order
    _, _, idx_Bond_order = read_mol_neighbors_bond_order(
        args.file)
    # ic(Bond_order)
    # ic(residualMols)
    return neighbor, circleMols, residualMols, idx_Bond_order, idx_atomsCN


def get_xyzSplit(neighbor: dict[int, npt.NDArray], circleMols: list[list[int]], residualMols: list[list[np.int64]], Bond_order: dict[int, int], atomsCN: dict[int, int], flattenCircleMols: list[int]) -> dict[int, int]:
    xyzSplit: dict[int, int] = {}
    for mol in residualMols:
        mol = list(map(int, mol))
        # ic(mol)
        flexibleMols: list[int] = [
            a for a in mol if a not in flattenCircleMols]
        nodeMols: list[int] = [a for a in mol if a in flattenCircleMols]

        if len(flexibleMols) == 1:
            continue
        # ic(flexibleMols, nodeMols)
        # ic()

        flexibleMolsCNis4: list = [
            a for a in flexibleMols if atomsCN[a] == 4]
        # ic(mol, flexibleMolsCNis4, nodeMols)
        for nodeMol in nodeMols:
            argmin: int = mol.index(nodeMol)
            argmax: int = len(mol)-2
            for arg in range(argmin, argmax+1):
                if arg+1 < len(mol) and mol[arg+1] in flexibleMolsCNis4:
                    if Bond_order[mol[arg+1]] != 3:
                        if Bond_order[mol[arg]] == 3:
                            num: int = 1
                            while (Bond_order[mol[arg-num]] == 3):
                                num += 1
                            xyzSplit[mol[arg-num]] = mol[arg+1]
                            # ic(mol[arg-num], mol[arg+1])
                        else:
                            xyzSplit[mol[arg]] = mol[arg+1]
                            # ic(mol[arg], mol[arg+1])
                    else:
                        # Bond_order[mol[arg+1]] == 3
                        # next Atoms is bond_order is 3. if Atom is Carbon, it is End of molecule.
                        # nothing needs to do, Only pass
                        pass
                else:
                    # if it have two nodes, second node is set as end point.
                    # and it must be in flexibleMolCNis4 list
                    pass
    return xyzSplit


def gen_GeometryXYZs(xyzSplitDict: dict[int, int], args: argparse.Namespace) -> None:

    ic(xyzSplitDict)
    if args.manual:
        print("Assign the first number of list : ", end="")
        for key, value in xyzSplitDict.items():
            print(key, " ", end="")
        print("")
        loop: bool = True
        idx_xyzSplit: list[int] = []
        while (loop):
            pos: str = input()
            loop = False
            for idx in pos.split():
                from censo_ext.Tools.utility import function_is_int
                if (function_is_int(idx)):
                    if int(idx) not in xyzSplitDict.keys():
                        print(" Error numbers and input the data again")
                        loop = True
                    else:
                        idx_xyzSplit.append(int(idx))
                else:
                    print(" Error word and input the data again")
                    loop = True

        x: dict = {key: value for key, value in xyzSplitDict.items()
                   if key in idx_xyzSplit}
        xyzSplitDict = x

    xyzfile: GeometryXYZs = GeometryXYZs(args.file)
    xyzfile.method_read_xyz()
    file_In: Path = Path(".in.xyz")
    file_Out: Path = Path(".out.xyz")
    xyzfile.set_filename(file_In)
    xyzfile.method_save_xyz([1])

    from censo_ext.Tools.utility import move_file
    for key, value in xyzSplitDict.items():
        # ic(key, value)
        import censo_ext.xyzSplit as xyzSplit
        args_x: dict = {"file": file_In,
                        "atoms": [key, value], "cuts": 3, "print": False, "out": file_Out}
        sys.stdout = open(os.devnull, 'w')
        xyzSplit.main(argparse.Namespace(**args_x))
        sys.stdout = sys.__stdout__
        move_file(file_Out, file_In)

    move_file(file_In, args.out)
    print(f" The data is saved to {args.out} !!!")


def main(args: argparse.Namespace = argparse.Namespace()) -> None:

    if args == argparse.Namespace():
        args = cml(descr)
    print(descr)  # Program description
    print("    provided arguments: {}".format(" ".join(sysargv)))

    neighbor, circleMols, residualMols, Bond_order, atomsCN = read_data(
        args)
    flattenCircleMols: list[int] = []
    for mol in circleMols:
        flattenCircleMols += mol

    xyzSplit: dict[int, int] = get_xyzSplit(neighbor, circleMols, residualMols,
                                            Bond_order, atomsCN, flattenCircleMols)
    gen_GeometryXYZs(xyzSplit, args)


if __name__ == "__main__":
    main()

    # test
    # python3 xyzGenFlexible.py -i tests/data/crest_conformers.xyz
    # python3 xyzGenFlexible.py -i tests/data/crest_conformers.xyz -m
