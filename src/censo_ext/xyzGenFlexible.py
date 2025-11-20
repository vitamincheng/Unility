#!/usr/bin/env python
import argparse
import numpy as np
import numpy.typing as npt
from icecream import ic
from censo_ext.Tools.utility import AtomID, print_arguments
from censo_ext.Tools.xyzfile import GeometryXYZs
from pathlib import Path

descr = """
________________________________________________________________________________
| For Generation of Flexible xyz molecule 
| Usages   : xyzGenFlexible.py <geometry> [options]
| Input    : -i xyz file (only for 1st xyz file) [default traj.xyz]
| Output   : -o output xyz file [default output.xyz]
| [options]
| Manual   : -m Manually check out the function (xtb/orca/thermo) [defalut False]
| nCuts    : -c or cut Number of cut to make 360 degrees around the roation axis [default 3]
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
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        help="Verbose mode [default False]",
    )

    parser.add_argument(
        "-c",
        "--nCuts",
        dest="cuts",
        action="store",
        type=int,
        default=3,
        required=False,
        help="Number of cuts to make in 360 degrees around the rotation axis [default 3]",
    )

    args: argparse.Namespace = parser.parse_args()
    return args


def read_data(args) -> tuple[dict[AtomID, npt.NDArray[np.int64]], list[list[int]], list[list[np.int64]], dict[AtomID, int], dict[int, int], dict]:
    from censo_ext.Tools.topo import Topo
    from censo_ext.Tools.ml4nmr import read_mol_neighbors_bond_order
    Sts_topo: Topo = Topo(args.file)
    _, neighbor, circleMols, residualMols, residualMols_all_pairs = Sts_topo.topology()
    idx_atomsCN: dict[int, int] = Sts_topo.get_cn()
    if args.verbose:
        ic(neighbor, circleMols, residualMols)
        ic(idx_atomsCN)
    *_, idx_Bond_order = read_mol_neighbors_bond_order(args.file)
    if args.verbose:
        ic(idx_Bond_order)
        ic(residualMols)
    return neighbor, circleMols, residualMols, idx_Bond_order, idx_atomsCN, residualMols_all_pairs


def get_xyzSplit(residualMols: list[list[np.int64]], Bond_order: dict[AtomID, int], atomsCN: dict[int, int], flattenCircleMols: list[int], residualMols_all_pairs) -> dict[int, int]:
    xyzSplit: dict[int, int] = {}
    for Mol in residualMols:
        mol: list[int] = list(map(int, Mol))
        flexibleMols: list[int] = [
            a for a in mol if a not in flattenCircleMols]
        nodeMols: list[int] = [a for a in mol if a in flattenCircleMols]
        # ic(flexibleMols, nodeMols)
        if len(flexibleMols) == 1:
            continue
        mol = nodeMols+flexibleMols
        # ic(flexibleMols, nodeMols)

        flexibleMolsCNis4: list = [
            a for a in flexibleMols if atomsCN[a] == 4]
        # ic(mol, flexibleMolsCNis4, nodeMols)
        # ic(nodeMols)
        if len(nodeMols) == 1:
            # ic(residualMols_all_pairs[nodeMols[0]])
            a = residualMols_all_pairs[nodeMols[0]].values()
            import math
            b = [x for x in a if not math.isinf(x)]
            for x in range(0, max(b)-1):
                # print(x, x+1)
                out_key: int = 0
                out_value: int = 0
                if x == 0:
                    # print("key: ", nodeMols[0])
                    out_key = nodeMols[0]
                for key, distance in residualMols_all_pairs[nodeMols[0]].items():
                    if distance == x:
                        if len([c for c in residualMols_all_pairs[key].values() if c == 1]) != 1:
                            # print("key: ", key)
                            out_key = key
                    if distance == x+1:
                        if len([c for c in residualMols_all_pairs[key].values() if c == 1]) != 1:
                            # print("value: ", key)
                            out_value = key
                if out_key not in flexibleMolsCNis4 and out_value not in flexibleMolsCNis4:
                    pass
                else:
                    # ic(out_key, out_value)
                    xyzSplit[out_key] = out_value

    return xyzSplit


def gen_GeometryXYZs(xyzSplitDict: dict[int, int], args: argparse.Namespace) -> None:

    if args.verbose:
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

    inFile = Path(args.file)
    outFile = Path(args.out)

    xyzFile: GeometryXYZs = GeometryXYZs(inFile)
    xyzFile.method_read_xyz()
    splitIn: Path = Path(".in.xyz")
    splitOut: Path = Path(".out.xyz")
    xyzFile.set_filename(splitIn)
    # only read first xyz file
    xyzFile.method_save_xyz([1])

    from censo_ext.Tools.utility import move_file
    for key, value in xyzSplitDict.items():
        if args.verbose:
            ic(key, value)
        import censo_ext.xyzSplit as xyzSplit
        args_x: dict = {"file": splitIn, "atoms": [key, value], "cuts": args.cuts,
                        "print": False, "out": splitOut}
        # sys.stdout = open(os.devnull, 'w')
        xyzSplit.main(argparse.Namespace(**args_x))
        # sys.stdout = sys.__stdout__
        move_file(splitOut, splitIn)

    move_file(splitIn, outFile)
    print(f" The data is saved to {outFile} !!!")


def main(args: argparse.Namespace = argparse.Namespace()) -> None:

    if args == argparse.Namespace():
        args = cml()
    print_arguments()

    _, circleMols, residualMols, Bond_order, atomsCN, residualMols_all_pairs = read_data(
        args)

    if args.verbose:
        ic(circleMols)

    flattenCircleMols: list[int] = []
    for mol in circleMols:
        flattenCircleMols += mol
    flattenCircleMols = list(set(flattenCircleMols))
    if args.verbose:
        ic(residualMols, flattenCircleMols)
    xyzSplit: dict[int, int] = get_xyzSplit(residualMols,
                                            Bond_order, atomsCN, flattenCircleMols, residualMols_all_pairs)
    gen_GeometryXYZs(xyzSplit, args)


if __name__ == "__main__":
    main()

    # test
    # python3 xyzGenFlexible.py -i tests/data/crest_conformers.xyz
    # python3 xyzGenFlexible.py -i tests/data/crest_conformers.xyz -m
