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


def read_data(_xyzFile: GeometryXYZs, _verbose: bool) -> tuple[dict[AtomID, npt.NDArray[np.int64]], list[list[AtomID]], list[set[int]], dict[AtomID, int], dict[AtomID, int], dict]:
    from censo_ext.Tools.topo import Topo
    from censo_ext.Tools.ml4nmr import read_mol_neighbors_bond_order

    Sts_topo: Topo = Topo(xyzFile=_xyzFile)

    neighbor, circleMols, residualMols, residualMols_all_pairs = Sts_topo.topology()
    idx_atomsCN: dict[AtomID, int] = Sts_topo.get_cn()
    if _verbose:
        ic(neighbor, circleMols, residualMols)
        ic(idx_atomsCN)
        # xyzFile: GeometryXYZs = GeometryXYZs(_file)
        # xyzFile.method_read_xyz()
    *_, idx_Bond_order = read_mol_neighbors_bond_order(xyzFile=_xyzFile)
    if _verbose:
        ic(idx_Bond_order)
        ic(residualMols)
    return neighbor, circleMols, residualMols, idx_Bond_order, idx_atomsCN, residualMols_all_pairs


def get_xyzSplit(residualMols: list[set[int]], atomsCN: dict[AtomID, int], flattenCircleMols: list[int], residualMols_all_pairs) -> dict[int, int]:
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

        flexibleMolsCNis4: list = [
            a for a in flexibleMols if atomsCN[AtomID(a)] == 4]
        # ic(mol, flexibleMolsCNis4, nodeMols)
        # ic(nodeMols)

        if len(nodeMols) == 0:

            max_distnce: int = 0
            atomIDs: int = 0
            for x in residualMols_all_pairs.values():
                for key, distance in x.items():
                    if distance > max_distnce:
                        max_distnce = distance
                        atomIDs = key

            # Neighbor_atomIDs = [
            #    key for key, value in residualMols_all_pairs[atomIDs].items() if value == 1]

            nodeMols.append(atomIDs)
            # temp_Mols = Neighbor_atomIDs[0]
            flexibleMols.remove(atomIDs)
            # ic(nodeMols)
            # ic(flexibleMols)
            # for key, value in residualMols_all_pairs.items():
            #    for x, y in value.items():
            #        if x == temp_Mols:
            #            value[x] = 0
            # ic(residualMols_all_pairs)

        # ic(flexibleMols, nodeMols)

        if len(nodeMols) == 1 or 2:
            # ic(residualMols_all_pairs[nodeMols[0]])
            # ic(flexibleMolsCNis4)
            a = residualMols_all_pairs[nodeMols[0]].values()
            import math
            b = [x for x in a if not math.isinf(x)]
            for x in range(0, max(b)-1):
                # ic(x, x+1)
                out_key: int = 0
                out_value: int = 0
                if x == 0:
                    # ic("key: ", nodeMols[0])
                    out_key = nodeMols[0]
                for key, distance in residualMols_all_pairs[nodeMols[0]].items():
                    if distance == x:
                        if len([c for c in residualMols_all_pairs[key].values() if c == 1]) != 1:
                            # ic("key: ", key)
                            out_key = key
                    if distance == x+1:
                        if len([c for c in residualMols_all_pairs[key].values() if c == 1]) != 1:
                            # ic("value: ", key)
                            out_value = key
                if out_key not in flexibleMolsCNis4 and out_value not in flexibleMolsCNis4:
                    pass
                # if out_key == 0 or out_value == 0:
                #    pass
                else:
                    # ic(out_key, out_value)
                    xyzSplit[out_key] = out_value

        else:
            print("len(nodeMols)= ", len(nodeMols))
            raise NotImplementedError("Under Construct")
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

    xyzFile: GeometryXYZs = GeometryXYZs(args.file)
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

    move_file(splitIn, args.out)
    print(f" The data is saved to {args.out} !!!")


def main(args: argparse.Namespace = argparse.Namespace()) -> None:

    if args == argparse.Namespace():
        args = cml()
    print_arguments()

    xyzFile = GeometryXYZs(args.file)
    xyzFile.method_read_xyz()
    _, circleMols, residualMols, Bond_order, atomsCN, residualMols_all_pairs = read_data(
        _xyzFile=xyzFile, _verbose=args.verbose)

    if args.verbose:
        ic(circleMols)

    flattenCircleMols: list[int] = []
    for mol in circleMols:
        flattenCircleMols += mol
    flattenCircleMols = list(set(flattenCircleMols))
    if args.verbose:
        ic(residualMols, flattenCircleMols)
    xyzSplit: dict[int, int] = get_xyzSplit(
        residualMols, atomsCN, flattenCircleMols, residualMols_all_pairs)
    gen_GeometryXYZs(xyzSplit, args)


if __name__ == "__main__":
    main()

    # test
    # python3 xyzGenFlexible.py -i tests/data/crest_conformers.xyz
    # python3 xyzGenFlexible.py -i tests/data/crest_conformers.xyz -m
