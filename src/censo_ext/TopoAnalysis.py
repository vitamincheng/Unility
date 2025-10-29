#!/usr/bin/env python
import argparse
from censo_ext.Tools.utility import print_arguments
from censo_ext.Tools.xyzfile import GeometryXYZs
from censo_ext.xyzGenFlexible import get_xyzSplit, read_data
descr = """
________________________________________________________________________________
| Usages   : TopoAnalysis.py <geometry> [options]
| Input    : -i input xyz file [default traj.xyz]
| [options]
| index     : -d index of reference structure in xyz file [defalut 1]
| limits    : -l limits of delta std in xyz file [default 0.10]
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
        type=str,
        default="traj.xyz",
        help="Provide input xyz file name [default traj.xyz]",
    )

    parser.add_argument(
        "-d",
        "--idx",
        dest="idx",
        action="store",
        required=False,
        type=int,
        default=1,
        help="Provide index of xyz file [default 1] ",
    )

    parser.add_argument(
        "-l",
        "--limits",
        dest="limits",
        action="store",
        required=False,
        type=float,
        default=0.1,
        help="Provide limits of delta std in xyz file [default 0.1] ",
    )

    args: argparse.Namespace = parser.parse_args()
    return args


def cal_RMSD(xyzfile, idx_p, idx_q, bond_broken) -> float:
    from censo_ext.Tools.calculate_rmsd import cal_RMSD_xyz
    x: dict = {"remove_idx": None, "add_idx": None,
               "bond_broken": bond_broken, "ignore_Hydrogen": True}
    _, RMSD = cal_RMSD_xyz(
        xyzfile, idx_p, idx_q, args=argparse.Namespace(**x))
    return RMSD


def TopoAnalysis(args) -> None:

    idx_p = args.idx
    args.verbose = False
    neighbor, circleMols, residualMols, Bond_order, atomsCN, residualMols_all_pairs = read_data(
        args)
    flattenCircleMols: list[int] = []
    for mol in circleMols:
        flattenCircleMols += mol
    flattenCircleMols = list(set(flattenCircleMols))

    xyzSplit: dict[int, int] = get_xyzSplit(residualMols,
                                            Bond_order, atomsCN, flattenCircleMols, residualMols_all_pairs)
    from icecream import ic
    ic(residualMols)
    ic(circleMols)
    ic(neighbor)
    ic(flattenCircleMols)
    # ic(residualMols_all_pairs)
    xyzFile: GeometryXYZs = GeometryXYZs(args.file)
    xyzFile.method_read_xyz()
    limits = args.limits
    print("  ===== Parameter of limits =====")
    print(f"  the delta limits of standard deviation = {limits}")

    # Check circleMols factor
    print("  ===== Check circle molecule =====")
    result_circle: list = []
    for resMol in residualMols:
        # ic(resMol)
        node_mols = [x for x in resMol if x in flattenCircleMols]
        for node_mol in node_mols:
            # ic(node_mol)
            x = set.intersection(
                set(map(int, neighbor[int(node_mol)])), resMol)
            if (len(x)) != 1:
                print("  Something wrong in your residualMols")
                print("  Exit and Close the program !!!")
                exit(1)
            res_node_mol = list(x)[0]
            # ic(node_mol, res_node_mol)
            nNums = len(xyzFile)
            for x in range(1, nNums+1):
                if x == idx_p:
                    continue
                res_left = cal_RMSD(xyzfile=xyzFile, idx_p=idx_p, idx_q=x,
                                    bond_broken=(node_mol, res_node_mol))
                res_right = cal_RMSD(xyzfile=xyzFile, idx_p=idx_p, idx_q=x,
                                     bond_broken=(res_node_mol, node_mol))
                if res_left <= limits and res_right <= limits:
                    result_circle.append(
                        (x, node_mol, res_node_mol, res_left, res_right))
    print("   idx1  node  res_node      res_left      res_right")
    for x in result_circle:
        print(f"{x[0]:6d} {x[1]:6d} {x[2]:8d} {x[3]:14.7f} {x[4]:14.7f}")
    # Check straight chain
    print("  ===== Check straight molecule =====")
    result_straight: list = []
    for key, value in xyzSplit.items():
        nNums = len(xyzFile)
        for x in range(1, nNums+1):
            if x == idx_p:
                continue
            res_left = cal_RMSD(xyzfile=xyzFile, idx_p=idx_p, idx_q=x,
                                bond_broken=(key, value))
            res_right = cal_RMSD(xyzfile=xyzFile, idx_p=idx_p, idx_q=x,
                                 bond_broken=(value, key))
            if res_left <= limits and res_right <= limits:
                # ic(x, key, value, res_left, res_right)
                result_straight.append((x, key, value, res_left, res_right))
    print("   idx1   key     value      res_left      res_right")
    for x in result_straight:
        print(f"{x[0]:6d} {x[1]:6d} {x[2]:8d} {x[3]:14.7f} {x[4]:14.7f}")


def main(args: argparse.Namespace = argparse.Namespace()) -> None:
    if args == argparse.Namespace():
        args = cml()
    print_arguments()

    TopoAnalysis(args)


if __name__ == "__main__":
    main()
