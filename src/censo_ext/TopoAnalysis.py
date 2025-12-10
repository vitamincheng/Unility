#!/usr/bin/env python
import argparse
from pathlib import Path
from icecream import ic
from censo_ext.Tools.utility import AtomID, print_arguments
from censo_ext.Tools.xyzfile import GeometryXYZs
from censo_ext.xyzGenFlexible import get_xyzSplit, read_data
type cell_reports = tuple[int, int, int, float, float]

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

    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        help="Verbose mode [default False] ",
    )

    parser.add_argument(
        "-c",
        "--check",
        dest="check",
        action="store_true",
        help="Check mode of chemical structures [default False] ",
    )

    args: argparse.Namespace = parser.parse_args()
    return args


def cal_RMSD(xyzfile: GeometryXYZs, idx_p: int, idx_q: int, bond_broken: tuple[int, int], check: bool) -> float:

    from censo_ext.Tools.calculate_rmsd import cal_RMSD_xyz
    _, RMSD = cal_RMSD_xyz(
        xyzfile, idx_p, idx_q, _remove_idx=None, _add_idx=[bond_broken[1]], _bond_broken=bond_broken, _ignore_Hydrogen=True, _check=check)
    return RMSD


def TopoAnalysis(_file: Path, _index: int, _verbose: bool, _limits: float, _check: bool) -> tuple[list[cell_reports], list[cell_reports]]:

    idx1_p: int = _index
    tmp: bool = _verbose
    neighbor, circleMols, residualMols, Bond_order, atomsCN, residualMols_all_pairs = read_data(
        _file=_file, _verbose=False, _check=_check)
    _verbose = tmp

    flattenCircleMols: list[int] = []
    for mol in circleMols:
        flattenCircleMols += mol
    flattenCircleMols = list(set(flattenCircleMols))
    xyzSplit: dict[int, int] = get_xyzSplit(
        residualMols, atomsCN, flattenCircleMols, residualMols_all_pairs)

    if _verbose:
        ic(xyzSplit)
        ic(residualMols)
        ic(circleMols)
        ic(neighbor)
        ic(flattenCircleMols)
        ic(residualMols_all_pairs)

    xyzFile: GeometryXYZs = GeometryXYZs(_file)
    xyzFile.method_read_xyz()
    limits = _limits
    print("  ===== Parameter of limits =====")
    print(f"  the delta limits of standard deviation = {limits}")

    # Check circleMols factor
    result_circle: list[cell_reports] = []
    for resMol in residualMols:
        if _verbose:
            ic(resMol)
        node_mols: list[int] = [
            int(x) for x in resMol if x in flattenCircleMols]
        for node_mol in node_mols:
            inter_x: set[int] = set.intersection(
                set(map(int, neighbor[AtomID(int(node_mol))])), resMol)
            if (len(inter_x)) != 1:
                print("  Something wrong in your residualMols")
                print("  Exit and Close the program !!!")
                exit(1)
            res_node_mol: int = list(inter_x)[0]
            if _verbose:
                ic(node_mol, res_node_mol)

            for x in range(1, len(xyzFile)+1):
                if x == idx1_p:
                    continue
                res_left: float = cal_RMSD(xyzfile=xyzFile, idx_p=idx1_p, idx_q=x,
                                           bond_broken=(node_mol, res_node_mol), check=_check)
                res_right: float = cal_RMSD(xyzfile=xyzFile, idx_p=idx1_p, idx_q=x,
                                            bond_broken=(res_node_mol, node_mol), check=_check)
                if res_left <= limits and res_right <= limits:
                    result_circle.append(
                        (x, node_mol, res_node_mol, res_left, res_right))

    # Check straight chain
    result_straight: list[cell_reports] = []
    for key, value in xyzSplit.items():
        for x in range(1, len(xyzFile)+1):
            if x == idx1_p:
                continue
            res_left = cal_RMSD(xyzfile=xyzFile, idx_p=idx1_p, idx_q=x,
                                bond_broken=(key, value), check=_check)
            res_right = cal_RMSD(xyzfile=xyzFile, idx_p=idx1_p, idx_q=x,
                                 bond_broken=(value, key), check=_check)
            if res_left <= limits and res_right <= limits:
                if _verbose:
                    ic(x, key, value, res_left, res_right)
                result_straight.append((x, key, value, res_left, res_right))

    return result_circle, result_straight


def print_report(result_circle: list[cell_reports], result_straight: list[cell_reports]) -> None:

    print("  ===== Check circle molecule =====")
    print("   idx1  node  res_node      res_left      res_right")
    for x in result_circle:
        print(f"{x[0]:6d} {x[1]:6d} {x[2]:8d} {x[3]:14.7f} {x[4]:14.7f}")

    print("  ===== Check straight molecule =====")
    print("   idx1   key     value      res_left      res_right")
    for x in result_straight:
        print(f"{x[0]:6d} {x[1]:6d} {x[2]:8d} {x[3]:14.7f} {x[4]:14.7f}")


def save_files(_index: int, _file: Path, result_circle: list[cell_reports], result_straight: list[cell_reports]) -> None:

    _file = Path(_file)

    if len(result_circle) >= 1:
        print("  ===== Save circle molecule =====")
        # print(result_circle)
        circleDir: Path = Path("Circle")
        if circleDir.is_dir():
            import shutil
            shutil.rmtree(circleDir, ignore_errors=True)
        circleDir.mkdir()

        pairs: set[tuple[int, int]] = {(x[1], x[2]) for x in result_circle}
        for x in pairs:
            index1: list[int] = [_index]
            for y in result_circle:
                if x == (y[1], y[2]):
                    index1.append(y[0])
            print(index1)
            inFile: Path = Path(_file)
            outFile: Path = Path('_'.join(str(x) for x in index1)+".xyz")
            xyzFile: GeometryXYZs = GeometryXYZs(inFile)
            xyzFile.method_read_xyz()
            xyzFile.set_filename(circleDir / outFile)
            xyzFile.method_save_xyz(index1)

    if len(result_straight) >= 1:
        print("  ===== Save straight molecule =====")
        # print(result_straight)
        straightDir: Path = Path("Straight")
        if straightDir.is_dir():
            import shutil
            shutil.rmtree(straightDir, ignore_errors=True)
        straightDir.mkdir()

        pairs = {(x[1], x[2]) for x in result_straight}
        for x in pairs:
            index1: list[int] = [_index]
            for y in result_straight:
                if x == (y[1], y[2]):
                    index1.append(y[0])
            print(index1)
            inFile: Path = Path(_file)
            outFile: Path = Path('_'.join(str(x) for x in index1)+".xyz")
            xyzFile: GeometryXYZs = GeometryXYZs(inFile)
            xyzFile.method_read_xyz()
            xyzFile.set_filename(straightDir / outFile)
            xyzFile.method_save_xyz(index1)

    print("  ===== Finished to save the files =====")


def main(args: argparse.Namespace = argparse.Namespace()) -> None:
    if args == argparse.Namespace():
        args = cml()
    print_arguments()

    result_circle, result_straight = TopoAnalysis(_file=args.file, _index=args.idx,
                                                  _verbose=args.verbose, _limits=args.limits, _check=args.check)
    print_report(result_circle=result_circle, result_straight=result_straight)
    save_files(_index=args.idx, _file=args.file, result_circle=result_circle,
               result_straight=result_straight)


if __name__ == "__main__":
    main()
