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
        help="Provide index of xyz file (-1 is all indexes) [default 1]",
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
        "--auto",
        dest="auto",
        action="store_true",
        help="Auto mode of saved files by use xyzReturnOandZ [default False]",
    )

    args: argparse.Namespace = parser.parse_args()
    return args


def cal_RMSD(_xyzFile: GeometryXYZs, idx_p: int, idx_q: int, bond_broken: tuple[int, int]) -> float:

    from censo_ext.Tools.calculate_rmsd import cal_RMSD_xyz
    _, RMSD = cal_RMSD_xyz(
        _xyzFile, idx_p, idx_q, _remove_idx=None, _add_idx=[bond_broken[1]], _bond_broken=bond_broken, _ignore_Hydrogen=True)
    return RMSD


def TopoAnalysis(_xyzFile: GeometryXYZs, _index: int, _verbose: bool, _limits: float) -> tuple[list[cell_reports], list[cell_reports]]:

    idx1_p: int = _index
    tmp: bool = _verbose
    neighbor, circleMols, residualMols, Bond_order, atomsCN, residualMols_all_pairs = read_data(
        _xyzFile=_xyzFile, _verbose=False)
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

    if len(_xyzFile) < idx1_p:
        print(
            f" Error: Index {idx1_p} is out of range in your xyz file{_xyzFile.get_fileName()}.")
        exit(0)

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

            for x in range(1, len(_xyzFile)+1):
                if x == idx1_p:
                    continue
                std_left: float = cal_RMSD(_xyzFile=_xyzFile, idx_p=idx1_p, idx_q=x,
                                           bond_broken=(node_mol, res_node_mol))
                std_right: float = cal_RMSD(_xyzFile=_xyzFile, idx_p=idx1_p, idx_q=x,
                                            bond_broken=(res_node_mol, node_mol))
                if std_left <= _limits and std_right <= _limits:
                    result_circle.append(
                        (x, node_mol, res_node_mol, std_left, std_right))

    # Check straight chain
    result_straight: list[cell_reports] = []
    for key, value in xyzSplit.items():
        for x in range(1, len(_xyzFile)+1):
            if x == idx1_p:
                continue
            std_left = cal_RMSD(_xyzFile=_xyzFile, idx_p=idx1_p, idx_q=x,
                                bond_broken=(key, value))
            std_right = cal_RMSD(_xyzFile=_xyzFile, idx_p=idx1_p, idx_q=x,
                                 bond_broken=(value, key))
            if std_left <= _limits and std_right <= _limits:
                if _verbose:
                    ic(x, key, value, std_left, std_right)
                result_straight.append((x, key, value, std_left, std_right))

    return result_circle, result_straight


def print_report(_index: int, _circle: list[cell_reports], _straight: list[cell_reports]) -> None:

    if len(_circle) != 0:
        print("")
        print("  ===== Check circle molecule =====")
        print("   idx1_p   idx1_q  node  res_node      res_left      res_right")
        for x in _circle:
            print(
                f"   {_index:6d}   {x[0]:6d} {x[1]:6d} {x[2]:8d} {x[3]:14.7f} {x[4]:14.7f}")

    if len(_straight) != 0:
        print("")
        print("  ===== Check straight molecule =====")
        print("   idx1_p   idx1_q    key    value       res_left      res_right")
        for x in _straight:
            print(
                f"   {_index:6d}   {x[0]:6d} {x[1]:6d} {x[2]:8d} {x[3]:14.7f} {x[4]:14.7f}")
        print("  [key,value] [fixed,rotation]")


def save_files(_index: int, _xyzFile: GeometryXYZs, _circle: list[cell_reports], _straight: list[cell_reports], _auto: bool) -> None:

    if len(_circle) >= 1:
        print("")
        print("  ===== Save circle molecule =====")
        # print(result_circle)
        circleDir: Path = Path("Circle")
        if circleDir.is_dir():
            import shutil
            shutil.rmtree(circleDir, ignore_errors=True)
        circleDir.mkdir()

        pairs: set[tuple[int, int]] = {(x[1], x[2]) for x in _circle}
        for x in pairs:
            index1: list[int] = [_index]
            for y in _circle:
                if x == (y[1], y[2]):
                    index1.append(y[0])
            print(index1)
            outFile: Path = Path('_'.join(str(x) for x in index1)+".xyz")
            _xyzFile.set_filename(circleDir / outFile)
            if _auto:
                _xyzFile.method_xyzReturnOandZ_auto()
            _xyzFile.method_save_xyz(index1)

    if len(_straight) >= 1:
        print("")
        print("  ===== Save straight molecule =====")
        straightDir: Path = Path("Straight")
        if straightDir.is_dir():
            import shutil
            shutil.rmtree(straightDir, ignore_errors=True)
        straightDir.mkdir()

        pairs = {(x[1], x[2]) for x in _straight}
        for x in pairs:
            index1: list[int] = [_index]
            for y in _straight:
                if x == (y[1], y[2]):
                    index1.append(y[0])
            print(index1)
            outFile: Path = Path('_'.join(str(x) for x in index1)+".xyz")
            _xyzFile.set_filename(straightDir / outFile)
            if _auto:
                _xyzFile.method_xyzReturnOandZ_auto()
            _xyzFile.method_save_xyz(index1)

    print("  ===== Finished to save the files =====")


def main(args: argparse.Namespace = argparse.Namespace()) -> None:
    if args == argparse.Namespace():
        args = cml()
    print_arguments()

    _xyzFile: GeometryXYZs = GeometryXYZs(args.file)
    _xyzFile.method_read_xyz()

    print("  ===== Parameter of limits =====")
    print(f"  The total numbers of the xyzfile       = {len(_xyzFile.Sts)}")
    print(f"  The delta limits of standard deviation = {args.limits}")

    if args.idx == -1:
        for x in range(1, len(_xyzFile.Sts)+1):
            _Circle, _Straight = TopoAnalysis(_xyzFile=_xyzFile, _index=x,
                                              _verbose=args.verbose, _limits=args.limits)
            print_report(_index=x, _circle=_Circle, _straight=_Straight)
    else:
        _Circle, _Straight = TopoAnalysis(_xyzFile=_xyzFile, _index=args.idx,
                                          _verbose=args.verbose, _limits=args.limits)
        print_report(_index=args.idx, _circle=_Circle, _straight=_Straight)
        save_files(_index=args.idx, _xyzFile=_xyzFile, _circle=_Circle,
                   _straight=_Straight, _auto=args.auto)


if __name__ == "__main__":
    main()
