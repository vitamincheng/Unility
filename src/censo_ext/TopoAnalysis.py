#!/usr/bin/env python
import argparse
from pathlib import Path
from icecream import ic
from censo_ext.Tools.utility import AtomID, print_arguments
from censo_ext.Tools.xyzfile import GeometryXYZs
from censo_ext.xyzGenFlexible import get_xyzSplit, read_data
type cell_reports = tuple[int, int, int, int, float, float]

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


def cal_RMSD(_xyzFile: GeometryXYZs, idx1_p: int, idx1_q: int, bond_broken: tuple[int, int]) -> float:

    from censo_ext.Tools.calculate_rmsd import cal_RMSD_xyz
    _, _RMSD = cal_RMSD_xyz(xyzFile=_xyzFile, idx1_p=idx1_p, idx1_q=idx1_q,
                            _remove_idx=None, _add_idx=[bond_broken[1]], _bond_broken=bond_broken, _ignore_Hydrogen=True)
    return _RMSD


def TopoAnalysis(_xyzFile: GeometryXYZs, _idx1: int, _verbose: bool, _limits: float) -> tuple[list[cell_reports], list[cell_reports]]:

    idx1_p: int = _idx1
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

            for idx1_q in range(1, len(_xyzFile)+1):
                if idx1_p == idx1_q:
                    continue
                std_left: float = cal_RMSD(_xyzFile=_xyzFile, idx1_p=idx1_p, idx1_q=idx1_q,
                                           bond_broken=(node_mol, res_node_mol))
                std_right: float = cal_RMSD(_xyzFile=_xyzFile, idx1_p=idx1_p, idx1_q=idx1_q,
                                            bond_broken=(res_node_mol, node_mol))
                if std_left <= _limits and std_right <= _limits:
                    result_circle.append(
                        (idx1_p, idx1_q, node_mol, res_node_mol, std_left, std_right))

    # Check straight chain
    result_straight: list[cell_reports] = []
    for key, value in xyzSplit.items():
        for idx1_q in range(1, len(_xyzFile)+1):
            if idx1_p == idx1_q:
                continue
            std_left = cal_RMSD(_xyzFile=_xyzFile, idx1_p=idx1_p, idx1_q=idx1_q,
                                bond_broken=(key, value))
            std_right = cal_RMSD(_xyzFile=_xyzFile, idx1_p=idx1_p, idx1_q=idx1_q,
                                 bond_broken=(value, key))
            if std_left <= _limits and std_right <= _limits:
                result_straight.append(
                    (idx1_p, idx1_q, key, value, std_left, std_right))

    return result_circle, result_straight


def print_report(_xyzFile: GeometryXYZs, _circle: list[cell_reports], _straight: list[cell_reports]) -> None:

    if len(_circle) != 0:
        print("")
        print("  ===== Check circle molecule =====")
        print(
            " idx1_p idx1_q    #_p    #_q     node res_node       res_left      res_right")
        for idx1_p, idx1_q, node, res_node, res_left, res_right in _circle:
            nClusters_p = _xyzFile.Sts[idx1_p-1].comment_nClusters
            nClusters_q = _xyzFile.Sts[idx1_q-1].comment_nClusters
            print(
                f"  {idx1_p:5d}  {idx1_q:5d}  {nClusters_p:5d}  {nClusters_q:5d}    {node:5d}    {res_node:5d} {res_left:14.7f} {res_right:14.7f}")

    if len(_straight) != 0:
        print("")
        print("  ===== Check straight molecule =====")
        print(
            " idx1_p idx1_q    #_p    #_q      key    value       res_left      res_right")
        for idx1_p, idx1_q, key, value, res_left, res_right in _straight:
            nClusters_p = _xyzFile.Sts[idx1_p-1].comment_nClusters
            nClusters_q = _xyzFile.Sts[idx1_q-1].comment_nClusters
            print(
                f"  {idx1_p:5d}  {idx1_q:5d}  {nClusters_p:5d}  {nClusters_q:5d}    {key:5d}    {value:5d} {res_left:14.7f} {res_right:14.7f}")
        print("  [key,value] [fixed,rotation]")


def save_files(_xyzFile: GeometryXYZs, _circle: list[cell_reports], _straight: list[cell_reports], _auto: bool) -> None:

    if len(_circle) >= 1:
        print("")
        print("  ===== Save circle molecule =====")
        # print(result_circle)
        circleDir: Path = Path("Circle")
        if circleDir.is_dir():
            import shutil
            shutil.rmtree(circleDir, ignore_errors=True)
        circleDir.mkdir()

        pairs: set[tuple[int, int]] = {(x[2], x[3]) for x in _circle}
        for x in pairs:
            index1: list[int] = [_circle[0][0]]
            nClusters: list[int] = [_circle[0][0]]
            for y in _circle:
                if x == (y[2], y[3]):
                    index1.append(y[1])
                    nClusters.append(_xyzFile.Sts[y[1]-1].comment_nClusters)

            print(Path('_'.join(str(x)
                  for x in index1)+".xyz"), "\t\t\t", nClusters)

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

        pairs = {(x[2], x[3]) for x in _straight}
        for x in pairs:
            index1: list[int] = [_straight[0][0]]
            nClusters: list[int] = [_straight[0][0]]
            for y in _straight:
                if x == (y[2], y[3]):
                    index1.append(y[1])
                    nClusters.append(_xyzFile.Sts[y[1]-1].comment_nClusters)

            print(Path('_'.join(str(x)
                  for x in index1)+".xyz"), "\t\t\t", nClusters)

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
    print(f"  The total numbers of the xyzfile         = {len(_xyzFile.Sts)}")
    print(f"  The delta limits of standard deviation   = {args.limits}")
    print(f"  The index1 of the structures in xyz file = {args.idx}")

    if args.idx == -1:
        for idx1 in range(1, len(_xyzFile.Sts)+1):
            _Circle, _Straight = TopoAnalysis(_xyzFile=_xyzFile, _idx1=idx1,
                                              _verbose=args.verbose, _limits=args.limits)
            print_report(_xyzFile=_xyzFile, _circle=_Circle,
                         _straight=_Straight)
    else:
        _Circle, _Straight = TopoAnalysis(_xyzFile=_xyzFile, _idx1=args.idx,
                                          _verbose=args.verbose, _limits=args.limits)
        print_report(_xyzFile=_xyzFile, _circle=_Circle, _straight=_Straight)
        save_files(_xyzFile=_xyzFile, _circle=_Circle,
                   _straight=_Straight, _auto=args.auto)


if __name__ == "__main__":
    main()
