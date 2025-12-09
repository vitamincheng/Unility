#!/usr/bin/env python
import argparse
# import os
from scipy.spatial.transform import Rotation as R
import numpy as np
import numpy.typing as npt
from censo_ext.Tools.utility import AtomID, print_arguments
from censo_ext.Tools.xyzfile import GeometryXYZs
# from icecream import ic
from pathlib import Path

descr = """
________________________________________________________________________________
| For Return to origin and lay on the xz plane
| Usage    : xyzReturnOandZ.py <geometry> [options]                  
| Input    : -i xyz file [default traj.xyz]
| Output   : -o Saved xyz file [default output.xyz] 
| [Options]
| Atom     : -a or --atom [1 2 3] idx of atom 
|              1 : Fixed atom and return origin
|              2 : Rotation atom and z axis
|              3 : Rotation atom and lay on xz plane
| Automatic: --auto Automatically search for the minimum deviation setting 
|                   origin of all atoms [default False] 
| Replace  : -r Replace the input file [default False]
| Print    : -p Print the final data on screen [default False]
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
        "-a",
        "--atom",
        dest="atom",
        action="store",
        type=int,
        nargs=3,
        help="Provide three idx of atom's nubmers [1(origin) 2(z axis) 3(xz plane)]",
    )

    parser.add_argument(
        "-r",
        "--replace",
        dest="replace",
        action="store_true",
        help="Replace the original input file [default False]",
    )

    parser.add_argument(
        "-p",
        "--print",
        dest="print",
        action="store_true",
        help="Print the final data on stdout (on screen) [default False]",
    )

    parser.add_argument(
        "--auto",
        dest="auto",
        action="store_true",
        help="Automatically search for the minimum deviation setting origin of all atoms [default False]",
    )

    args: argparse.Namespace = parser.parse_args()
    return args


def idx_3atom_opt(inFile: Path) -> tuple[AtomID, AtomID, AtomID]:
    from censo_ext.Tools.factor import method_factor_analysis
    args_x: dict = {"file": inFile,
                    "factor": 0.5, "debug": False, "opt": False}
    _LowFactor: list[AtomID]
    _Deviation: dict[AtomID, float]
    _LowFactor, _Deviation = method_factor_analysis(
        args=argparse.Namespace(**args_x))

    _Bonding: list[list[AtomID]] = []
    for x in _LowFactor:
        from censo_ext.Tools.topo import Topo
        # Sts_topo: Topo = Topo(inFile)
        _Bonding.append(Topo(inFile).method_bonding(_bonding=x, _print=False))

    _3AtomID: list[list[AtomID]] = []
    for idx0, x in enumerate(_LowFactor):
        # total numbers >=3 or >2 (one of total numbers is )
        if len(_Bonding[idx0]) > 1:
            tmp: list[AtomID] = []
            tmp.append(x)
            for y in _Bonding[idx0]:
                tmp.append(y)
            _3AtomID.append(tmp)

    from itertools import combinations
    Combined_3AtomID: list[tuple[AtomID, AtomID, AtomID]] = []
    for x in _3AtomID:
        for y in list(combinations(x, 3)):
            Combined_3AtomID.append(y)

    idx1_Atoms: list[AtomID] = list(_Deviation.keys())
    STD_Atoms: list[float] = list(_Deviation.values())

    intp_minTotalDev: int = 0
    minTotalDev: float = 100
    for idx0, x in enumerate(Combined_3AtomID):
        TotalDevAtoms: float = 0.0
        for y in x:
            TotalDevAtoms += (STD_Atoms[idx1_Atoms.index(AtomID(y))])
        if minTotalDev > TotalDevAtoms:
            minTotalDev = TotalDevAtoms
            intp_minTotalDev: int = idx0

    print("")
    print(f" 3 atom idx of lowest total factor {Combined_3AtomID[intp_minTotalDev]}")  # nopep8
    print("")
    return (Combined_3AtomID[intp_minTotalDev])


def main(args: argparse.Namespace = argparse.Namespace()) -> None:

    if args == argparse.Namespace():
        args = cml()
    print_arguments()

    from censo_ext.Tools.utility import IsExist
    inFile = Path(args.file)
    outFile = Path(args.out)
    IsExist(args.file)

    if not args.atom and not args.auto:
        raise ValueError(" No any sepific atom in your provided arguments ")

    p_idx1: AtomID
    q_idx1: AtomID
    r_idx1: AtomID
    if not args.atom and args.auto:
        print("\n Automated to set the 3 atoms to return origin and lay on XZ plane")
        print(" First FactorAnalysis.py will executive and second continue the RetrunOandZ.py ")
        p_idx1, q_idx1, r_idx1 = idx_3atom_opt(inFile)
    else:
        p_idx1, q_idx1, r_idx1 = args.atom

    xyzFile: GeometryXYZs = GeometryXYZs(inFile)
    xyzFile.method_read_xyz()

    # Process xyz file
    for idx0_St in range(len(xyzFile)):

        dxyz: npt.NDArray[np.float64] = xyzFile.Sts[idx0_St].coord[p_idx1-1].copy()
        xyzFile.Sts[idx0_St].coord -= dxyz  # type: ignore
        z_axis = (0, 0, np.sqrt(
            np.sum(np.square(xyzFile.Sts[idx0_St].coord[q_idx1-1]))))

        rotation_axis = xyzFile.Sts[idx0_St].coord[q_idx1-1] + z_axis

        Normalized_RotationAxis: npt.NDArray[np.float64] = np.array([
            0, 1, 0]) if np.linalg.norm(rotation_axis) == 0 else rotation_axis / np.linalg.norm(rotation_axis)

        R_pq = R.from_rotvec(np.pi*Normalized_RotationAxis)
        xyzFile.Sts[idx0_St].coord = R_pq.apply(
            xyzFile.Sts[idx0_St].coord)  # type: ignore

        Angle_qr = np.angle(complex(xyzFile.Sts[idx0_St].coord[r_idx1-1][0], complex(
            xyzFile.Sts[idx0_St].coord[r_idx1-1][1])))
        R_qr = R.from_euler('z', -Angle_qr)
        xyzFile.Sts[idx0_St].coord = R_qr.apply(
            xyzFile.Sts[idx0_St].coord)  # type: ignore

    # Save or print result
    if args.print:
        xyzFile.method_print([])
    else:
        fileName: Path = inFile if args.replace else outFile
        print(f"    Saved to {fileName}")
        xyzFile.set_filename(fileName)
        xyzFile.method_save_xyz([])


if __name__ == "__main__":
    main()

    # test
    # python3 xyzReturnOandZ.py -i ../tests/crest_conformers.xyz -a 30 45 47
    # python3 xyzReturnOandZ.py -i ../tests/crest_conformers.xyz --auto
