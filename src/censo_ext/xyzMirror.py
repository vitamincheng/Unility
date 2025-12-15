#!/usr/bin/env python
import argparse
from scipy.spatial.transform import Rotation as R
import numpy as np
import numpy.typing as npt
from censo_ext.Tools.utility import AtomID, IntpID, print_arguments
from censo_ext.Tools.xyzfile import GeometryXYZs
from pathlib import Path

descr = """
________________________________________________________________________________
| For Return to origin and mirror on the xz plane
| Usage    : xyzMirror.py <geometry> [options]                  
| Input    : -i xyz file [default traj.xyz]
| Output   : -o Saved xyz file [default output.xyz] 
| [Options]
| Atom     : -a or --atom [1 2 3] idx of atom [required] 
|              1 : Fixed atom and return origin (mirror)
|              2 : Rotation atom and z axis (not mirror)
|              3 : Rotation atom and lay on xz plane (not mirror)
| Replace  : -r Replace the input file [default False]
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
        required=True,
        help="Provide three idx of atom's nubmers (mirror xz plane, active 1 but not include 2,3) \
            if empty will mirror total atoms [1(origin) 2(z axis) 3(xz plane)]",
    )

    parser.add_argument(
        "-r",
        "--replace",
        dest="replace",
        action="store_true",
        help="Replace the original input file [default False]",
    )

    args: argparse.Namespace = parser.parse_args()
    return args


def main(args: argparse.Namespace = argparse.Namespace()) -> None:

    if args == argparse.Namespace():
        args = cml()
    print_arguments()

    from censo_ext.Tools.utility import IsExist
    inFile = Path(args.file)
    outFile = Path(args.out)
    IsExist(inFile)

    if not args.atom:
        print("  No any sepific atom in your provided arguments ")
        print("  xyzMirror.py will mirror all atoms !!!")
        xyzFile: GeometryXYZs = GeometryXYZs(inFile)
        xyzFile.method_read_xyz()
        Nums: int = len(xyzFile.Sts[0].coord)
        idx0_H: list[IntpID] = [IntpID(x) for x in [*range(Nums)]]
        p_idx1: AtomID = AtomID(1)
        q_idx1: AtomID = AtomID(1)
        r_idx1: AtomID = AtomID(1)
    else:
        p_idx1: AtomID
        q_idx1: AtomID
        r_idx1: AtomID
        p_idx1, q_idx1, r_idx1 = args.atom

        xyzFile: GeometryXYZs = GeometryXYZs(inFile)
        xyzFile.method_read_xyz()
        from censo_ext.Tools.topo import Topo
        _topo: Topo = Topo(xyzFile, check=True)

        idx1_H: list[AtomID] = _topo.method_broken_bond_H(
            _bond_broken=(p_idx1, q_idx1), _print=False)
        idx0_H: list[IntpID] = [IntpID(x-1) for x in idx1_H]

    # Process xyz file
    for St in xyzFile.Sts:

        dxyz: npt.NDArray[np.float64] = St.coord[p_idx1-1].copy()
        St.coord -= dxyz  # type: ignore
        z_axis = (0, 0, np.sqrt(np.sum(np.square(St.coord[q_idx1-1]))))

        rotation_axis = St.coord[q_idx1-1] + z_axis

        Normalized_RotationAxis: npt.NDArray[np.float64] = np.array([
            0, 1, 0]) if np.linalg.norm(rotation_axis) == 0 else rotation_axis / np.linalg.norm(rotation_axis)

        R_pq = R.from_rotvec(np.pi*Normalized_RotationAxis)
        St.coord = R_pq.apply(St.coord)  # type: ignore

        Angle_qr = np.angle(complex(St.coord[r_idx1-1][0], complex(
            St.coord[r_idx1-1][1])))
        R_qr = R.from_euler('z', -Angle_qr)
        St.coord = R_qr.apply(St.coord)  # type: ignore

        for idx0, x in enumerate(St.coord):
            if idx0 in idx0_H:
                x[1] = -x[1]

    fileName: Path = inFile if args.replace else outFile
    print(f"    Saved to {fileName}")
    xyzFile.set_filename(fileName)
    xyzFile.method_save_xyz([])
    from censo_ext.Tools.topo import Topo
    _topo = Topo(xyzFile, check=True)


if __name__ == "__main__":
    main()

    # test
    # python3 xyzReturnOandZ.py -i ../tests/crest_conformers.xyz -a 30 45 47
    # python3 xyzReturnOandZ.py -i ../tests/crest_conformers.xyz --auto
