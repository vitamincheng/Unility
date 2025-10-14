#!/usr/bin/env python
import argparse
from scipy.spatial.transform import Rotation as R
import numpy as np
import numpy.typing as npt
from censo_ext.Tools.utility import print_arguments
from censo_ext.Tools.xyzfile import GeometryXYZs
from pathlib import Path

descr = """
________________________________________________________________________________
| For Return to origin and mirror on the xz plane
| Usage    : xyzMirror.py <geometry> [options]                  
| Input    : -i xyz file [default traj.xyz]
| Output   : -o Saved xyz file [default output.xyz] 
| [Options]
| Atom     : -a or --atom [1 2 3] idx of atom  
|              1 : Fixed atom and return origin (mirror)
|              2 : Rotation atom and z axis (not mirror)
|              3 : Rotation atom and lay on xz plane (not mirror)
| Replace  : -r Replace the input file
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
        required=False,
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
        idx0_H = [*range(Nums)]
        p_idx = 1
        q_idx = 1
        r_idx = 1
    else:
        p_idx: int
        q_idx: int
        r_idx: int
        p_idx, q_idx, r_idx = args.atom

        xyzFile: GeometryXYZs = GeometryXYZs(inFile)
        xyzFile.method_read_xyz()

        from censo_ext.Tools.topo import Topo
        args_x: dict = {"file": inFile, "bond_broken": (p_idx, q_idx),
                        "print": False, "debug": False}
        Sts_topo: Topo = Topo(args_x["file"])
        idx1_H: list[int] = Sts_topo.method_broken_bond_H(
            args=argparse.Namespace(**args_x))
        idx0_H: list[int] = [x-1 for x in idx1_H]

    # Process xyz file
    for idx_St in range(len(xyzFile)):

        dxyz: npt.NDArray[np.float64] = xyzFile.Sts[idx_St].coord[p_idx-1].copy()
        xyzFile.Sts[idx_St].coord -= dxyz  # type: ignore
        z_axis = (0, 0, np.sqrt(
            np.sum(np.square(xyzFile.Sts[idx_St].coord[q_idx-1]))))

        rotation_axis = xyzFile.Sts[idx_St].coord[q_idx-1] + z_axis

        Normalized_RotationAxis: npt.NDArray[np.float64] = np.array([
            0, 1, 0]) if np.linalg.norm(rotation_axis) == 0 else rotation_axis / np.linalg.norm(rotation_axis)

        R_pq = R.from_rotvec(np.pi*Normalized_RotationAxis)
        xyzFile.Sts[idx_St].coord = R_pq.apply(
            xyzFile.Sts[idx_St].coord)  # type: ignore

        Angle_qr = np.angle(complex(xyzFile.Sts[idx_St].coord[r_idx-1][0], complex(
            xyzFile.Sts[idx_St].coord[r_idx-1][1])))
        R_qr = R.from_euler('z', -Angle_qr)
        xyzFile.Sts[idx_St].coord = R_qr.apply(
            xyzFile.Sts[idx_St].coord)  # type: ignore

        # from icecream import ic
        # ic(idx0_H)
        # ic(xyzFile.Sts[idx_St].coord)

        for idx0, x in enumerate(xyzFile.Sts[idx_St].coord):
            if idx0 in idx0_H:
                x[1] = -x[1]

        # from icecream import ic
        # ic(xyzFile.Sts[idx_St].coord)

    fileName: Path = inFile if args.replace else outFile
    print(f"    Saved to {fileName}")
    xyzFile.set_filename(fileName)
    xyzFile.method_save_xyz([])


if __name__ == "__main__":
    main()

    # test
    # python3 xyzReturnOandZ.py -i ../tests/crest_conformers.xyz -a 30 45 47
    # python3 xyzReturnOandZ.py -i ../tests/crest_conformers.xyz --auto
