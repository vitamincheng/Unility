#!/usr/bin/env python
from scipy.spatial.transform import Rotation as R
import argparse
from pathlib import Path
import numpy as np
import numpy.typing as npt
from censo_ext.Tools.utility import IntpID, print_arguments
from censo_ext.Tools.xyzfile import GeometryXYZs
descr = """
________________________________________________________________________________
| For search the confomrers from various angles of cleavage specifying two atoms                        
| Usage    : xyzSplit.py [options]                  
| Input    : -i Read xyz file [default traj.xyz]
| Output   : -o Save xyz file [default output.xyz] 
| [Options]
| Atom     : -a or --atom [1 2] idx of atom's number  
|              1 : Fixed atom
|              2 : Rotation axis atom (360 degrees) 
| nCut     : -c or cut Number of cut to make 360 degrees around the roation axis [default 3]
| Print    : -p Print output to screen [default False]
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
        help="Provide one output xyz file [defalut output.xyz]",
    )

    parser.add_argument(
        "-a",
        "--atom",
        dest="atoms",
        action="store",
        type=int,
        nargs=2,
        default=None,
        required=True,
        metavar=('FIXED', 'ROTATION'),
        help="two atom indics: first is fixed, second is rotation axis"
    )

    parser.add_argument(
        "-c",
        "--cut",
        dest="cuts",
        action="store",
        type=int,
        default=3,
        required=False,
        help="Number of cuts to make in 360 degrees around the rotation axis [default 3]",
    )

    parser.add_argument(
        "-p",
        "--print",
        dest="print",
        action="store_true",
        help="Print output to screen [default False]",
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

    if not args.cuts or not args.atoms:
        print("  Please input your atoms that you want to split ")
        print("  Exit and Close the program !!!")
        exit(0)

    from censo_ext.Tools.utility import delete_all_files
    if not args.print:
        delete_all_files(outFile)

    idx1_p, idx1_q = args.atoms
    nCutters: int = args.cuts

    xyzFile: GeometryXYZs = GeometryXYZs(inFile)
    xyzFile.method_read_xyz()
    from censo_ext.Tools.topo import Topo
    _topo = Topo(xyzFile)
    broken_bond_H: list[IntpID] = [
        IntpID(x-1) for x in _topo.method_broken_bond_H(_bond_broken=(idx1_q, idx1_p), _print=False)]

    for idx1_St, St in enumerate(xyzFile.Sts, 1):

        dxyz: npt.NDArray[np.float64] = St.coord[idx1_p-1].copy()
        inital: list[npt.NDArray[np.float64]] = St.coord.copy()

        for nCutter in range(nCutters):

            St.coord = inital.copy()
            St.coord -= dxyz  # type: ignore

            rotation_axis: npt.NDArray[np.float64] = St.coord[idx1_q-1]
            rotation_vector: npt.NDArray[np.float64] = rotation_axis / \
                np.linalg.norm(rotation_axis)

            r_pq = R.from_rotvec(2*np.pi*(nCutter/nCutters)*rotation_vector)

            for idx0 in broken_bond_H:
                St.coord[idx0] = r_pq.apply(St.coord[idx0])

            St.coord += dxyz

            if args.print:
                xyzFile.method_print([idx1_St])
            else:
                xyzFile.set_filename(outFile)
                xyzFile.method_save_xyz_append([idx1_St])
                _topo = Topo(xyzFile)

    if not args.print:
        print(f"    Save to the file : {outFile}")


if __name__ == "__main__":
    main()

    # For example
    #   xyzSplit.py -i tests/data/crest_conformers1.xyz -a 52 55 -c 3
    #   xyzSplit.py -i tests/data/crest_conformers1.xyz -a 52 55 -c 3 -p
