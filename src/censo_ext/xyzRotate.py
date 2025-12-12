#!/usr/bin/env python
import argparse
from pathlib import Path
from censo_ext.Tools.utility import print_arguments
from censo_ext.Tools.xyzfile import GeometryXYZs
descr = """
________________________________________________________________________________
| For rotation the confomrer from specific angle of cleavage specifying two atoms                        
| Usage    : xyzRotate.py [options]                  
| Input    : -i Read xyz file [default traj.xyz]
| Output   : -o Save xyz file [default output.xyz] 
| [Options]
| Atom     : -a or --atom [1 2] idx of atom's number  
|              1 : Fixed atom
|              2 : Rotation axis atom (360 degrees) 
| nCut     : -c or cut Number of cut to make 360 degrees around the roation axis 
| sPecific : -s Specific number of the cut's number from 1 to (nCuts-1) 
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
        "-s",
        "--specific",
        dest="spec",
        action="store",
        type=int,
        default=1,
        required=False,
        help="Specific number of cut number the rotation axis [default 1]",
    )

    parser.add_argument(
        "-p",
        "--print",
        dest="print",
        action="store_true",
        help="Print output to screen [default False]",
    )
    parser.add_argument(
        "--check",
        dest="check",
        action="store_true",
        help="Check hydrogen atom is only one bond [default False]",
    )

    args: argparse.Namespace = parser.parse_args()
    return args


def main(args: argparse.Namespace = argparse.Namespace()) -> None:

    if args == argparse.Namespace():
        args = cml()
    print_arguments()

    from censo_ext.Tools.utility import IsExist
    inFile: Path = Path(args.file)
    outFile = Path(args.out)
    IsExist(inFile)
    xyzFile: GeometryXYZs = GeometryXYZs(inFile)

    if not args.cuts or not args.atoms:
        print("  Please input your atoms that you want to split ")
        print("  Exit and Close the program !!!")
        exit(0)

    from censo_ext.Tools.utility import delete_all_files
    if not args.print:
        delete_all_files(outFile)

    idx1_p, idx1_q = args.atoms

    if args.spec >= args.cuts or args.spec <= 0:
        print("  Specific number of cut' number is error !!!")
        print("  Exit and Close the program !!!")
        exit(1)

    xyzFile.method_read_xyz()
    xyzFile.Method_xyzRotate(_check=args.check, idx1_p=idx1_p,
                             idx1_q=idx1_q, _cuts=args.cuts, _nspec=args.spec)
    xyzFile.set_filename(outFile)
    xyzFile.method_save_xyz([])
    from censo_ext.Tools.topo import Topo
    _topo = Topo(xyzFile, check=args.check)

    if not args.print:
        print(f"    Save to the file : {outFile}")


if __name__ == "__main__":
    main()

    # For example
    #   xyzSplit.py -i tests/data/crest_conformers1.xyz -a 52 55 -c 3
    #   xyzSplit.py -i tests/data/crest_conformers1.xyz -a 52 55 -c 3 -p
