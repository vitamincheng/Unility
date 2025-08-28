#!/usr/bin/env python
import argparse
from sys import argv as sysargv
import subprocess

descr = """
________________________________________________________________________________
|                                          [07.18.2023] vitamin.cheng@gmail.com
| Modification of cregen using rthr(0.175) bthr(0.03) ethr (0.15) ewin(4.0)
| Usage : cregen.py <geometry> [options]
| [Options]
| Input    : -i one xyz file [default isomers.xyz]
| Output   : -o one xyz file [default cluster.xyz]
| rthr     : --rthr set RMSD threshold Angstrom [default 0.175]
| bthr     : --bthr set lower bound for the rotatoional constant threshold 
|                   [default 0.03]
| ethr     : --ethr set energy threshold between conformer pair in kcal/mol 
|                   [default 0.15]
| ewin     : --ewin set the energy threshold in kcal/mol [default 4]
| Package  : Tools
|______________________________________________________________________________
"""


def cml(descr) -> argparse.Namespace:
    """ Get args object from commandline interface.
        Needs argparse module."""
    parser = argparse.ArgumentParser(
        description="",
        #        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage=argparse.SUPPRESS,
    )  # argparse.RawDescriptionHelpFormatter) #,

    parser.add_argument(
        "-i",
        "--input",
        dest="file",
        action="store",
        required=False,
        default="isomers.xyz",
        help="Provide one xyz file to reorganize the serial numbers [default isomers.xyz]",
    )

    parser.add_argument(
        "-o",
        "--output",
        dest="out",
        action="store",
        required=False,
        default="cluster.xyz",
        help="Provide one xyz file to save the data [default cluster.xyz]",
    )

    parser.add_argument(
        "--rthr",
        dest="rthr",
        action="store",
        required=False,
        default=0.175,
        type=float,
        help="set RMSD threshold Angstrom [default 0.175]",
    )

    parser.add_argument(
        "--bthr",
        dest="bthr",
        action="store",
        required=False,
        default=0.03,
        type=float,
        help="set lower bound for the rotatoional constant threshold [default 0.03]",
    )

    parser.add_argument(
        "--ethr",
        dest="ethr",
        action="store",
        required=False,
        default=0.15,
        type=float,
        help="set energy threshold between conformer pair in kcal/mol [default 0.15]",
    )

    parser.add_argument(
        "--ewin",
        dest="ewin",
        action="store",
        required=False,
        default=4,
        type=float,
        help="ewin set the energy threshold in kcal/mol [default 4]",
    )

    args: argparse.Namespace = parser.parse_args()
    return args


def main(args: argparse.Namespace = argparse.Namespace()) -> None:

    if args == argparse.Namespace():
        args = cml(descr)
    # if not args.print :
    print(descr)  # Program description
    print(f"    provided arguments: {" ".join(sysargv)}")

    from censo_ext.Tools.utility import IsExist
    IsExist(args.file)

    if args.file != "isomers.xyz":
        subprocess.call(f"cp {args.file} isomers.xyz", shell=True)
        print(f"  cp {args.file} isomers.xyz")

    crest_cmd: str = f"crest isomers.xyz --cregen isomers.xyz --rthr {str(args.rthr)} -- bthr {str(args.bthr)} --ethr {str(args.ethr)} --ewin {str(args.ewin)} > isomers.out"  # nopep8

    from censo_ext.Tools.utility import program_IsExist
    program_IsExist("crest")

    subprocess.call(crest_cmd, shell=True)
    print("  " + crest_cmd)
    subprocess.call("mv -f crest_ensemble.xyz cluster.xyz", shell=True)
    print("  mv -f crest_ensemble.xyz cluster.xyz")
    subprocess.call(
        "xyzSerial.py -i cluster.xyz --new --print > tmp && mv -f tmp cluster.xyz", shell=True)
    if args.out != "cluster.xyz":
        subprocess.call(f"mv -f cluster.xyz {args.out}", shell=True)
        print(f"  mv -f cluster.xyz {args.out}")
    if args.file != "isomers.xyz":
        subprocess.call("rm -rf isomers.xyz", shell=True)

    subprocess.call(
        "rm -rf coord* cre_members crest.energies crest_best.xyz scoord.1 struc.xyz isomers.xyz.sorted crest_input_copy.xyz crest.restart", shell=True)


if __name__ == "__main__":
    main()
