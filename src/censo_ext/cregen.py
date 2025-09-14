#!/usr/bin/env python
import argparse
from sys import argv as sysargv
import subprocess
from pathlib import Path

descr = """
________________________________________________________________________________
| Modification of cregen using rthr(0.175) bthr(0.03) ethr (0.15) ewin(4.0)
| Usage    : cregen.py <geometry> [options]
|
| Input    : -i one xyz file [default isomers.xyz]
| Output   : -o one xyz file [default cluster.xyz]
| [Options]
| rthr     : --rthr set RMSD threshold Angstrom [default 0.175]
| bthr     : --bthr set lower bound for the rotatoional constant threshold 
|                   [default 0.03]
| ethr     : --ethr set energy threshold between conformer pair in kcal/mol 
|                   [default 0.15]
| ewin     : --ewin set the energy threshold in kcal/mol [default 4]
|______________________________________________________________________________
"""


def cml() -> argparse.Namespace:
    """ Get args object from commandline interface.
        Needs argparse module."""
    parser = argparse.ArgumentParser(
        description=descr,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage=argparse.SUPPRESS,
    )

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
        args = cml()
    # if not args.print :
    print(descr)  # Program description
    print(f"    provided arguments: {" ".join(sysargv)}")

    inFile: Path = Path(args.file)
    outFile: Path = Path(args.out)
    isomers: Path = Path("isomers.xyz")
    cluster: Path = Path("cluster.xyz")

    from censo_ext.Tools.utility import IsExist
    IsExist(inFile)

    if args.file != isomers:
        subprocess.call(f"cp {inFile} {isomers}", shell=True)
        print(f"  cp {inFile} {isomers}")

    prog = "crest"
    from censo_ext.Tools.utility import prog_IsExist
    prog_IsExist(prog)

    crest_cmd: str = f"{prog} {isomers} --cregen {isomers} --rthr {args.rthr} -- bthr {args.bthr} --ethr {args.ethr} --ewin {args.ewin} > isomers.out"  # nopep8

    subprocess.call(crest_cmd, shell=True)
    print(f"  {crest_cmd}")

    subprocess.call(f"mv -f crest_ensemble.xyz {cluster}", shell=True)
    print(f"  mv -f crest_ensemble.xyz {cluster}")
    subprocess.call(
        f"xyzSerial.py -i {cluster} --new --print > tmp && mv -f tmp {cluster}", shell=True)
    if outFile != cluster:
        subprocess.call(f"mv -f {cluster} {outFile}", shell=True)
        print(f"  mv -f {cluster} {outFile}")
    if inFile != isomers:
        subprocess.call(f"rm -rf {isomers}", shell=True)

    subprocess.call(
        "rm -rf coord* cre_members crest.energies crest_best.xyz scoord.1 struc.xyz isomers.xyz.sorted crest_input_copy.xyz crest.restart", shell=True)


if __name__ == "__main__":
    main()
