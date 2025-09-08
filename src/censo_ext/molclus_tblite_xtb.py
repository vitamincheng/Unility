#!/usr/bin/env python.12
from censo_ext.Tools.xyzfile import GeometryXYZs
import argparse
from pathlib import Path
import tblite.interface as tb
from berny import Berny, geomlib, angstrom
import numpy as np
descr = """
________________________________________________________________________________
| For GFN2-xTB by using tblite libary (python)  
| Usages   : molclus_tblite_xtb.py <geometry> [options]
| [options]
| Input    : -i input file [default traj.xyz]
| Output   : -o output file [default isomers.xyz]
| Opt      : --opt To optimize energy [default False]
| Method   : --method To set the method GFN1-xTB/GFN2-xTB [default GFN2-xTB]
| Solvent  : --alpb To set the solvent (for alpb mode and prefered choice)
|             CHCl3/DMSO/H2O(water)  
| Solvent  : --gbsa To set the solvent (for gbsa mode)
|             methanol/CHCl3/DMSO/H2O(water)  
| Charge   : --chrg to set the charge on the molecule [default 0]
| UHF      : --uhf to set the number of unpaired electrons [default 1]
|______________________________________________________________________________
"""


def cml(descr) -> argparse.Namespace:
    """ Get args object from commandline interface.
        Needs argparse module."""
    parser = argparse.ArgumentParser(
        description="",
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
        help="Provide input_file name [default traj.xyz]",
    )

    parser.add_argument(
        "-o",
        "--output",
        dest="out",
        action="store",
        required=False,
        default="isomers.xyz",
        help="Provide output_file name [default isomers.xyz]",
    )

    parser.add_argument(
        "--method",
        dest="method",
        action="store",
        required=False,
        default="GFN2-xTB",
        help="Method : GFN1-xTB/GFN2-xTB [default GFN2-xTB]",
    )

    parser.add_argument(
        "--alpb",
        dest="alpb",
        action="store",
        required=False,
        help="Provide used the solvnet chloroform/water ",
    )

    parser.add_argument(
        "--gbsa",
        dest="gbsa",
        action="store",
        required=False,
        help="Provide used the solvnet methanol/CHCl3/DMSO/H2O(water) [pending]",
    )

    parser.add_argument(
        "--chrg",
        dest="chrg",
        action="store",
        type=int,
        required=False,
        default=0,
        help="to set the charge on the molecule [default 0] [pending]",
    )

    parser.add_argument(
        "--uhf",
        dest="uhf",
        action="store",
        type=int,
        required=False,
        default=1,
        help="to set the number of unpaired electrons [default 1] [pending]",
    )

    parser.add_argument(
        "--opt",
        dest="opt",
        action="store_true",
        help="Optimize energy [default False]",
    )

    args: argparse.Namespace = parser.parse_args()
    return args


def main(args: argparse.Namespace = argparse.Namespace()) -> None:

    if args == argparse.Namespace():
        args = cml(descr)

    single_traj_Name = ".single_traj.xyz"
    infile: GeometryXYZs = GeometryXYZs(args.file)
    infile.method_read_xyz()

    print(f" Inputted geometry file: {args.file}")
    print(" Loading basic information from the inputted geometry file ...")
    print(f" There are totally       {len(infile)} geometries in the inputted geometry file\n")  # nopep8
    print(f" Setting method :  {args.method}")
    cmd_solvent = "vacuum"
    if args.alpb:
        cmd_solvent = args.alpb
    elif args.gbsa:
        cmd_solvent = args.gbsa
    print(f" Setting solvent :  {cmd_solvent}")
    print(" Loading setting data ...")

    print(" All conformer in the inputted geometry file will be processed")
    print(" Cleaning old input and temporary files ...")
    print(" Running: rm isomers.xyz *.tmp")

    for idx in range(1, len(infile)+1, 1):
        # idx_str : str = f"{[idx]:05d}"
        infile.set_filename(Path(single_traj_Name))
        infile.method_save_xyz([idx])
        print(f"                          *** Configuration         {idx}  ****")  # nopep8
        print(f" Loading geometry	 {idx}  from the inputted geometry file")      # nopep8
        print(" Generating  file...")

        optimizer = Berny(geomlib.readfile(single_traj_Name))
        geom = next(optimizer)
        elements = [symbol for symbol, _ in geom]
        initial_coordinates = np.asarray(
            [coordinate for _, coordinate in geom])
        xtb = tb.Calculator(args.method, tb.symbols_to_numbers(elements),
                            initial_coordinates * angstrom)
        xtb.set("verbosity", 0)
        if args.alpb:
            xtb.add("alpb-solvation", cmd_solvent)
        elif args.gbsa:
            xtb.add("gbsa-solvation", cmd_solvent)
        results = xtb.singlepoint()

        print(" Running:  xtb_tblite_cmd")

        if args.opt:
            initial_energy = results["energy"]
            initial_gradient = results["gradient"]
            optimizer.send((initial_energy, initial_gradient / angstrom))

            trajectory = [
                (initial_energy, initial_gradient, initial_coordinates)]
            for geom in optimizer:
                coordinates = np.asarray(
                    [coordinate for _, coordinate in geom])
                xtb.update(positions=coordinates * angstrom)
                results = xtb.singlepoint(results)

                energy = results["energy"]
                gradient = results["gradient"]
                optimizer.send((energy, gradient / angstrom))

                trajectory.append((energy, gradient, coordinates))
            infile.Sts[idx-1].comment_energy = trajectory[-1][0]
            infile.Sts[idx-1].coord = list(trajectory[-1][2])

        else:
            infile.Sts[idx-1].comment_energy = results["energy"]

    infile.method_rewrite_comment()
    infile.method_comment_new()
    infile.set_filename(args.out)
    infile.method_save_xyz([])

    from censo_ext.Tools.utility import delete_all_files
    delete_all_files(single_traj_Name)


if __name__ == "__main__":
    main()

#   test
#   python molclus_tblite_xtb.py  -i ../../tests/data/06.EthylAcetate/02.ORCA_r2SCAN_3C/traj.xyz --alpb chloroform --opt
#   python molclus_tblite_xtb.py  -i ../../tests/data/06.EthylAcetate/02.ORCA_r2SCAN_3C/traj.xyz --alpb chloroform
#   python molclus_tblite_xtb.py  -i ../../tests/data/06.EthylAcetate/02.ORCA_r2SCAN_3C/traj.xyz
