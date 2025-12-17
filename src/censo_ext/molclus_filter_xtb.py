#!/usr/bin/env python
from censo_ext.Tools.Parameter import Eh
from censo_ext.Tools.spectra import Boltzmann_Weighting
from censo_ext.Tools.xyzfile import GeometryXYZs
import argparse
import subprocess
from pathlib import Path
from censo_ext.Tools.utility import print_arguments
descr = """
________________________________________________________________________________
| For Filter of single point of GFN-xTB of molecules of xyz file
| Usages   : molclus_filter_xtb.py <geometry> [options]
| Input    : -i input file [default traj.xyz]
| Output   : -o output file [default isomers.xyz]
| [options]
| Method   : --method To set the method gfn0/gfn1/gfn2/gfnff [default gfn2]
| Solvent  : --alpb To set the solvent (for alpb mode and prefered choice)
|             CHCl3/DMSO/H2O(water)
| Solvent  : --gbsa To set the solvent (for gbsa mode)
|             methanol/CHCl3/DMSO/H2O(water)
| Charge   : --chrg to set the charge on the molecule [default 0]
| UHF      : --uhf to set the number of unpaired electrons [default 1]
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
        default="gfn2",
        help="Method : gfn0/gfn1/gfn2/gfnff [default gfn2]",
    )

    parser.add_argument(
        "--alpb",
        dest="alpb",
        action="store",
        required=False,
        help="Provide used the solvnet CHCl3/DMSO/H2O(water) ",
    )

    parser.add_argument(
        "--gbsa",
        dest="gbsa",
        action="store",
        required=False,
        help="Provide used the solvnet methanol/CHCl3/DMSO/H2O(water) ",
    )

    parser.add_argument(
        "--chrg",
        dest="chrg",
        action="store",
        type=int,
        required=False,
        default=0,
        help="to set the charge on the molecule [default 0]",
    )

    parser.add_argument(
        "--uhf",
        dest="uhf",
        action="store",
        type=int,
        required=False,
        default=1,
        help="to set the number of unpaired electrons [default 1]",
    )
    parser.add_argument(
        "--thr",
        dest="thr",
        action="store",
        type=float,
        required=False,
        default=100,
        help="to set the threshold of the electron energy of molecules (Kcal/mol) [default 100]",
    )
    parser.add_argument(
        "--check",
        dest="check",
        action="store_true",
        help="to Check the molbar is the same with first xyz structure [default False]",
    )

    args: argparse.Namespace = parser.parse_args()
    return args


def main(args: argparse.Namespace = argparse.Namespace()) -> None:

    if args == argparse.Namespace():
        args = cml()

    print_arguments()

    single_traj_Name = Path(".solo.xyz")
    temp_isomer_Name = Path(".isomers.xyz")
    xtb_cmd: str = ""

    # default to read the file
    inFile = Path(args.file)
    outFile = Path(args.out)
    xyzFile: GeometryXYZs = GeometryXYZs(inFile)
    xyzFile.method_read_xyz()

    # Check the molbar
    if args.check:
        print("  ===== Check the molbar =====")
        from molbar.barcode import get_molbars_from_coordinates
        inCoords: list = [x.coord for x in xyzFile.Sts]
        inNames: list = [x.names for x in xyzFile.Sts]
        molbars = get_molbars_from_coordinates(inCoords, inNames)

        std = molbars[0]
        print(len(molbars))
        for idx0, molbar in enumerate(molbars):
            if std != molbar:
                print(f"Index {idx0+1} in your xyz file have different molbar")
                for idx1, bar in enumerate(molbars):
                    print("")
                    print(f"Index of {idx1+1} : ")
                    print(f"{bar}")
                exit(1)

        print("  In your xyz file have the same molbar")
        print(f"{molbars[0]}")
        print("")

    # Default to xtb command: singel point energy
    print("  ===== single point energy of xtb =====")
    from censo_ext.Tools.utility import prog_IsExist
    prog = "xtb"
    prog_IsExist(prog)
    xtb_cmd += prog

    print(f" Inputted geometry file: {inFile}")
    xtb_cmd += f" {single_traj_Name}"
    print(" Loading basic information from the inputted geometry file ...")
    print(f" There are totally       {len(xyzFile)} geometries in the inputted geometry file\n")  # nopep8
    print(f" Setting method :  {args.method}")
    cmd_solvent = "vacuum"
    if args.alpb and (not args.gbsa):
        cmd_solvent = args.alpb
    elif args.gbsa:
        cmd_solvent = args.gbsa
    print(f" Setting solvent :  {cmd_solvent}")
    print(" Loading setting data ...")
    xtb_cmd += f" --{args.method}"
    if args.alpb and (not args.gbsa):
        xtb_cmd += f" --alpb {args.alpb}"
    elif args.gbsa:
        xtb_cmd += f" --gbsa {args.gbsa}"

    xtb_cmd += f" --chrg {args.chrg}  --uhf {args.uhf}"

    print(" All conformer in the inputted geometry file will be processed")
    print(" Cleaning old input and temporary files ...")
    print(" Running: rm isomers.xyz *.tmp")

    Energy: list[float] = []
    for idx1 in range(1, len(xyzFile)+1, 1):
        xyzFile.set_filename(single_traj_Name)
        xyzFile.method_save_xyz([idx1])
        print(f"                          *** Configuration         {idx1}  ****")  # nopep8
        print(f" Loading geometry	 {idx1}  from the inputted geometry file")      # nopep8
        print(" Generating  file...")
        subprocess.call(f"{xtb_cmd} > xtb.out", shell=True)
        idx1_str = (f"{idx1:05d}")
        print(f" Running:  {xtb_cmd} > {idx1_str}.out")

        intp_energy_lines: int | None = None
        # print("singe point")
        lines: list[str] = open("xtb.out", "r").readlines()
        for idy0, y in enumerate(lines):
            if r"TOTAL ENERGY" in y:
                intp_energy_lines = idy0
        if intp_energy_lines:
            xyzFile.Sts[idx1 -
                        1].comment_energy = float(lines[intp_energy_lines].split()[3])
            Energy.append(float(lines[intp_energy_lines].split()[3])*Eh)

    # save singe point energy of xtb
    xyzFile.method_rewrite_comment()
    xyzFile.method_comment_new()
    xyzFile.set_filename(outFile)

    # print the Boltzmann weighting
    print("")
    print("  ===== Boltzmann Distribution =====")
    print(f"  threshold energy = {args.thr} (kcal/mol)")
    print("")
    print("  Boltzmann Weighting Table")

    import numpy as np
    import numpy.typing as npt
    np_Energy = np.array(Energy)*Eh
    np_Energy = (np_Energy-np_Energy.min())
    intp_Energy: npt.NDArray[np.intp] = np.argsort(np_Energy)
    for idx0, x in enumerate(intp_Energy.copy()):
        if np_Energy[x] >= args.thr:
            intp_Energy = np.delete(intp_Energy, idx0)

    BW: npt.NDArray[np.float64] = Boltzmann_Weighting(
        np_Energy[intp_Energy], TEMP=298.15)

    zip_energy: zip[tuple[npt.NDArray[np.intp], npt.NDArray[np.float64], npt.NDArray[np.float64]]] = zip(
        intp_Energy+1, np_Energy[intp_Energy], BW)

    print("  index1           Energy (kcal/mol)             BW")
    for x, y, z in zip_energy:
        print(f"{x:8d}           {y:17.10f}       {z:8.4f}")
    xyzFile.method_save_xyz((intp_Energy+1).tolist())

    subprocess.call(
        "rm -rf charges wbo xtb.out xtbrestart xtbtopo.mol xtbopt* .xtboptok", shell=True)
    from censo_ext.Tools.utility import delete_all_files
    delete_all_files(temp_isomer_Name, single_traj_Name)


if __name__ == "__main__":
    main()

#   test
#   python3 molclus_xtb.py -i ../tests/crest_conformers.xyz --alpb CHCl3
