#!/usr/bin/env python
from pathlib import Path
import argparse
import numpy as np
import numpy.typing as npt
from censo_ext.Tools.utility import print_arguments
from censo_ext.Tools.xyzfile import GeometryXYZs
descr = """
________________________________________________________________________________
|                                          [08.18.2024] vitamin.cheng@gmail.com
| Transform from ANMR .out file to .dat file (Modification)
| Usage    : TsCNMR.py [options]
| [options]
| Input    : -i input out file [default anmrc.out] 
| Output   : -o output dat file [default anmrc.dat]    
| Extra    : --extra xyz file for carbon's height 
|            (depend on different numbers of C-H and will show in nmrplot.py)  
| Verbose  : -v verbose mode and show the detail 
| Print    : -p print the final result on screen
|______________________________________________________________________________
"""


def cml():
    """ Get args object from commandline interface.
        Needs argparse module."""
    parser = argparse.ArgumentParser(
        description=f"{descr}",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage=argparse.SUPPRESS,
    )

    parser.add_argument(
        "--extra",
        dest="extra",
        action="store",
        required=False,
        help="Provide the name of input xyz file for carbon's height",
    )

    parser.add_argument(
        "-o",
        "--output",
        dest="out",
        action="store",
        default="anmrc.dat",
        required=False,
        help="Provide the name of output file [default anmrc.dat]",
    )

    parser.add_argument(
        "-i",
        "--input",
        dest="file",
        action="store",
        default="anmrc.out",
        required=False,
        help="Provide the name of input file [default anmrc.out]",
    )

    parser.add_argument(
        "-p",
        "--print",
        dest="print",
        action="store_true",
        help="Print the final result of calculation on screen [default False]",
    )

    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        help="Verbose mode and show the detail of calculation [default False]",
    )

    args = parser.parse_args()


########## GLOBAL DECLARATIONS ##########


Scan_width = 32 * 1024
Size_Frq = 238.8955
Width_pp = Size_Frq / Scan_width
TMS_ppm = 14.5154
TMS_points = TMS_ppm / Width_pp
Line_width = 1
FWHM = Line_width / 1000
HMHM = FWHM / 2

########## END GLOBAL DECLARATIONS ##########


def Set_Peak_full_Spectra(x, multi, npData) -> None:
    for idx in range(Scan_width):
        npData[idx] = npData[idx] + \
            Lorentzian_Distribution(x-(idx*Width_pp-TMS_ppm), HMHM)*((multi+1))


def Lorentzian_Distribution(dx, mu):
    return 1 / 3.1415926 * (mu / (dx*dx + mu+mu))


def main(args: argparse.Namespace = argparse.Namespace()) -> None:
    if args == argparse.Namespace():
        args = cml()
    print_arguments()

    if not Path(args.file).exists:
        print(f"    {args.file} , the file is not exist ...")
        print("  Exit and Close the program !!!")
        exit(0)

    if not args.print:
        print("")
        print(f" Reading the {args.file} file ")

    import censo_ext.Tools.anmrfile as anmrfile
    inAnmr = anmrfile.Anmr()
    inAnmr.method_read_anmrSJ(args.file)

    if args.extra:
        import censo_ext.Tools.ml4nmr as ml4nmr
        xyzFile: GeometryXYZs = GeometryXYZs(args.extra)
        xyzFile.method_read_xyz()
        mol, neighbors, bond_order = ml4nmr.read_mol_neighbors_bond_order(
            xyzFile=xyzFile)
        List_nProton: list[int] = [x for x in bond_order.values()]

    List_ppm: list[float] = [ppm[3] for ppm in inAnmr.anmrS]
    npData: npt.NDArray = np.zeros((Scan_width, 2))

    for idx in range(len(List_ppm)):
        if args.extra:
            Set_Peak_full_Spectra(
                List_ppm[idx], List_nProton[idx], npData.T[1])  # type: ignore
        else:
            Set_Peak_full_Spectra(List_ppm[idx], 1, npData.T[1])

    # npData.T[1] = npData.T[1]*0.1

    npData.T[0] = np.arange(Scan_width)*Width_pp-TMS_ppm

    threshold: float = 0.001
    outData: npt.NDArray = npData[np.logical_not(npData[:, 1] < threshold)]
    outData: npt.NDArray = np.insert(
        outData, 0, (npData[0][0], threshold), axis=0)
    outData: npt.NDArray = np.insert(outData, len(outData),
                                     (npData[-1][0], threshold), axis=0)

    if args.out:
        from censo_ext.Tools.utility import save_simulation_spectra_file
        save_simulation_spectra_file(args.out, outData)

    if args.print:
        print(outData)


if __name__ == "__main__":
    main()
