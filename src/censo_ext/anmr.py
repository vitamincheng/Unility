#!/usr/bin/env python
import copy
import numpy as np
import numpy.typing as npt
import argparse
from icecream import ic
from pathlib import Path

from censo_ext.Tools.anmrfile import Anmr

descr = """
________________________________________________________________________________
| anmr.py
| Usages   : anmr.py [options]
| [options]
| Output   : -o output file [default output.dat]
| Dir      : -D the directory of input files(CONF) [default .]
| mf       : -mf magnetic frequency of scan nmr [default 500.0]
| lw       : -lw line width of scan nmr [1.0 for H, 20 for C]
| auto     : -auto --auto automated to adjust the threshold of J and AB quartet [default False]
| decoups  : -de --decoups deCoupings mode [defalut False]
| average  : -av load the average folder data to plot spectra [default False]
| thr      : -t -thr threshold of coupling constant (J) [default 0.30]
| thrtab   : -tab -thrab threshold of AB quartet (JCoups / diff chemical shift) [default 0.020]
| tbpent   : -tb threshold of AB quartet bond pententration distance [default 4]
| mss      : -mss max of spin numbers [default 10]
|            if your computer have slowly CPU, try to use mss 4 or 5.
| Cutoff   : --cutoff cutoff limits in quantum chemistry in nmr [default 0.001]
| Start    : -start start ppm of plotting spectra [default from data]
| End      : -end end ppm of plotting spectra [default from data]
| ref      : reference standard - see .anmrrc file
| BOBYQA   : -b --bobyqa BOBYQA mode [default False]
| JSON     : -j --json Read the raw data of every single peak [if is -1(All)]
| ascal    : -ascal chemical shift scaling a if the reference is absent [pending]
| bscal    : -bcsal chemical shift scaling b if the reference is absent [pending]
|______________________________________________________________________________
"""


def cml() -> argparse.Namespace:
    """ Get args object from commandline interface.
        Needs argparse module."""
    parser = argparse.ArgumentParser(
        description="",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage=argparse.SUPPRESS,
    )

    parser.add_argument(
        "-o",
        "--output",
        dest="out",
        action="store",
        required=False,
        default="output.dat",
        help="Provide output_file name [default output.dat]",
    )

    parser.add_argument(
        "-D",
        "--dir",
        dest="dir",
        action="store",
        required=False,
        default=".",
        help="Provide the directory name [default .]",
    )

    parser.add_argument(
        "-mf",
        "--magnfreq",
        dest="mf",
        action="store",
        type=float,
        required=False,
        default=500.0,
        help="magnetic frequency of scan nmr [default 500.0]",
    )
    parser.add_argument(
        "-lw",
        "-linewidth",
        dest="lw",
        action="store",
        type=float,
        required=False,
        help="line width of scan nmr [default 1.0 for H, 20.0 for C]",
    )

    parser.add_argument(
        "-ascal",
        dest="ascal",
        action="store",
        type=float,
        required=False,
        help="chemical shift scaling a, "
        "For   "
        "For   ",
    )

    parser.add_argument(
        "-bscal",
        dest="bscal",
        action="store",
        type=float,
        required=False,
        help="chemical shift scaling b, "
        "For  "
        "For  ",
    )

    parser.add_argument(
        "-t",
        "--thr",
        dest="thr",
        action="store",
        type=float,
        required=False,
        help="threshold of J coupling constant Hz [default 0.30]",
    )

    parser.add_argument(
        "-tab",
        "--thrab",
        dest="thrab",
        action="store",
        type=float,
        required=False,
        default=0.020,
        help="threshold of AB quartet (J/chemical shift) [default 0.020]",
    )

    parser.add_argument(
        "-tb",
        "--tbpent",
        dest="tb",
        action="store",
        type=int,
        required=False,
        default=4,
        help="threshold of AB quartet bond pententration distance [default 4]",
    )

    parser.add_argument(
        "-cut",
        "--cutoff",
        dest="cutoff",
        action="store",
        type=float,
        default=0.001,
        help="cutoff limits in quantum chemistry in nmr [default 0.001]",
    )

    parser.add_argument(
        "-start",
        "--start",
        dest="start",
        action="store",
        type=float,
        required=False,
        help="start ppm of plotting spectra [default from data]",
    )

    parser.add_argument(
        "-end",
        "--end",
        dest="end",
        action="store",
        type=float,
        required=False,
        help="end ppm of plotting spectra [default from data]",
    )

    parser.add_argument(
        "-av",
        "--average",
        dest="average",
        action="store_true",
        help="Load the averager foler data to plot spectra when use auto argument [default False]",
    )

    parser.add_argument(
        "-b",
        "--bobyqa",
        dest="bobyqa",
        action="store_true",
        help="BOBYQA mode [default False]",
    )

    parser.add_argument(
        "-de",
        "--decoups",
        dest="decoups",
        action="store_true",
        help="Decouplings mode [important] for carbon spectra [default False]",
    )

    parser.add_argument(
        "-j",
        "--json",
        dest="json",
        action="store",
        type=int,
        nargs="+",
        required=False,
        help="Read the jason file of raw single peak [if is -1(All)]",
    )

    parser.add_argument(
        "-mss",
        "--mss",
        dest="mss",
        action="store",
        type=int,
        required=False,
        default=10,
        help="max of spin numbers [default 10]",
    )

    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        help="Verbose mode and show the detail of data [defalut False]",
    )

    parser.add_argument(
        "-auto",
        "--auto",
        dest="auto",
        action="store_true",
        help="Automated to adjust the threshold of AB quartet [defalut False]",
    )

    args: argparse.Namespace = parser.parse_args()
    return args


def setup_anmr(args: argparse.Namespace) -> Anmr:
    """Set up and return an Anmr object by reading required files.

    This function creates an Anmr object and reads the necessary configuration
    files including the .anmrrc file and nucleotide information file.

    Args:
        args (argparse.Namespace): Command line arguments containing:
            - dir (str): Directory path for the Anmr object
            - verbose (bool): Verbose output flag

    Returns:
        Anmr: Configured Anmr object with loaded files

    Example:
        >>> args = parser.parse_args()
        >>> anmr_obj = setup_anmr(args)
    """
    # Create an Anmr object and read required files
    # Import necessary modules for processing
    from censo_ext.Tools.anmrfile import Anmr
    inAnmr: Anmr = Anmr(Dir=args.dir, verbose=args.verbose)
    inAnmr.method_read_anmrrc()
    inAnmr.method_read_nucinfo()
    return inAnmr


def process_average_data(inAnmr: Anmr, args: argparse.Namespace) -> None:
    """Process average ORCA SJ data based on existing files or generate new ones.

    This function checks if average ORCA SJ data exists and loads it if requested.
    If not, it processes all ORCA files to generate the average data.

    Args:
        inAnmr (Anmr): Anmr object containing the data to process
        args (argparse.Namespace): Command line arguments containing:
            - average (bool): Flag to indicate if averaging should be done
            - bobyqa (bool): Flag to indicate if bobyqa method should be used

    Raises:
        FileNotFoundError: If there's an issue with loading existing average ORCA SJ data

    Example:
        >>> process_average_data(anmr_obj, args)
    """
    if inAnmr.get_avg_orcaSJ_Exist() and args.average:
        if not inAnmr.method_load_avg_orcaSJ(bobyqa_bool=args.bobyqa):
            raise FileNotFoundError(
                "Something wrong in your Average orcaSJ data !!!")
    else:
        # Process all ORCA files and generate average data
        inAnmr.method_read_enso()
        inAnmr.method_read_folder_orcaSJ()
        for idx1, Active in enumerate(inAnmr.get_Anmr_Active(), 1):
            if idx1 == 1:  # only one Active nuclear element
                inAnmr.method_filter_active_orcaSJ(Active)
            elif idx1 > 1:
                print("  Only for ONE Active Nuclear element, waiting to build")
                print("  Exit and Close the program !!!")
                exit(0)
        inAnmr.method_update_equiv_orcaSJ()
        inAnmr.method_avg_orcaSJ()
        inAnmr.method_save_avg_orcaSJ()


def preprocess_spin_system(inAnmr: Anmr, args: argparse.Namespace) \
    -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], list[int],
             int | None, int | None, argparse.Namespace]:
    """Preprocess spin system data for NMR analysis based on active nuclear element.

    This function extracts spin parameters and coupling constants from the Anmr object
    and processes them according to the active nuclear element (C or H). It handles
    different processing workflows for carbon and hydrogen nuclei, including:
    - Reading molecular structure for carbon
    - Identifying equivalent hydrogens for hydrogen
    - Setting appropriate ranges and DPI values
    - Validating data consistency

    Args:
        inAnmr (Anmr): Anmr object containing the NMR data to process
        args (argparse.Namespace): Command line arguments containing:
            - mf (float): Multiplication factor for spin parameters
            - lw (float): Line width parameter
            - thr (float): Threshold parameter

    Returns:
        tuple: A tuple containing:
            - inSParams (npt.NDArray[np.float64]): Spin parameters array
            - inJCoups (npt.NDArray[np.float64]): Coupling constants array
            - inHydrogen (list[int]): List of hydrogen atom indices
            - Active_range (int | None): Active range for the spectrum
            - dpi (int | None): Dots per inch for plotting
            - args (argparse.Namespace): Updated arguments

    Raises:
        ValueError: If there are inconsistencies in the ORCA data files
        SystemExit: If unsupported active nuclear elements are encountered

    Example:
        >>> s_params, j_coups, hydrogen, range_val, dpi_val, updated_args = preprocess_spin_system(anmr_obj, args)
    """

    inAnmr.avg_orcaSJ.method_print_orcaS()
    inAnmr.avg_orcaSJ.method_print_orcaJ()

    # Extract spin parameters and coupling constants from the Anmr object
    inSParams: npt.NDArray[np.float64] = np.array(
        list(inAnmr.avg_orcaSJ.SParams.values()))*args.mf
    inJCoups: npt.NDArray[np.float64] = np.array(inAnmr.avg_orcaSJ.JCoups)

    # Initialize variables for processing based on active nuclear element
    dpi: int | None = None
    Active_range: int | None = None
    inHydrogen: list[int] = []
    inFile: Path = Path("crest_conformers.xyz")

    # Process different nuclear element (C or H)
    for idx1, Active in enumerate(inAnmr.get_Anmr_Active(), 1):
        if idx1 == 1:  # only one Active nuclear element
            if Active == 'C':
                # Carbon processing - read molecular structure using ML4NMR tool
                inHydrogen, Active_range, dpi = _preprocess_carbon_spin_system(
                    inAnmr, args, inFile)

                # Set all coupling constants to zero for carbon processing
                inJCoups = np.zeros_like(inJCoups)

                # validate that the atom types match between ORCA files
                if set(inAnmr.avg_orcaSJ.idx1Atoms.values()) != set(Active):
                    raise ValueError(
                        "  Yours Average orcaS.out and orcaJ.out have something errors !!!")

            elif Active == 'H':
                # Hydrogen processing - identify equivalent hydrogens from magnetization data
                inHydrogen, Active_range, dpi = _preprocess_hydrogen_spin_system(
                    inAnmr, args, inFile)

                # Validate that the atom types match between ORCA files
                if set(inAnmr.avg_orcaSJ.idx1Atoms.values()) != set(Active):
                    raise ValueError(
                        "  Your orcaS.out have Something error !!!")
            else:
                print("  Other Active Nuclear element, waiting to build")
                print("  Exit and Close the program !!!")
                exit(0)

        elif idx1 > 1:
            print("  Only for ONE Active Nuclear element, waiting to build")
            print("  Exit and Close the program !!!")
            exit(0)

    return inSParams, inJCoups, inHydrogen, Active_range, dpi, args


def _preprocess_carbon_spin_system(inAnmr: Anmr, args: argparse.Namespace, inFile: Path) -> tuple[list[int], int, int]:
    """Process carbon spin system data.

    Args:
        inAnmr (Anmr): Anmr object containing the NMR data
        args (argparse.Namespace): Command line arguments
        inFile (Path): Input file path for molecular structure

    Returns:
        tuple[list[int], int, int]: Hydrogen counts, Active range, and DPI values
    """
    from censo_ext.Tools.ml4nmr import read_mol_neighbors_bond_order
    *_, bond_order = read_mol_neighbors_bond_order(
        inAnmr.get_Dir() / inFile)

    # For C/CH/CH2/CH3 from 0 1 2 3 to 1 2 3 4 for Carbon spectra
    # Convert bond order values to hydrogen counts (add 1 to each value)
    inHydrogen = [value+1 for value in bond_order.values()]
    Active_range = 200
    dpi = 500
    if not args.lw:
        args.lw = 20
    if not args.thr:
        args.thr = args.lw * 0.3

    return inHydrogen, Active_range, dpi


def _preprocess_hydrogen_spin_system(inAnmr: Anmr, args: argparse.Namespace, inFile: Path) -> tuple[list[int], int, int]:
    """Process hydrogen spin system data.

    Args:
        inAnmr (Anmr): Anmr object containing the NMR data
        args (argparse.Namespace): Command line arguments
        inFile (Path): Input file path for molecular structure

    Returns:
        tuple[list[int], int, int]: Hydrogen counts, Active range, and DPI values
    """
    idx1_nMagEqvHydrogens: dict[int, int] = {}
    for key in inAnmr.avg_orcaSJ.idx1Atoms.keys():
        idx1_nMagEqvHydrogens[key] = inAnmr.nMagnetEqvs[key]

    # Remove hydrogen atoms that are part of acid groups
    for y in inAnmr.get_idx1_acid_atoms_NoShow_RemoveH(
            inAnmr.get_Dir()/inFile):
        if y in idx1_nMagEqvHydrogens.keys():
            del idx1_nMagEqvHydrogens[y]

    # Extract hydrogen counts for each equivalent group
    inHydrogen = [
        value for value in idx1_nMagEqvHydrogens.values()]

    Active_range = 10
    dpi = 10000
    if not args.lw:
        args.lw = 1
    if not args.thr:
        args.thr = args.lw * 0.3

    return inHydrogen, Active_range, dpi


def _process_qm_hydrogen_spin_system(inParameter, idx0_ab_group_sets, mat_filter_multi, inAnmr, args) -> tuple[list[int], list[list[tuple[float, float]]]]:
    """Process QM hydrogen spin system data.

    Args:
        inParameter: Input parameters for the calculation
        idx0_ab_group_sets: Sets of AB group indices
        mat_filter_multi: Multiplier matrix filter
        inAnmr: Anmr object containing the NMR data
        args: Command line arguments

    Returns:
        tuple[list[int], list[list[tuple[float, float]]]]: Index range and peak lists
    """

    from censo_ext.Tools.qm import qm_base, qm_full, qm_multiplet
    accPeaks: list[list[tuple[float, float]]] = []
    inHydrogen: list[int]
    inSParams, inJCoups, inHydrogen = inParameter
    print("")
    print(" ===== Processing =====")
    print(" the group of calculate spectra :", len(idx0_ab_group_sets))
    print("  idx len(x) {x's AB quartet} {x's all - x's AB quartet} ")

    if len(inSParams*inHydrogen) <= args.mss:
        idx0_ab_group = list(idx0_ab_group_sets[0])
        v: npt.NDArray[np.float64] = inSParams[idx0_ab_group]
        J: npt.NDArray[np.float64] = inJCoups[idx0_ab_group].T[idx0_ab_group]

        for idx0, idx0_ab_group_set in enumerate(idx0_ab_group_sets):
            mat_multi_idx0: list[int] = mat_filter_multi[idx0].astype(
                int).tolist()
            idx0_ab_group: list[int] = list(idx0_ab_group_set)
            idx1_ab_group: set[int] = set(a+1 for a in idx0_ab_group_set)
            mat_multi_x_idx0: list[int] = [
                idx0_set*a for a, idx0_set in enumerate(mat_multi_idx0)if idx0_set != 0]
            print(f'{(idx0+1):>5d}{len(idx0_ab_group):>5d}', f'{idx1_ab_group}', set(
                a+1 for a in mat_multi_x_idx0).difference(idx1_ab_group))

        QM_Bases: list[tuple[float, float]] = qm_full(v=list(
            v), J=J, nIntergals=len(inHydrogen), args=args)
        from nmrsim.math import normalize_peaklist
        # QM_Multiplet = normalize_peaklist(QM_Base, len(inHydrogen))
        # accPeaks.append(QM_Multiplet)
        accPeaks.append(QM_Bases)

    else:
        for idx0, idx0_ab_group_set in enumerate(idx0_ab_group_sets):

            mat_multi_idx0: list[int] = mat_filter_multi[idx0].astype(
                int).tolist()
            idx0_ab_group: list[int] = list(idx0_ab_group_set)
            idx1_ab_group: set[int] = set(a+1 for a in idx0_ab_group_set)
            mat_multi_x_idx0: list[int] = [
                idx0_set*a for a, idx0_set in enumerate(mat_multi_idx0)if idx0_set != 0]
            print(f'{(idx0+1):>5d}{len(idx0_ab_group):>5d}', f'{idx1_ab_group}', set(
                a+1 for a in mat_multi_x_idx0).difference(idx1_ab_group))

            v: npt.NDArray[np.float64] = inSParams[idx0_ab_group]
            J: npt.NDArray[np.float64] = inJCoups[idx0_ab_group].T[idx0_ab_group]

            QM_Bases: list[tuple[float, float]] = qm_base(v=list(
                v), J=J, nIntergals=inHydrogen[idx0_ab_group.index(idx0)], idx0_nspins=idx0_ab_group.index(idx0), args=args)

            QM_Multiplet: list[tuple[float, float]] = []
            for QM_base in QM_Bases:
                idx0_multiplicity: list[int] = list(
                    set(mat_multi_x_idx0).difference(idx0_ab_group_set))
                # Chemical Shift, the numbers of Hydrogen in inJ
                inJCoups_multi: list[tuple[float, int]] = []
                for idx_m in idx0_multiplicity:
                    # this is not necessary, but the mat_multi_x_idx0 and idx_ab_group_set is OK
                    if np.fabs(inSParams[idx0]-inSParams[idx_m]) > 0.1:
                        inJCoups_multi.append(
                            (inJCoups[idx0][idx_m], inHydrogen[idx_m]))
                    else:
                        print("something wrong!!!!")
                        exit(0)

                if len(inJCoups_multi) >= 1:
                    tmp: npt.NDArray[np.float64] = np.array(
                        qm_multiplet(QM_base[0], nIntergals=1, J=inJCoups_multi))
                    tmp.T[1] *= QM_base[1]
                    QM_Multiplet += tmp.tolist()
                elif len(inJCoups_multi) == 0:
                    QM_Multiplet = QM_Bases
                else:
                    raise ValueError(
                        "  inJCoups_multi,  Exit and Close the program !!!")
            # normalize_peaklist is necessary,beacuse QM_Multiplet not only one
            from nmrsim.math import normalize_peaklist
            QM_Multiplet = normalize_peaklist(
                QM_Multiplet, inHydrogen[idx0])

            if len(accPeaks) == 0 and len(QM_Multiplet) == 0:
                pass
            else:
                accPeaks.append(QM_Multiplet)

    import json
    with open(inAnmr.get_Dir()/Path("peaks.json"), "w") as jsonFile:
        json.dump(accPeaks, jsonFile)

    idx0_peaks_range = [*range(len(accPeaks))]
    return idx0_peaks_range, accPeaks


def _process_qm_carbon_spin_system(inParameter, inAnmr) -> list[list[tuple[float, float]]]:
    """Process QM carbon spin system data.

    Args:
        inParameter: Input parameters for the calculation
        inAnmr: Anmr object containing the NMR data

    Returns:
        list[list[tuple[float, float]]]: List of peak lists
    """

    inSParams, _, inHydrogen = inParameter
    accPeaks: list[list[tuple[float, float]]] = []
    for idx0, ppm in enumerate(inSParams):
        dat: list = []
        dat.append((float(ppm*(-1)), float(inHydrogen[idx0])))
        accPeaks.append(dat)

    import json
    with open(inAnmr.get_Dir()/Path("peaks.json"), "w") as jsonFile:
        json.dump(accPeaks, jsonFile)
    return accPeaks


def _process_qm_json_spin_system(inAnmr, args) -> tuple[list[int], list[list[tuple[float, float]]]]:
    """Process QM JSON spin system data from ANMR.

    This function reads peak data from a JSON file and returns the appropriate
    range of peaks based on the provided arguments.

    Args:
        inAnmr: ANMR object containing directory information
        args: Arguments object containing json parameter

    Returns:
        tuple[list[int], list[list[tuple[float, float]]]]: A tuple containing
        the peak indices range and the actual peak data
    """

    accPeaks: list[list[tuple[float, float]]] = []
    import json
    with open(inAnmr.get_Dir()/Path("peaks.json"), "r") as jsonFile:
        accPeaks = json.load(jsonFile)

    if args.json[0] == -1:
        idx0_peaks_range = [*range(len(accPeaks))]
    else:
        idx0_peaks_range = args.json
    return idx0_peaks_range, accPeaks


def process_AB_quartet(inParameter: list[npt.NDArray[np.float64] | list[int]], inAnmr: Anmr, args: argparse.Namespace) \
        -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], list[int], list[set[int]], npt.NDArray[np.uint8]]:
    """Process AB quartet spin systems for NMR spectrum simulation.

    This function identifies and categorizes spin systems into AB quartets and multiplets
    based on chemical shift differences and coupling constants. It handles the complex
    logic of merging overlapping spin systems and adjusting thresholds automatically.

    Args:
        inParameter (list[npt.NDArray[np.float64] | list[int]]): Input parameters including
            chemical shifts, coupling constants, and hydrogen counts
        inAnmr (Anmr): ANMR object containing directory and magnetization information
        args (argparse.Namespace): Command line arguments with processing parameters

    Returns:
        tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], list[int], list[set[int]], npt.NDArray[np.uint8]]:
        - Chemical shifts array
        - Coupling constants matrix  
        - Hydrogen counts list
        - Spin system group sets
        - Multiplet filter matrix

    Example:
        >>> params = [shifts, couplings, hydrogens]
        >>> result = process_AB_quartet(params, anmr_obj, args)
    """

    # Initialize variables for AB quartet detection and processing
    idx0_ab_group_sets: list[set[int]] = []
    mat_filter_multi: npt.NDArray[np.uint8] = np.array([])

    # Main processing loop for identifying and categorizing spin systems
    inSParams: npt.NDArray[np.float64]
    inJCoups: npt.NDArray[np.float64]
    inHydrogen: list[int]
    inSParams, inJCoups, inHydrogen = inParameter  # type: ignore

    if not args.json and inAnmr.get_Anmr_Active()[0] == 'H':
        inJCoups_origin: npt.NDArray[np.float64] = copy.deepcopy(inJCoups)

        while (1):

            # the numbers of args.mss is low, all nucleus will be computated as AB quartet
            if len(inSParams*inHydrogen) <= args.mss:
                inJCoups = copy.deepcopy(inJCoups_origin)

                # Step 1: Filter out small coupling constants based on threshold
                # Delete Too Small JCoups J = args.lw*(-0.3) ~ args.lw*(0.3) use matrix Filter
                # 1: keep and 0: neglect
                # reset all inJCoups to zero
                mat_filter_low_factor: npt.NDArray[np.uint8] = np.zeros(
                    (inSParams.size, inSParams.size), dtype=np.uint8)
                mat_filter_low_factor = (
                    np.abs(inJCoups) > args.thr).astype(np.uint8)
                inJCoups[np.logical_not(mat_filter_low_factor)] = 0
                del mat_filter_low_factor

                # Step 2: add more one hydrogen in the same position as the first hydrogen
                idx_repeat = (np.array(inHydrogen)-1).nonzero()[0]
                for idx in idx_repeat[::-1]:
                    for idy in range(1, inHydrogen[idx], 1):
                        inSParams = np.insert(inSParams, idx, inSParams[idx])
                        inJCoups = np.insert(inJCoups, idx, inJCoups[idx], axis=1)  # nopep8
                        inJCoups = np.insert(inJCoups, idx, inJCoups[idx], axis=0)  # nopep8
                inHydrogen = [1]*len(inSParams)
                mat_filter_low_factor = np.ones(
                    (inSParams.size, inSParams.size), dtype=np.uint8)
                np.fill_diagonal(mat_filter_low_factor, 0)
                mat_filter_ab_quartet: npt.NDArray[np.uint8] = copy.deepcopy(
                    mat_filter_low_factor)
            else:
                # Step 1: Filter out too small coupling constants based on threshold
                # Delete too small JCoups J = args.lw*(-0.3) ~ args.lw*(0.3) use matrix Filter
                # 1: keep and 0: neglect
                mat_filter_low_factor: npt.NDArray[np.uint8] = np.zeros(
                    (inSParams.size, inSParams.size), dtype=np.uint8)
                inJCoups = copy.deepcopy(inJCoups_origin)
                mat_filter_low_factor = (
                    np.abs(inJCoups) > args.thr).astype(np.uint8)
                inJCoups[np.logical_not(mat_filter_low_factor)] = 0

                # Step 2: Identify potential AB quartet systems
                mat_filter_ab_quartet: npt.NDArray[np.uint8] = np.zeros(
                    (inSParams.size, inSParams.size), dtype=np.uint8)
                import math
                for idx0, x in enumerate(inSParams):
                    for idy0, y in enumerate(inSParams):

                        # Check if two chemical shifts are very close (AB quartet condition)
                        # if two chemical shift is very close, will perform AB quartet
                        # if x-y == 0 the Ratio_J_Hz will crash
                        # if the JCoups is negative, the peaks will prioritize to use normal QM calculation.
                        # if the peaks which is not in AB Quartet and have positive nubmers will use Multiplet (nmrsim)
                        # Normal the maximum of 3-JCoups is 18 Hz. If JCoups is set to 20 Hz, the delta Chemical Shift is set to 1.0 ppm
                        # it will have 0.0004 ppm < 0.001 ppm (lw = 1)
                        # If diff of two chemical shift AB quart is 0.006 ppm and the JCouping of J=20 Hz, it will produce 0.22 Hz
                        # limits_min_ab = 0.005 is enough for almost condtions

                        limits_min_ab: float = 0.005
                        if (math.fabs(x-y) < limits_min_ab and mat_filter_low_factor[idx0][idy0] == 1) or (inJCoups[idx0][idy0] <= -args.thr):
                            mat_filter_ab_quartet[idx0][idy0] = 1
                        else:
                            if math.fabs(x-y) < limits_min_ab:
                                Ratio_J_Hz: float = 10000
                            else:
                                Ratio_J_Hz: float = math.fabs(
                                    inJCoups[idx0][idy0]/(math.fabs(x-y)))
                            if Ratio_J_Hz < args.thrab:
                                mat_filter_ab_quartet[idx0][idy0] = 0
                            elif Ratio_J_Hz >= args.thrab and mat_filter_low_factor[idx0][idy0] == 1:
                                mat_filter_ab_quartet[idx0][idy0] = 1
                            else:
                                if mat_filter_low_factor[idx0][idy0] == 1:
                                    raise ValueError(
                                        f"{idx0} {x} {idy0} {y} was not found or is a directory")

            if args.verbose:
                ic(mat_filter_ab_quartet)

            # Calculate which couplings are NOT part of AB quartets (multiplets)
            mat_filter_multi = mat_filter_low_factor - mat_filter_ab_quartet

            # the atom connect relation of AB quaret
            idx0_ab_connect: list[list[int | set[int]]] = []
            for idx0, x in enumerate(mat_filter_ab_quartet):
                group: set[int] = set(np.array(x*(idx0+1)).nonzero()[0].tolist())  # return arg # nopep8 #idx0+1 is only for nozero
                group.add(idx0)
                idx0_ab_connect.append([idx0, group])
            if args.verbose:
                ic(idx0_ab_connect)

            # Merge overlapping spin system groups
            idx0_ab_group_sets = []
            for _, x in idx0_ab_connect:
                if len(x) == 1:                      # type: ignore
                    idx0_ab_group_sets.append(x)     # type: ignore
                elif len(x) > 1:                     # type: ignore
                    group: set[int] = x              # type: ignore
                    loop: bool = True
                    bond_penetration: int = 1
                    while loop:
                        loop = False
                        for y in group:
                            if (not group.issuperset(idx0_ab_connect[y][1])) and bond_penetration <= args.tb:  # type: ignore # nopep8
                                loop = True
                                group = group.union(idx0_ab_connect[y][1])  # type: ignore # nopep8
                            bond_penetration += 1
                    idx0_ab_group_sets.append(group)
                else:
                    raise ValueError("  idx0_ab_group_sets is bugs !!!")

            if len(inSParams*inHydrogen) > args.mss:
                # Handle CH3 equivalent groups manually (symmetry considerations)
                # So if chemical shift in AB quartet region need to move to multiplet
                list_Equivalent3: list[int] = []
                for key in inAnmr.nMagnetEqvs.keys():
                    if inAnmr.nMagnetEqvs[key] == 3:
                        for idy0, y in enumerate(inAnmr.avg_orcaSJ.idx1Atoms):
                            if y == min(inAnmr.NeighborMangetEqvs[key]):
                                list_Equivalent3.append(idy0)
                set_Equivalent3: set[int] = set(list_Equivalent3)
                del list_Equivalent3

                # Adjust groups to account for equivalent protons
                # Equivalent3 is idx0 numbers
                for idx0, group_set in enumerate(idx0_ab_group_sets):
                    set_move: set[int] = group_set.intersection(
                        set_Equivalent3)
                    if not len(set_move) == 0:
                        idx0_ab_group_sets[idx0] = set(group_set).difference(
                            set_move).union(set([idx0]))
                        set_move = set_move.difference(set([idx0]))
                    for y in set_move:
                        mat_filter_multi[idx0][y] = 1

            if args.verbose:
                ic(idx0_ab_group_sets)

            #  show every step of threshold
            if args.verbose:
                print(" ===== Processing =====")
                print(f" threshold of JCoupling  : {args.thr:>3.5f}")
                print(f" threshold of AB quartet : {args.thrab:>3.5f}")
                print(
                    "  idx len(x) {x's AB quartet} {x's all - x's AB quartet} ")

            # calculation the maximum length of the slice of AB quartet
            max_len_AB: int = 0
            for idx0, group_set in enumerate(idx0_ab_group_sets):
                mat_multi_idx0: list[int] = [
                    idx0_set*x for x, idx0_set in enumerate(mat_filter_multi[idx0].tolist())if idx0_set != 0]
                if args.verbose:
                    print(f'{(idx0+1):>5d}{len(group_set):>5d}', {a+1 for a in group_set}, set(
                        a+1 for a in mat_multi_idx0).difference({a+1 for a in group_set}))
                if len(group_set) > max_len_AB:
                    max_len_AB = len(group_set)

            # Adjust thresholds automatically if maximum spin system exceeds limit
            if (max_len_AB > args.mss):
                if args.auto:
                    args.thrab = args.thrab + 0.0025
                    args.thr = args.thr + 0.030
                else:
                    print("  Need to tidy the nspins of AB quartet and use cmd -thrab ")
                    print("  Exit and Close the program !!!")
                    exit(0)
            else:
                break

        # Additional processing for AB quartets with identical chemical shifts
        # AB quartet if more than two peaks of AB quartet, added closed peaks (not in AB quartet)
        print(" ===== Modification AB quartet =====")
        for idx0, idx0_ab_group_set in enumerate(idx0_ab_group_sets):
            idx0_ab_group: list[int] = list(idx0_ab_group_set)
            if len(set(list(inSParams[idx0_ab_group]))) == 1:
                from censo_ext.Tools.spectra import find_nearest
                _, Move_idx0 = find_nearest(list(inSParams[idx0_ab_group]),
                                            inSParams[idx0_ab_group[0]])
                arg = np.argwhere(
                    inSParams[:] == inSParams[int(Move_idx0)])[0]
                idx0_ab_group_sets[idx0] = idx0_ab_group_set.union(
                    set(int(x) for x in arg))

        # Display the parameter of Full Spectra
        for idx0, idx0_ab_group_set in enumerate(idx0_ab_group_sets):
            mat_multi_idx0: list[int] = mat_filter_multi[idx0].astype(
                int).tolist()
            idx0_ab_group: list[int] = list(idx0_ab_group_set)
            idx1_ab_group: set[int] = set(a+1 for a in idx0_ab_group_set)
            mat_multi_x_idx0: list[int] = [
                idx0_set*x for x, idx0_set in enumerate(mat_multi_idx0)if idx0_set != 0]
            print(f'{(idx0+1):>5d}{len(idx0_ab_group):>5d}', f'{idx1_ab_group}', set(
                a+1 for a in mat_multi_x_idx0).difference(idx1_ab_group))
        print(" Use this parameter to calculate the Full Spectra")
    return inSParams, inJCoups, inHydrogen, idx0_ab_group_sets, mat_filter_multi


def generate_final_spectrum(finalPeaks: list[tuple[float, float]], inAnmr: Anmr, dpi: int | None,
                            Active_range: int | None, args: argparse.Namespace) -> npt.NDArray[np.float64]:
    """Generate and plot the final spectrum from peaks data.

    This function processes peak data to create a spectrum plot. It handles
    the plotting configuration based on provided parameters and returns the
    resulting spectrum array.

    Args:
        finalPeaks (list[tuple[float, float]]): List of peak data as (frequency, intensity) tuples
        inAnmr (Anmr): ANMR object containing directory information
        dpi (int | None): DPI setting for the plot output
        Active_range (int | None): Active range parameter for plotting
        args (argparse.Namespace): Command line arguments namespace

    Returns:
        npt.NDArray[np.float64]: Array containing the final spectrum data

    Example:
        >>> peaks = [(100.0, 0.5), (200.0, 0.8)]
        >>> result = generate_final_spectrum(peaks, anmr_obj, 300, 10, args)
    """

    print("")
    print(" ===== Processing to plotting spectra =====")
    print(" Wait a minutes ...")
    args.out = str(inAnmr.get_Dir()/Path(args.out))

    if dpi and Active_range:
        print(" All done ...")
        from censo_ext.Tools.qm import print_plot
        return print_plot(in_plist=finalPeaks, dpi=dpi, args=args, Active_range=Active_range)
    else:
        print("  dpi and Active_range is wrong")
        print("  Exit and Close the program !!!")
        exit(0)


def main(args: argparse.Namespace = argparse.Namespace()) -> npt.NDArray[np.float64]:

    if args == argparse.Namespace():
        args = cml()

    # Setup
    inAnmr: Anmr = setup_anmr(args=args)

    # Handle average data loading if specified
    process_average_data(inAnmr=inAnmr, args=args)

    # preprocessing_spin_system
    inParameter: list[npt.NDArray[np.float64] | list[int]]
    *inParameter, Active_range, dpi, args = preprocess_spin_system(inAnmr=inAnmr, args=args)

    if args.verbose:
        ic(*inParameter)
    *inParameter, idx0_ab_group_sets, mat_filter_multi = process_AB_quartet(
        inParameter, inAnmr=inAnmr, args=args)

    idx0_peaks_range: list[int] = []
    accPeaks: list[list[tuple[float, float]]]
    if not args.json and inAnmr.get_Anmr_Active()[0] == 'H':
        idx0_peaks_range, accPeaks = _process_qm_hydrogen_spin_system(
            inParameter, idx0_ab_group_sets, mat_filter_multi, inAnmr, args)
    elif not args.json and inAnmr.get_Anmr_Active()[0] == 'C':
        accPeaks = _process_qm_carbon_spin_system(inParameter, inAnmr)
    elif args.json:
        idx0_peaks_range, accPeaks = _process_qm_json_spin_system(inAnmr, args)
    else:
        raise ValueError("  Something Wrong in your get_anmr_Active()")

    finalPeaks: list[tuple[float, float]] = []
    if inAnmr.get_Anmr_Active()[0] == 'H':
        for idx0, peak in enumerate(accPeaks):
            if idx0 in idx0_peaks_range:
                finalPeaks += peak
    elif inAnmr.get_Anmr_Active()[0] == 'C':
        for peak in accPeaks:
            finalPeaks += peak
    else:
        raise ValueError("  Something Wrong in your get_anmr_Active()")

    return generate_final_spectrum(finalPeaks=finalPeaks, inAnmr=inAnmr, dpi=dpi,
                                   Active_range=Active_range, args=args)


if __name__ == "__main__":
    main()
