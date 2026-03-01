#!/usr/bin/env python
#
from matplotlib.axes import Axes
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
from censo_ext.Tools.anmrfile import AD_Normal
from censo_ext.Tools.xyzfile import GeometryXYZs
import argparse
import numpy as np
import numpy.typing as npt
# from icecream import ic
from pathlib import Path
from censo_ext.Tools.utility import IntpID, print_arguments
from censo_ext.anmr import AtomID
from matplotlib.figure import Figure
descr = """
________________________________________________________________________________
| plot_2D_ditance.py   
| Usages   : plot_2D_distance.py <geometry> [options]
| Input    : -i input file [default crest_conformers.xyz]
| Dir      : -d input directory [default .]
|______________________________________________________________________________
"""


def cml() -> argparse.Namespace:
    """ Get args object from commandline interface.
        Needs argparse module."""
    parser = argparse.ArgumentParser(
        # description=f"{descr}",
        # formatter_class=argparse.RawDescriptionHelpFormatter,
        usage=argparse.SUPPRESS,
    )

    parser.add_argument(
        "-i",
        "--input",
        dest="file",
        action="store",
        required=False,
        default="crest_conformers.xyz",
        help="Provide input_file name [default crest_conformers.xyz]",
    )

    parser.add_argument(
        "-d",
        "--dir",
        dest="dir",
        action="store",
        required=False,
        default=".",
        help="Provide output_file name [default .]",
    )

    parser.add_argument(
        "-c",
        "--contour",
        dest="contour",
        action="store",
        required=False,
        default=4,
        help="Provide minimum of contour [default 4]",
    )

    parser.add_argument(
        "-g",
        "--gamma",
        dest="gamma",
        action="store",
        required=False,
        default=0.01,
        help="Provide gamma parameter of lorentzian [default 0.01]",
    )

    parser.add_argument(
        "-p",
        "--pts",
        dest="pts",
        action="store",
        required=False,
        default=1024,
        help="Provide the points of 2D spectra [default 1024]",
    )

    parser.add_argument(
        "-start",
        dest="start",
        action="store",
        required=False,
        default=None,
        help="Provide start of ppm [default from data]",
    )
    parser.add_argument(
        "-end",
        dest="end",
        action="store",
        required=False,
        default=None,
        help="Provide end of ppm [default from data]",
    )

    return parser.parse_args()


def main(args: argparse.Namespace = argparse.Namespace()) -> None:

    if args == argparse.Namespace():
        args = cml()

    print_arguments()

    xyzFile: GeometryXYZs = GeometryXYZs(Path(args.file))
    xyzFile.method_read_xyz()

    names_str: dict[AtomID, str] = xyzFile.Sts[0].names
    names_xyzFile: list[AtomID] = [key for key,
                                   value in names_str.items() if value == "H"]
    # ic(names_xyzFile)
    intp_Atoms: list[IntpID] = list(np.array(names_xyzFile)-1)
    nShapes: int = len(intp_Atoms)
    from censo_ext.Tools.anmrfile import Anmr
    inAnmr: Anmr = Anmr(Dir=args.dir)
    inAnmr.method_read_anmrrc()
    inAnmr.method_read_nucinfo()
    # inAnmr.method_print_nucinfo()
    inAnmr.avg_Data_AD = AD_Normal(Dir=args.dir)
    inAnmr.get_avg_orcaSJ_Exist()
    inAnmr.method_BOBYQA_load_avg_orcaSJ()
    inAnmr.avg_orcaSJ.method_load_anmrrc_linear(inAnmr.get_Anmrrc_linear())
    inAnmr.avg_orcaSJ.method_setup_ChemicalShifts()
    inAnmr.avg_orcaSJ.method_print_av_orcaS()
    args.mf = 500
    inSParams: dict[AtomID, float] = inAnmr.avg_orcaSJ.ChemicalShifts.copy()
    inAnmr.avg_orcaSJ.method_teardown_ChemicalShifts()
    Elements_str: dict[AtomID, str] = inAnmr.avg_orcaSJ.Element.copy()
    Elements: list[AtomID] = [x for x in Elements_str.keys()]

    inAnmr.method_read_enso()
    # ic(inAnmr.enso['ONOFF'])
    # ic(inAnmr.enso['BW'])
    atomIDs_H_AD: list[AtomID] = [key for key,
                                  value in Elements_str.items() if value == "H"]
    nShapes_AD: int = len(atomIDs_H_AD)
    intp_Atoms_AD: list[IntpID] = (np.array(atomIDs_H_AD)-1).tolist()
    Result_AD: npt.NDArray[np.float64] = np.zeros(
        (nShapes_AD, nShapes_AD), dtype=np.float64)
    for idx0, St in enumerate(xyzFile.Sts):
        St_local = np.array(St.coord)[intp_Atoms]
        from scipy.spatial.distance import cdist
        Distances: npt.NDArray[np.float64] = np.zeros(
            (nShapes, nShapes), dtype=np.float64)
        Distances = cdist(St_local, St_local)  # type: ignore

        np.seterr(divide='ignore')
        Distances = (10**6) / (Distances ** 6)
        np.seterr(divide='warn')

        diag_indices = np.diag_indices_from(Distances)
        Distances[diag_indices] = 0
        # ic(Distances)
        # ic(Result_AD)
        # ic(inAnmr.NeighborMangetEqvs)
        # ic(len(intp_Atoms))
        # ic(len(intp_Atoms_AD))
        # ic(intp_Atoms)
        for idx0, atomid_AD in enumerate(intp_Atoms_AD):
            # print(idx0, atomid_AD)
            if atomid_AD in intp_Atoms:
                # print(idx0, x+1, intp_Atoms_AD.index(int(x)))
                if len(inAnmr.NeighborMangetEqvs[AtomID(int(atomid_AD)+1)]) != 1:
                    # print(inAnmr.NeighborMangetEqvs[AtomID(int(x)+1)])
                    list_AtomsID: list[AtomID] = inAnmr.NeighborMangetEqvs[AtomID(
                        atomid_AD+1)]
                    list_intpID = [intp_Atoms.index(
                        IntpID(x-1)) for x in list_AtomsID]
                    Distances[list_intpID] = np.average(
                        Distances[list_intpID], axis=0)
                    Distances.T[list_intpID] = np.average(
                        Distances.T[list_intpID], axis=0)
                # else:
                #    print(inAnmr.NeighborMangetEqvs[AtomID(int(x)+1)])
            else:
                print(" something wrong in your data")
                exit(0)
        wait_remove_AtomsID: list[IntpID] = [
            x for x in intp_Atoms if x not in intp_Atoms_AD]
        wait_remove_AtomsID = sorted(wait_remove_AtomsID, reverse=True)
        # print(wait_remove_AtomsID)
        for a in wait_remove_AtomsID:
            z = intp_Atoms.index(a)
            Distances = np.delete(Distances, z, axis=1)
            Distances = np.delete(Distances, z, axis=0)

        # ic(Distances)
        Result_AD += inAnmr.enso['BW'][idx0]*Distances

    diag_indices = np.diag_indices_from(Result_AD)
    Result_AD[diag_indices] = 0
    data_x: npt.NDArray[np.float64] = Load_Directory(args)

    if args.start is None or args.end is None:
        start = np.min(data_x.T[0])
        end = np.max(data_x.T[0])
    else:
        start = float(args.start)
        end = float(args.end)

    # Create a grid of x and y values
    atomsID_AD: npt.NDArray[np.float64] = np.linspace(
        start, end, args.pts).astype(np.float64)
    y: npt.NDArray[np.float64] = np.linspace(
        start, end, args.pts).astype(np.float64)
    X, Y = np.meshgrid(atomsID_AD, y)
    Z: npt.NDArray[np.float64]

    # Define Lorentzian parameters
    # Calculate the 2D Lorentzian
    gamma_x: float = args.gamma
    gamma_y: float = args.gamma

    for idx0, atomsID_AD in enumerate(Result_AD):
        for idy0, amp in enumerate(atomsID_AD):
            amplitude: float = amp
            try:
                Z += lorentzian_2d(X, Y, amplitude, inSParams[Elements[idx0]], inSParams[Elements[idy0]],  # type: ignore
                                   gamma_x, gamma_y)
            except NameError:
                Z = lorentzian_2d(X, Y, amplitude, inSParams[Elements[idx0]], inSParams[Elements[idy0]],
                                  gamma_x, gamma_y)

    # Plotting the result
    plot_diagram(args.contour, Result_AD, data_x,
                 start, end, X, Y, Z)  # type: ignore

    # plt.colorbar(_plot, ax=ax, label='Intensity')
    plt.show()


def plot_diagram(contour: float, Result: npt.NDArray[np.float64], data_x: npt.NDArray[np.float64],
                 start: float, end: float, X: npt.NDArray[np.float64], Y: npt.NDArray[np.float64],
                 Z: npt.NDArray[np.float64]) -> None:
    fig: Figure = plt.figure(figsize=(11.7, 8.3), dpi=100)
    gs: GridSpec = fig.add_gridspec(2, 2,  width_ratios=(1, 19), height_ratios=(1, 9),
                                    left=0.03, right=0.97, bottom=0.03, top=0.97,
                                    wspace=0.1, hspace=0.1)
    ax: Axes = fig.add_subplot(gs[1, 1])

    ax_histx: Axes = fig.add_subplot(gs[0, 1], sharex=ax)
    ax_histy: Axes = fig.add_subplot(gs[1, 0], sharey=ax)
    ax_histx.get_xaxis().set_visible(False)
    ax_histx.get_yaxis().set_visible(False)
    ax_histx.axis('off')
    ax_histy.get_xaxis().set_visible(False)
    ax_histy.get_yaxis().set_visible(False)
    ax_histy.axis('off')

    x_axis_data: npt.NDArray[np.float64] = data_x.T[1]
    y_axis_data: npt.NDArray[np.float64] = data_x.T[1]

    ax_histx.plot(data_x.T[0], x_axis_data)
    ax_histy.plot(-y_axis_data, data_x.T[0])
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')
    fig.suptitle(r"$10^{6}$ / $r^{6}$", fontsize=12, x=0.10, y=0.98)

    import math
    min_contour: int = int(contour)*2**6
    int_ratio_2: int = math.ceil(math.log2(np.max(Result)/min_contour))
    init: npt.NDArray[np.int64] = np.arange(int_ratio_2+1).astype(np.int64)
    lv: npt.NDArray[np.float64] = (
        np.power(2, init)*min_contour).astype(np.float64)
    # _plot = ax.contour(X, Y, Z, levels=lv, cmap=matplotlib.cm.Blues_r)
    _plot = ax.contour(X, Y, Z, levels=lv, cmap='seismic')  # type: ignore
    ax.set_xlim(start, end)
    ax.set_ylim(start, end)
    ax_histx.set_xlim(start, end)
    ax_histy.set_ylim(start, end)
    ax.invert_xaxis()
    ax.invert_yaxis()
    ax.clabel(_plot, fontsize=6)


def Load_Directory(args: argparse.Namespace) -> npt.NDArray[np.float64]:
    import censo_ext.anmr as anmr
    import os
    import sys
    dir_H: str = args.dir
    args.average = True
    args_x: dict = {"auto": True, "average": args.average, "bobyqa": False, "mf": 500,
                    "dir": dir_H, "thr": None, "json": None, "thrab": 0.025,
                    "verbose": False, "lw": 1, "tb": 4, "mss": 10, "cutoff": 0.001,
                    "show": False, "start": None, "end": None, "out": "output.npz"}
    sys.stdout = open(os.devnull, 'w')
    data_x: npt.NDArray[np.float64] = anmr.main(argparse.Namespace(**args_x)).T
    return data_x


def lorentzian_2d(x: npt.NDArray[np.float64], y: npt.NDArray[np.float64], amplitude: float, center_x: float, center_y: float, gamma_x: float, gamma_y: float, rotation_angle=0) -> npt.NDArray[np.float64]:
    """
    Calculates a 2D Lorentzian function.

    Parameters:
        x (numpy.ndarray): x-coordinates.
        y (numpy.ndarray): y-coordinates.
        amplitude (float): The peak amplitude.
        center_x (float): The x-coordinate of the peak center.
        center_y (float): The y-coordinate of the peak center.
        gamma_x (float): The half-width at half-maximum (HWHM) in the x-direction.
        gamma_y (float): The half-width at half-maximum (HWHM) in the y-direction.
        rotation_angle (float, optional): Rotation angle in radians. Defaults to 0.

    Returns:
        numpy.ndarray: The 2D Lorentzian values.
    """
    # Translate and rotate coordinates
    xp: npt.NDArray[np.float64] = (x - center_x) * np.cos(rotation_angle) - \
        (y - center_y) * np.sin(rotation_angle)
    yp: npt.NDArray[np.float64] = (x - center_x) * np.sin(rotation_angle) + \
        (y - center_y) * np.cos(rotation_angle)

    # Calculate the Lorentzian
    return amplitude / (1 + (xp / gamma_x)**2 + (yp / gamma_y)**2)


if __name__ == "__main__":
    main()

#   test
#   python3 molclus_xtb.py -i ../tests/crest_conformers.xyz --alpb CHCl3
