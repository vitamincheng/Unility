#!/usr/bin/env python
# ToDo
# CH3 Equiv problem
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
from censo_ext.Tools.utility import print_arguments
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

    args: argparse.Namespace = parser.parse_args()
    return args


def main(args: argparse.Namespace = argparse.Namespace()) -> None:

    if args == argparse.Namespace():
        args = cml()

    print_arguments()

    inFile = Path(args.file)
    xyzFile: GeometryXYZs = GeometryXYZs(inFile)
    xyzFile.method_read_xyz()

    from censo_ext.Tools.anmrfile import Anmr
    inAnmr: Anmr = Anmr(Dir=args.dir)
    inAnmr.method_read_anmrrc()
    inAnmr.avg_Data_AD = AD_Normal(Dir=args.dir)
    inAnmr.get_avg_orcaSJ_Exist()
    inAnmr.method_BOBYQA_load_avg_orcaSJ()
    inAnmr.avg_orcaSJ.method_load_anmrrc_linear(inAnmr.get_Anmrrc_linear())
    inAnmr.avg_orcaSJ.method_setup_ChemicalShifts()
    inAnmr.avg_orcaSJ.method_print_av_orcaS()
    args.mf = 500
    inSParams: dict[AtomID, float] = inAnmr.avg_orcaSJ.ChemicalShifts.copy()
    inAnmr.avg_orcaSJ.method_teardown_ChemicalShifts()
    Element: dict[AtomID, str] = inAnmr.avg_orcaSJ.Element.copy()
    idx0_Element: list[AtomID] = [x for x in Element.keys()]

    inAnmr.method_read_enso()
    # ic(inAnmr.enso['ONOFF'])
    # ic(inAnmr.enso['BW'])
    H_Atoms: list[AtomID] = [key for key,
                             value in Element.items() if value == "H"]
    nShapes: int = len(H_Atoms)
    idx0_Atoms: list = list(np.array(H_Atoms)-1)
    Distances: npt.NDArray[np.float64] = np.zeros(
        (nShapes, nShapes), dtype=np.float64)
    Result: npt.NDArray[np.float64] = np.zeros(
        (nShapes, nShapes), dtype=np.float64)
    for idx0, St in enumerate(xyzFile.Sts):
        St_local = np.array(St.coord)[idx0_Atoms]
        from scipy.spatial.distance import cdist
        Distances = cdist(St_local, St_local)  # type: ignore
        np.seterr(divide='ignore')
        Result += inAnmr.enso['BW'][idx0]*(10**6) / (Distances ** 6)
        np.seterr(divide='warn')

    diag_indices = np.diag_indices_from(Result)
    Result[diag_indices] = 0
    data_x = Load_Directory(args)

    if args.start is None or args.end is None:
        start = np.min(data_x.T[0])
        end = np.max(data_x.T[0])
    else:
        start = float(args.start)
        end = float(args.end)

    # Create a grid of x and y values
    x: npt.NDArray = np.linspace(start, end, args.pts)
    y: npt.NDArray = np.linspace(start, end, args.pts)
    X, Y = np.meshgrid(x, y)
    Z: npt.NDArray

    # Define Lorentzian parameters
    gamma_x = args.gamma
    gamma_y = args.gamma

    # Calculate the 2D Lorentzian
    for idx0, x in enumerate(Result):
        for idy0, amp in enumerate(x):
            amplitude: float = amp
            try:
                Z += lorentzian_2d(X, Y, amplitude, inSParams[idx0_Element[idx0]], inSParams[idx0_Element[idy0]],  # type: ignore
                                   gamma_x, gamma_y)
            except NameError:
                Z = lorentzian_2d(X, Y, amplitude, inSParams[idx0_Element[idx0]], inSParams[idx0_Element[idy0]],
                                  gamma_x, gamma_y)

    # Plotting the result
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
    min_contour: int = int(args.contour)*2**6
    int_ratio_2: int = math.ceil(math.log2(np.max(Result)/min_contour))
    init = np.arange(int_ratio_2+1)
    lv = (np.power(2, init)*min_contour)
    # _plot = ax.contour(X, Y, Z, levels=lv, cmap=matplotlib.cm.Blues_r)
    _plot = ax.contour(X, Y, Z, levels=lv, cmap='seismic')  # type: ignore
    ax.set_xlim(start, end)
    ax.set_ylim(start, end)
    ax_histx.set_xlim(start, end)
    ax_histy.set_ylim(start, end)
    ax.invert_xaxis()
    ax.invert_yaxis()
    ax.clabel(_plot, fontsize=6)

    # plt.colorbar(_plot, ax=ax, label='Intensity')
    plt.show()


def Load_Directory(args) -> npt.NDArray:
    import censo_ext.anmr as anmr
    import os
    import sys
    in_dir = args.dir
    args.average = True
    dir_H = in_dir
    args_x: dict = {"auto": True, "average": args.average, "bobyqa": False, "mf": 500,
                    "dir": dir_H, "thr": None, "json": None, "thrab": 0.025,
                    "verbose": False, "lw": 1, "tb": 4, "mss": 10, "cutoff": 0.001,
                    "show": False, "start": None, "end": None, "out": "output.npz"}
    sys.stdout = open(os.devnull, 'w')
    data_x: npt.NDArray[np.float64] = anmr.main(argparse.Namespace(**args_x)).T
    return data_x


def lorentzian_2d(x: npt.NDArray, y: npt.NDArray, amplitude: float, center_x: float, center_y: float, gamma_x: float, gamma_y: float, rotation_angle=0) -> npt.NDArray[np.float64]:
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
    xp: npt.NDArray = (x - center_x) * np.cos(rotation_angle) - \
        (y - center_y) * np.sin(rotation_angle)
    yp: npt.NDArray = (x - center_x) * np.sin(rotation_angle) + \
        (y - center_y) * np.cos(rotation_angle)

    # Calculate the Lorentzian
    return amplitude / (1 + (xp / gamma_x)**2 + (yp / gamma_y)**2)


if __name__ == "__main__":
    main()

#   test
#   python3 molclus_xtb.py -i ../tests/crest_conformers.xyz --alpb CHCl3
