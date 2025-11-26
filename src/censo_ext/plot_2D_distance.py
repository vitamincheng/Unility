#!/usr/bin/env python
from matplotlib.axes import Axes
import matplotlib.pyplot as plt
from censo_ext.Tools.anmrfile import AD_Normal
from censo_ext.Tools.xyzfile import GeometryXYZs
import argparse
import numpy as np
import numpy.typing as npt
from icecream import ic
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
    args: argparse.Namespace = parser.parse_args()
    return args


def main(args: argparse.Namespace = argparse.Namespace()) -> None:

    if args == argparse.Namespace():
        args = cml()

    print_arguments()

    inFile = Path(args.file)
    xyzFile: GeometryXYZs = GeometryXYZs(inFile)
    xyzFile.method_read_xyz()
    # xyzFile.method_print(idx1_St=[])

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
    inSParams = inAnmr.avg_orcaSJ.ChemicalShifts.copy()
    inAnmr.avg_orcaSJ.method_teardown_ChemicalShifts()
    Element = inAnmr.avg_orcaSJ.Element.copy()
    idx0_Element = [x for x in Element.keys()]

    ic(inSParams)
    ic(Element)
    ic(len(Element))
    inAnmr.method_read_enso()
    ic(inAnmr.enso['ONOFF'])
    ic(inAnmr.enso['BW'])
    H_Atoms: list[AtomID] = [key for key,
                             value in Element.items() if value == "H"]
    ic(H_Atoms)
    nShapes = len(H_Atoms)
    ic(nShapes)
    idx0_Atoms = list(np.array(H_Atoms)-1)
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

    # for x in Result:
    #    print(x)
    # Your code with potential division by zero

    # Create a grid of x and y values
    x: npt.NDArray = np.linspace(0, 8, 1024)
    y: npt.NDArray = np.linspace(0, 8, 1024)
    X, Y = np.meshgrid(x, y)
    Z: npt.NDArray

    # Define Lorentzian parameters
    gamma_x = 0.03
    gamma_y = 0.03

    # Calculate the 2D Lorentzian
    for idx0, x in enumerate(Result):
        for idy0, amp in enumerate(x):
            amplitude: float = amp
            # print(amplitude, inSParams[idx0_Element[idx0]],
            #      inSParams[idx0_Element[idy0]])
            try:
                Z += lorentzian_2d(X, Y, amplitude, inSParams[idx0_Element[idx0]], inSParams[idx0_Element[idy0]],
                                   gamma_x, gamma_y)
            except NameError:
                Z = lorentzian_2d(X, Y, amplitude, inSParams[idx0_Element[idx0]], inSParams[idx0_Element[idy0]],
                                  gamma_x, gamma_y)

    # Plotting the result
    _fig: Figure = plt.figure(figsize=(11.7, 8.3), dpi=100)
    ax: Axes = plt.subplot()
    #  import matplotlib.ticker as ticker
    # CS2 = plt.contour(X, Y, Z, locator=plt.LogLocator())
    # fmt = ticker.LogFormatterMathtext()
    # fmt.create_dummy_axis()
    # plt.clabel(CS2, CS2.levels, fmt=fmt)
    lv = np.linspace(np.min(Result), np.max(Result), 10)
    # plt.contourf(X, Y, Z, levels=lv, cmap='viridis')
    contour_plot = ax.contourf(
        X, Y, Z, levels=lv, cmap='viridis', extend='both')
    ax.contour(X, Y, Z, levels=lv)
    plt.colorbar(contour_plot, ax=ax, label='Intensity')
    ax.set_title('2D Lorentzian Function')
    ax.set_xlabel('X-axis')
    ax.set_ylabel('Y-axis')
    plt.show()


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
