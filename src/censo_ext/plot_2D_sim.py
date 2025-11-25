#!/usr/bin/env python
import argparse
from matplotlib.axes import Axes
from matplotlib.text import Text
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec
from censo_ext.Tools.utility import AtomID, print_arguments
import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
import os
import sys
from pathlib import Path
descr = """
________________________________________________________________________________
| For plot_2D_Sim.py
| Usages    : plot_2D_sim.py <geometry> [options]
| [options]
| Directory : -d two input directory folder [required] 
| Proton    : -p limits of proton spectra [default from data]
| Carbon    : -c limits of carbon spectra [default from data]
|______________________________________________________________________________
"""


def cml() -> argparse.Namespace:
    """ Get args object from commandline interface. Needs argparse module."""
    parser = argparse.ArgumentParser(
        description="descr",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage=argparse.SUPPRESS)
    parser.add_argument(
        "-d",
        "--dir",
        dest="dir",
        action="store",
        type=str,
        required=True,
        nargs=2,
        help="Provide two input directory folder [dir of Hydrogen, dir of Carbon]",
    )
    parser.add_argument(
        "-av",
        "--average",
        dest="average",
        action="store_true",
        help="Use the Average Directory of nmr [default False]",
    )
    parser.add_argument(
        "-p",
        "--proton",
        dest="h_limits",
        action="store",
        required=False,
        default=None,
        nargs=2,
        type=float,
        help="Start plotting from '<start>' ppm and End plotting from '<end>' ppm in H spectra",
    )
    parser.add_argument(
        "-c",
        "--carbon",
        dest="c_limits",
        action="store",
        required=False,
        default=None,
        nargs=2,
        type=float,
        help="Start plotting from '<start>' ppm and End plotting from '<end>' ppm in C spectra",
    )
    args: argparse.Namespace = parser.parse_args()
    return args


def Load_Directory(args) \
        -> tuple[tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]],
                 tuple[float, float], tuple[float, float]]:

    import censo_ext.anmr as anmr
    in_dir: tuple[Path, Path] = args.dir[0], args.dir[1]
    h_limits: tuple[float, float] = args.h_limits
    c_limits: tuple[float, float] = args.c_limits
    directory_H, directory_C = in_dir
    args_x: dict = {"auto": True, "average": args.average, "bobyqa": False, "mf": 500,
                    "dir": directory_H, "thr": None, "json": None, "thrab": 0.025,
                    "verbose": False, "lw": 1, "tb": 4, "mss": 10, "cutoff": 0.001,
                    "show": False, "start": None, "end": None, "out": "output.npz"}
    sys.stdout = open(os.devnull, 'w')
    data_x: npt.NDArray[np.float64] = anmr.main(argparse.Namespace(**args_x)).T
    sys.stdout = sys.__stdout__

    args_y: dict = {"auto": True, "average": args.average, "bobyqa":  False, "mf": 500,
                    "dir": directory_C, "thr": None, "json": None, "thrab": 0.025,
                    "verbose": False, "lw": 1, "tb": 4, "mss": 10, "cutoff": 0.001,
                    "show": False, "start": None, "end": None, "out": "output.npz"}
    sys.stdout = open(os.devnull, 'w')
    data_y: npt.NDArray[np.float64] = anmr.main(argparse.Namespace(**args_y)).T
    sys.stdout = sys.__stdout__

    if h_limits is None:
        h_limits = float(min(data_x.T[0])), float(max(data_x.T[0]))
    else:
        start, end = h_limits
        if start > end:
            h_limits = end, start
    if c_limits is None:
        c_limits = float(min(data_y.T[0])), float(max(data_y.T[0]))
    else:
        start, end = c_limits
        if start > end:
            c_limits = end, start

    return (data_x, data_y), h_limits, c_limits


def draw_2D_basic(data_xy: tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]) -> tuple[Axes, tuple[Axes, float]]:
    data_x, data_y = data_xy
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
    y_axis_data: npt.NDArray[np.float64] = data_y.T[1]

    fig.suptitle("$^{1}$J (C,H)", fontsize=12, x=0.10, y=0.98)
    ax_histx.plot(data_x.T[0], x_axis_data)
    ax_histy.plot(-y_axis_data, data_y.T[0])
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')

    y_lowest = (-data_y).min()
    return ax, (ax_histy, y_lowest)


def plot_2D_slice(ax: Axes, ax_histy_float, in_dir: tuple[Path, Path], h_limits: tuple[float, float], c_limits: tuple[float, float]) \
        -> tuple[dict[AtomID, int], dict[AtomID, npt.NDArray], dict[AtomID, float], dict[AtomID, float]]:

    ax_histy, y_lowest = ax_histy_float
    from censo_ext.Tools.ml4nmr import read_mol_neighbors_bond_order
    from ase.atoms import Atoms
    Directory_H, Directory_C = in_dir
    mol: Atoms | list[Atoms]
    neighbor: dict[AtomID, npt.NDArray[np.int64]]
    bond_order: dict[AtomID, int]
    mol, neighbor, bond_order = read_mol_neighbors_bond_order(
        Directory_H/Path("crest_conformers.xyz"))
    idx_H_atom: list[int] = [idx+1 for idx,
                             i in enumerate(mol) if i.symbol == "H"]  # type: ignore # nopep8
    idx_C_atom: list[int] = [idx+1 for idx,
                             i in enumerate(mol) if i.symbol == "C"]  # type: ignore # nopep8

    neighbor: dict[AtomID, npt.NDArray[np.int64]] = {key: value for key,
                                                     value in neighbor.items() if key in idx_C_atom}
    for key, value in neighbor.items():
        neighbor[key] = np.array([x for x in value if x in idx_H_atom])

    tmp_c: list = list(np.genfromtxt(
        Directory_C / Path("Average/NMR/orcaS.out"), usecols=[0, 1]))

    Atoms_C: dict[AtomID, float] = {AtomID(int(x)): y for x, y in tmp_c}

    tmp_h: list = list(np.genfromtxt(
        Directory_H/Path("Average/NMR/orcaS.out"), usecols=[0, 1]))
    Atoms_H: dict[AtomID, float] = {AtomID(int(x)): y for x, y in tmp_h}
    ax.set_xlim(h_limits[1], h_limits[0])
    ax.set_ylim(c_limits[1], c_limits[0])
    x_lowest, x_highest = h_limits
    x_lowest = abs(x_lowest-x_highest)*0.03 + x_lowest

    for idx_C, C_ppm in Atoms_C.items():

        idx0_neighbor: dict = {}
        for idx_neighbor_Atoms_H in neighbor[idx_C]:
            for idx, value in enumerate(Atoms_H.keys()):
                if idx_neighbor_Atoms_H == value:
                    if value in idx0_neighbor:
                        idx0_neighbor[idx] = (idx0_neighbor[idx], value)
                    else:
                        idx0_neighbor[idx] = value

        if len(idx0_neighbor) != 0:
            for idx0, value in idx0_neighbor.items():

                import censo_ext.anmr as anmr
                x: dict = {'out': 'output.npz', 'mf': 500.0, "dir": Directory_H, 'lw': None,
                           'thr': None, 'thrab': 0.025, "verbose": False, 'tb': 4,
                           'cutoff': 0.001, 'start': None, 'end': None, 'show': False, 'mss': 10, 'auto': True,
                           'average': True, 'bobyqa': False, 'json': [idx0]}
                sys.stdout = open(os.devnull, 'w')
                np_dat: npt.NDArray[np.float64] = anmr.main(
                    args=argparse.Namespace(**x))
                sys.stdout = sys.__stdout__
                maximum: np.float64 = np.max(np_dat)
                ax.plot(np_dat[0], -np_dat[1]/maximum*5 + C_ppm, linewidth=1)
                ax_histy.text(y_lowest, C_ppm, f"{C_ppm:12.2f}",
                              ha="right", va="center", fontsize=6)
                if len(idx0_neighbor.values()) == 1:
                    text: Text = ax.text(
                        x_lowest, C_ppm, f"({idx_C}C,", ha="right", va="center", fontsize=8)
                    text = ax.annotate(f" {tuple(idx0_neighbor.values())[0]}H)",
                                       xycoords=text, xy=(1.00, 0.5), ha="left", va="center", color="blue", fontsize=8)
                else:
                    text = ax.text(
                        x_lowest, C_ppm, f"({idx_C}C,", ha="right", va="center", fontsize=8)
                    text = ax.annotate(f" {tuple(idx0_neighbor.values())}H)",
                                       xycoords=text, xy=(1.00, 0.5), ha="left", va="center", color="blue", fontsize=8)

    plt.subplots_adjust(hspace=0.5, wspace=0.5)
    plt.show()
    return bond_order, neighbor, Atoms_H, Atoms_C


def print_report(bond_order: dict[AtomID, int], neighbor: dict[AtomID, npt.NDArray], Atoms_H: dict[AtomID, float], Atoms_C: dict[AtomID, float]) -> None:

    print("   #C   Bond_Order   13C(HSQC)      1H(HSQC)        #H ")
    for idxAtom_C, C_ppm in Atoms_C.items():
        if bond_order[idxAtom_C] == 0:
            print(f"{idxAtom_C:>5d}     C   {C_ppm:>15.4f}", end="")
        else:
            print(f"{idxAtom_C:>5d}     CH{bond_order[idxAtom_C]:>1d} {C_ppm:>15.4f}", end="")  # nopep8

        idx1_neighbor: list[AtomID] = []
        for idx1_neighbor_Atoms_H in neighbor[idxAtom_C]:
            for idx, value in enumerate(Atoms_H.keys()):
                if idx1_neighbor_Atoms_H == value:
                    idx1_neighbor.append(value)

        if len(idx1_neighbor) == 0:
            print("")
        elif len(idx1_neighbor) == 1:
            for _, x in enumerate(idx1_neighbor):
                print(f"{Atoms_H[x]:>15.4f} {int(x):>10d}", end="")
                print("")
        else:
            for idx, x in enumerate(idx1_neighbor):
                if idx == 0:
                    print(f"{Atoms_H[x]:>15.4f} {int(x):>10d}", end="")
                else:
                    print("\n", " "*27,
                          f"{Atoms_H[x]:>15.4f} {int(x):>10d}", end="")
            print("")
    return


def main(args: argparse.Namespace = argparse.Namespace()) -> None:
    if args == argparse.Namespace():
        args = cml()
    print_arguments()

    in_dir: tuple[Path, Path] = args.dir[0], args.dir[1]

    data_xy, h_limits, c_limits = Load_Directory(args)
    ax, ax_histy_float = draw_2D_basic(data_xy)
    bond_order, neighbor, idxAtoms_H, idxAtoms_C = plot_2D_slice(
        ax, ax_histy_float, in_dir, h_limits, c_limits)
    print_report(bond_order, neighbor, idxAtoms_H, idxAtoms_C)


if __name__ == "__main__":
    main()
