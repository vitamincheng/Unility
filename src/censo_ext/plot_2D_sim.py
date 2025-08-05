#! /usr/bin/env python3
import argparse
from icecream import ic
# import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
import os
import sys


def Load_dat(fileName_H, FileName_C) -> tuple[npt.NDArray, npt.NDArray]:

    data_x: npt.NDArray[np.float64] = np.genfromtxt(fileName_H)
    data_y: npt.NDArray[np.float64] = np.genfromtxt(FileName_C)
    return data_x, data_y


def Load_Directory(directory_H, directory_C) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:

    import censo_ext.anmr as anmr
    args_x: dict = {"auto": True, "average": True,
                    "bobyqa": True, "mf": 500, "dir": directory_H, "thr": None, "json": [-1], "thrab": 0.025,
                    "tb": 4, "mss": 9, "cutoff": 0.001, "show": False, "start": None, "end": None, "out": "output.dat"}
    data_x = anmr.main(argparse.Namespace(**args_x)).T
    args_y: dict = {"auto": True, "average": True,
                    "bobyqa": True, "mf": 500, "dir": directory_C, "thr": None, "json": [-1], "thrab": 0.025,
                    "tb": 4, "mss": 9, "cutoff": 0.001, "show": False, "start": None, "end": None, "out": "output.dat"}
    data_y = anmr.main(argparse.Namespace(**args_y)).T

    return data_x, data_y


def plot_2D_basic(data_x, data_y):
    fig = plt.figure(figsize=(11.7, 8.3), dpi=100)
    gs = fig.add_gridspec(2, 2,  width_ratios=(1, 19), height_ratios=(1, 9),
                          left=0.03, right=0.97, bottom=0.03, top=0.97,
                          wspace=0.1, hspace=0.1)

    ax = fig.add_subplot(gs[1, 1])
    ax_histx = fig.add_subplot(gs[0, 1], sharex=ax)
    ax_histy = fig.add_subplot(gs[1, 0], sharey=ax)
    ax_histx.get_xaxis().set_visible(False)
    ax_histx.get_yaxis().set_visible(False)
    ax_histx.axis('off')
    ax_histy.get_xaxis().set_visible(False)
    ax_histy.get_yaxis().set_visible(False)
    ax_histy.axis('off')

    x_axis_data = data_x.T[1]
    y_axis_data = data_y.T[1]

    # height_x_axis_data =np.max(x_axis_data)
    # height_y_axis_data =np.max(y_axis_data)
    ax_histx.plot(data_x.T[0], x_axis_data)
    ax_histy.plot(-y_axis_data, data_y.T[0])
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')
    return ax


def plot_2D_slice(ax, data_x, data_y) -> tuple[dict[int, int], dict[int, npt.NDArray], dict, dict]:

    from censo_ext.Tools.ml4nmr import read_mol_neighbors_bond_order
    from pathlib import Path
    from ase.atoms import Atoms
    mol: Atoms | list[Atoms]
    neighbor: dict[int, npt.NDArray[np.int64]]
    bond_order: dict[int, int]
    mol, neighbor, bond_order = read_mol_neighbors_bond_order(
        Path("Test/34.Ergocalciferol/04.Hydrogen/crest_conformers.xyz"))
    idx_H_atom: list[int] = [idx+1 for idx,
                             i in enumerate(mol) if i.symbol == "H"]  # type: ignore # nopep8
    idx_C_atom: list[int] = [idx+1 for idx,
                             i in enumerate(mol) if i.symbol == "C"]  # type: ignore # nopep8

    neighbor = {key: value for key,
                value in neighbor.items() if key in idx_C_atom}
    for key, value in neighbor.items():
        neighbor[key] = np.array([x for x in value if x in idx_H_atom])

    tmp = list(np.genfromtxt(
        "Test/34.Ergocalciferol/07.Carbon/Average/NMR/orcaS-BOBYQA.out", usecols=[0, 1]))

    idxAtoms_C: dict = {int(x): -y for x, y in tmp}
    # list_idxAtoms_C:list = list(idxAtoms_C.keys())

    tmp = list(np.genfromtxt(
        "Test/34.Ergocalciferol/04.Hydrogen/Average/NMR/orcaS-BOBYQA.out", usecols=[0, 1]))
    idxAtoms_H: dict = {int(x): -y for x, y in tmp}
    # list_idxAtoms_H:list = list(idxAtoms_H.keys())
    ax.set_xlim(max(data_x.T[0]),  min(data_x.T[0]))
    ax.set_ylim(max(data_y.T[0]), min(data_y.T[0]))

    for idxAtom_C, C_ppm in idxAtoms_C.items():
        idx0_neighbor: list = []
        for idx_neighbor_Atoms_H in neighbor[idxAtom_C]:
            for idx, value in enumerate(idxAtoms_H.keys()):
                if idx_neighbor_Atoms_H == value:
                    idx0_neighbor.append(idx)

        if len(idx0_neighbor) != 0:
            for idx0 in idx0_neighbor:

                import censo_ext.anmr as anmr
                x = {'out': 'output.dat', 'mf': 500.0, "dir": "Test/34.Ergocalciferol/04.Hydrogen", 'lw': None, 'ascal': None, 'bscal': None, 'thr': None, 'thrab': 0.025,
                     'tb': 4, 'cutoff': 0.001, 'start': None, 'end': None, 'show': False, 'mss': 9, 'auto': True, 'average': True, 'bobyqa': True, 'json': [idx0]}
                sys.stdout = open(os.devnull, 'w')
                np_dat = anmr.main(args=argparse.Namespace(**x))
                sys.stdout = sys.__stdout__
                maximum = np.max(np_dat)
                ax.plot(np_dat[0], -np_dat[1]/maximum*5 + C_ppm, linewidth=0.5)
                ax.text(0, C_ppm, f"{C_ppm:12.3f}",
                        ha="right", va="center", fontsize=6)

    plt.subplots_adjust(hspace=0.5, wspace=0.5)
    return bond_order, neighbor, idxAtoms_H, idxAtoms_C


def print_report(bond_order, neighbor, idxAtoms_H, idxAtoms_C) -> None:

    print("   #C   Bond_Order   13C(HSQC)      1H(HSQC)        #H ")

    for idxAtom_C, C_ppm in idxAtoms_C.items():

        if bond_order[idxAtom_C] == 0:
            print(f"{idxAtom_C:>5d}     C   {C_ppm:>15.4f}", end="")
        else:
            print(f"{idxAtom_C:>5d}     CH{bond_order[idxAtom_C]:>1d} {C_ppm:>15.4f}", end="")  # nopep8

        idx0_neighbor: list = []
        for idx_neighbor_Atoms_H in neighbor[idxAtom_C]:
            for idx, value in enumerate(idxAtoms_H.keys()):
                if idx_neighbor_Atoms_H == value:
                    idx0_neighbor.append(value)

        if len(idx0_neighbor) == 0:
            print("")
        elif len(idx0_neighbor) == 1:
            for idx, x in enumerate(idx0_neighbor):
                print(f"{idxAtoms_H[x]:>15.4f} {int(x):>10d}", end="")
                print("")
        else:
            for idx, x in enumerate(idx0_neighbor):
                if idx == 0:
                    print(f"{idxAtoms_H[x]:>15.4f} {int(x):>10d}", end="")
                else:
                    print("\n", " "*27,
                          f"{idxAtoms_H[x]:>15.4f} {int(x):>10d}", end="")
            print("")
    return


def plot_target(ax) -> None:
    data: npt.NDArray[np.float64] = np.genfromtxt("out")
    ic(data)
    for i in data:
        ax.scatter(i[1], i[0], marker="o", color="r",
                   s=30, alpha=0.5)  # type: ignore
    return


def main(args=argparse.Namespace()) -> None:
    if args == argparse.Namespace():
        pass
    # hidden = True

    directory_H: str = "Test/34.Ergocalciferol/04.Hydrogen"
    directory_C: str = "Test/34.Ergocalciferol/07.Carbon"
    data_x, data_y = Load_Directory(directory_H, directory_C)

    # data_x, data_y = Load_dat(
    #    directory_H + "/output.dat", directory_C+"/output.dat")

    ax = plot_2D_basic(data_x, data_y)
    bond_order, neighbor, idxAtoms_H, idxAtoms_C = plot_2D_slice(
        ax, data_x, data_y)
    # plot_target(ax)
    print_report(bond_order, neighbor, idxAtoms_H, idxAtoms_C)
    from censo_ext.Tools.utility import save_figure
    save_figure()
    # if hidden == False:
    ax.set_xlim(max(data_x.T[0]), min(data_x.T[0]))
    ax.set_ylim(max(data_y.T[0]), min(data_y.T[0]))
    plt.show()


if __name__ == "__main__":
    main()
