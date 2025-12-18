#!/usr/bin/env python
# from censo_ext.Tools.ml4nmr import read_mol_neighbors
from censo_ext import xyzReturnOandZ
from censo_ext.Tools.Parameter import Eh
from censo_ext.Tools.spectra import Boltzmann_Weighting
from censo_ext.Tools.utility import print_arguments
import argparse

from censo_ext.Tools.xyzfile import GeometryXYZs
descr = """
    ________________________________________________________________________________
    | For Generation of xyz molecule
    | Usages   : xyz.py <geometry> [options]
    | [options]
    |______________________________________________________________________________
    """


def cml() -> argparse.Namespace:
    """ Get args object from commandline interface. Needs argparse module."""
    parser = argparse.ArgumentParser(
        description="descr",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage=argparse.SUPPRESS)
    parser.add_argument(
        "-i",
        "--input",
        dest="file",
        action="store",
        required=False,
        default="isomers.xyz",
        help="Provide one input xyz file [default isomers.xyz]",
    )
    parser.add_argument(
        "-o",
        "--out",
        dest="out",
        action="store",
        required=False,
        default="clusters.xyz",
        help="Provide one input xyz file [default clusters.xyz]",
    )
    parser.add_argument(
        "--ewin",
        dest="ewin",
        action="store",
        required=False,
        type=float,
        default=100,
        help="the electron energy threshold [default 100 (Kcal/mol)]",
    )
    parser.add_argument(
        "-rthr",
        "--rthr",
        dest="rthr",
        action="store",
        required=False,
        type=float,
        default=1.0,
        help="the threshold of interia [default 1.0 (amu/A^2)]",
    )
    parser.add_argument(
        "--temp",
        dest="temp",
        action="store",
        required=False,
        type=float,
        default=298.15,
        help="the temperature [default 298.15 K]",
    )
    args: argparse.Namespace = parser.parse_args()
    return args


def main(args: argparse.Namespace = argparse.Namespace()) -> None:
    if args == argparse.Namespace():
        args = cml()
    print_arguments()

    x: dict = {"file": args.file, "auto": True,
               "atom": None, "print": False, "replace": True, "out": None}
    x_args = argparse.Namespace(**x)
    xyzReturnOandZ.main(x_args)

    xyzFile: GeometryXYZs = GeometryXYZs(args.file)
    xyzFile.method_read_xyz()
    nSts_origin = len(xyzFile.Sts)
    xyzFile.method_compute_COM()
    xyzFile.method_compute_Inertia()
    import numpy as np
    import numpy.typing as npt
    list_inertia: list | npt.NDArray = []
    for St in xyzFile.Sts:
        list_inertia.append(St.inertia)
    _inertia = np.array(list_inertia)

    idx0_St_remove: list = []
    for index0, x in enumerate(_inertia):
        std = _inertia[index0].copy()

        for idx0, x in enumerate(_inertia):
            _inertia[idx0] = x - std

        idx0_diff: list[int] = []
        for idx0, x in enumerate(_inertia):

            # print(index0, idx0)
            The_Same_St = xyzFile.method_compare_the_same_core(
                index0, idx0)
            if np.sum(np.square(np.array(x))) <= args.rthr and The_Same_St:
                idx0_diff.append(idx0)
                # print(index0, idx0, end="")
                # print(" ===")
        # print(idx0_diff)
        if len(idx0_diff) != 1:
            idx0_diff = [x for x in idx0_diff if x > index0]
            if len(idx0_diff) >= 1:
                idx0_St_remove.append(idx0_diff)
        # print(idx0_St_remove)
        # print("")

    import itertools
    idx0_St_remove = list(itertools.chain.from_iterable(idx0_St_remove))

    idx0_index: list[int] = [* range(len(_inertia))]
    idx0_index = [x for x in idx0_index if x not in idx0_St_remove]

    xyzFile.method_xyzExtract(idx0_index)
    Energy: list[float] = [St._comment_energy for St in np.array(xyzFile.Sts)]

    import numpy as np
    import numpy.typing as npt
    np_Energy = np.array(Energy)*Eh
    np_Energy = np_Energy - np_Energy.min()
    intp_Energy: npt.NDArray[np.intp] = np.argsort(np_Energy)
    intp_remove_Energy = np.argwhere(np_Energy[intp_Energy] > args.ewin)

    intp_Energy = np.delete(intp_Energy, intp_remove_Energy)

    BW: npt.NDArray[np.float64] = Boltzmann_Weighting(
        np_Energy[intp_Energy], TEMP=args.temp)

    zip_energy: zip[tuple[npt.NDArray[np.intp], npt.NDArray[np.float64], npt.NDArray[np.float64]]] = zip(
        intp_Energy+1, np.array(np_Energy[intp_Energy]), BW)

    # print the parameter of the Boltzmann weighting
    print("")
    print("  ===== Boltzmann Distribution =====")
    print(f"  threshold energy              = {args.ewin} (kcal/mol)")
    print(f"  threshold of inertia          = {args.rthr} (amu/A^2)")
    print(f"  Temperature                   = {args.temp} (K)")
    print(f"  The numbers of Start Clusters = {nSts_origin} ")
    print(f"  The numbers of Final Clusters = {len(intp_Energy)} ")
    print(f"  Saved File Name               = {args.out} ")

    print("")
    print("  ===== Boltzmann Weighting Table =====")
    print("  index1           Energy (kcal/mol)             BW")
    for x, y, z in zip_energy:
        print(f"{x:8d}           {y:17.10f}       {z:8.4f}")
    print("  ===== Finished =====")
    print("")

    xyzFile.set_filename(args.out)
    xyzFile.method_rewrite_comment()
    xyzFile.method_comment_new()
    xyzFile.method_save_xyz((intp_Energy+1).tolist())


if __name__ == "__main__":
    main()
