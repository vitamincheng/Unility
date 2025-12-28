#!/usr/bin/env python
# from censo_ext.Tools.ml4nmr import read_mol_neighbors
from censo_ext import xyzReturnOandZ
from censo_ext.Tools.Parameter import Eh
from censo_ext.Tools.spectra import Boltzmann_Weighting
# from censo_ext.Tools.symmetry import method_get_point_group
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
        "-b",
        "--boltz",
        dest="boltz",
        action="store_true",
        required=False,
        help="only use Boltzmann distribution, not use interia to sort [default False]",
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
        "-bthr",
        "--bthr",
        dest="bthr",
        action="store",
        required=False,
        type=float,
        default=100.0,
        help="the threshold of interia [default 100.0 (amu/A^2)]",
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

    xyzFile: GeometryXYZs = GeometryXYZs(args.file)
    xyzFile.method_read_xyz()
    nSts_origin: int = len(xyzFile)

    moment: dict = {"file": args.file, "auto": True,
                    "atom": None, "print": False, "replace": True, "out": None}
    x_args = argparse.Namespace(**moment)
    xyzReturnOandZ.main(x_args)

    xyzFile.method_compute_COM()
    xyzFile.method_compute_Inertia()
    import numpy as np
    import numpy.typing as npt
    list_inertia: list | npt.NDArray = []
    for St in xyzFile.Sts:
        list_inertia.append(St.inertia)
    _inertia = np.array(list_inertia)

    idx0_St_remove: list[list[int]] = []
    print("")
    print(" idx0_p     idx0_q")
    for idx0_p, moment in enumerate(_inertia):
        std = _inertia[idx0_p].copy()

        #
        for idx0_q, moment in enumerate(_inertia):
            _inertia[idx0_q] = moment - std

        idx0_diff: list[int] = []
        for idx0_q, moment in enumerate(_inertia):

            if idx0_p != idx0_q and np.sum(np.square(np.array(moment))) <= args.bthr and not args.boltz:
                idx0_diff.append(idx0_q)
                print(f"{idx0_p+1:5d}    | {idx0_q+1:5d}")

        if len(idx0_diff) != 0:
            idx0_diff = [x for x in idx0_diff if x > idx0_p]
            if len(idx0_diff) >= 0:
                idx0_St_remove.append(idx0_diff)
        # print(idx0_St_remove)
        # print("")

    import itertools
    idx0_St_remove_flat: list[int] = list(set(
        itertools.chain.from_iterable(idx0_St_remove)))
    idx0_index: list[int] = [* range(len(_inertia))]
    idx0_index = [x for x in idx0_index if x not in idx0_St_remove_flat]

    xyzFile.method_xyzExtract(idx0_index)
    Energy: list[float] = [St._comment_energy for St in np.array(xyzFile.Sts)]
    nClusters: list[int] = [
        St.comment_nClusters for St in np.array(xyzFile.Sts)]

    import numpy as np
    import numpy.typing as npt
    np_Energy = np.array(Energy)*Eh
    np_Energy = np_Energy - np_Energy.min()
    intp_Energy: npt.NDArray[np.intp] = np.argsort(np_Energy)
    intp_remove_Energy = np.argwhere(np_Energy[intp_Energy] > args.ewin)

    intp_Energy = np.delete(intp_Energy, intp_remove_Energy)
    # print(np.array(nClusters)[intp_Energy])

    BW: npt.NDArray[np.float64] = Boltzmann_Weighting(
        np_Energy[intp_Energy], TEMP=args.temp)

    zip_energy: zip[tuple[npt.NDArray[np.intp], npt.NDArray[np.float64], npt.NDArray[np.float64]]] = zip(
        intp_Energy, np.array(nClusters)[intp_Energy], np.array(np_Energy[intp_Energy]), BW)

    # print the parameter of the Boltzmann weighting
    print("")
    print("  ===== Boltzmann Distribution =====")
    print(f"  threshold energy              = {args.ewin} (kcal/mol)")
    if args.boltz:
        print("  Using only Boltzmann Distribution (not sort by inertia)")
    else:
        print(f"  threshold of inertia          = {args.bthr} (amu/A^2)")

    print(f"  Temperature                   = {args.temp} (K)")
    print(f"  The numbers of Start Cluster  = {nSts_origin} ")
    print(f"  The numbers of Final Cluster  = {len(intp_Energy)} ")
    if not args.boltz:
        print(
            f"  The indexes of remove         = {[x+1 for x in idx0_St_remove_flat]}")
    print(f"  Saved file                    = {args.out} ")

    print("")
    print("  ===== Boltzmann Weighting Table =====")
    print("  index1           Energy (kcal/mol)             BW")
    for intp, moment, y, z in zip_energy:

        print(f"{intp+1:5d} {moment:8d}           {y:17.10f}       {z:8.4f}")
    print("  ===== Finished =====")
    print("")

    xyzFile.set_filename(args.out)
    xyzFile.method_save_xyz((intp_Energy+1).tolist())


if __name__ == "__main__":
    main()
