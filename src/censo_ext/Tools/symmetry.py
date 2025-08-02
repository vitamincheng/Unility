#!/usr/bin/env python3

from pointgroup import PointGroup
import numpy as np
import numpy.typing as npt
from censo_ext.Tools.xyzfile import Geometry


def method_get_point_group(Sts: list[Geometry], idx: int, Hydrogen: bool) -> str:
    """
    Calculate the point group symmetry of a given xyz file.

    Args:
        Sts(list[Geometry]): A list of Geometry objects.
        idx(int): The index of the xyz file to process.
        Hydrogen(bool): Whether to include hydrogen atoms in the calculation.

    Returns:
        str: The point group symmetry as a string.
    """
    # Sts : the list of xyz files
    # idx : the number of xyz file
    pos: npt.NDArray[np.float64]
    sym: npt.NDArray[np.float64]

    if not Hydrogen:
        idx0_names: list[int] = [key-1 for key, value in Sts[idx].names.items()
                                 if value != 'H']

        pos = np.array(Sts[idx].coord)[idx0_names]
        sym = np.array(list(Sts[idx].names.values()))[idx0_names]

    else:
        pos = np.array([a.tolist() for a in Sts[idx].coord])
        sym = np.array([a for a in Sts[idx].names.values()])

    return PointGroup(pos, sym).get_point_group()
