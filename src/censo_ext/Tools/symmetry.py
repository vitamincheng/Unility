#!/usr/bin/env python

from pointgroup import PointGroup
import numpy as np
import numpy.typing as npt
from censo_ext.Tools.xyzfile import Geometry


def method_get_point_group(Sts: list[Geometry], idx: int, Hydrogen: bool) -> str:
    """
    Determine the point group symmetry for a specified molecular geometry.

    This function calculates the point group of a molecule by analyzing the 
    atomic coordinates and symbols, with optional exclusion of hydrogen atoms.

    Args:
        Sts(list[Geometry]): List of molecular geometry objects containing 
                             atomic coordinates and element names
        idx(int): Index specifying which geometry in the list to analyze
        Hydrogen(bool): Flag indicating whether to include hydrogen atoms 
                         in the symmetry calculation (True = include, False = exclude)

    Returns:
        str: Point group symbol (e.g., 'C2v', 'D3h', 'Oh') representing the 
             molecular symmetry

    Example:
        >>> Sts = [geometry1, geometry2, ...]
        >>> point_group = method_get_point_group(Sts, 0, False)
        >>> print(f"Molecular point group: {point_group}")
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
