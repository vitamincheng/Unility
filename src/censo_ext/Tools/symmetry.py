#!/usr/bin/env python3

from pointgroup import PointGroup
from icecream import ic
import numpy as np
import numpy.typing as npt
from censo_ext.Tools.xyzfile import Geometry


def method_get_point_group(Sts: list[Geometry], idx: int, Hydrogen: bool) -> str:
    # Sts : the list of xyz files
    # idx : the number of xyz file
    pos: npt.NDArray[np.float64]
    sym: npt.NDArray[np.float64]

    if Hydrogen == False:
        idx0_names: list[int] = [key-1 for key, value in Sts[idx].names.items()
                                 if value != 'H']

        pos = np.array(Sts[idx].coord)[idx0_names]
        sym = np.array(list(Sts[idx].names.values()))[idx0_names]

    else:
        pos = np.array([a.tolist() for a in Sts[idx].coord])
        sym = np.array([a for a in Sts[idx].names.values()])

    return PointGroup(pos, sym).get_point_group()
