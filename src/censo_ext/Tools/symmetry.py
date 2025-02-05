#!/usr/bin/env python3

from pointgroup import PointGroup
from icecream import ic
import numpy as np


def method_get_point_group(Sts: list, idx: int, Hydrogen: bool) -> str:
    # Sts : the list of xyz files
    # idx : the number of xyz file

    if Hydrogen == False:
        idx_names: list = [key-1 for key, value in Sts[idx].names.items()
                           if value != 'H']

        pos: np.ndarray = np.array(Sts[idx].coordinate)[idx_names]
        sym: np.ndarray = np.array(list(Sts[idx].names.values()))[idx_names]

    else:
        pos: np.ndarray = np.array([a.tolist() for a in Sts[idx].coordinates])
        sym: np.ndarray = np.array([a for a in Sts[idx].names.values()])

    pg1 = PointGroup(pos, sym)
    return pg1.get_point_group()
