#!/usr/bin/env python3

from pointgroup import PointGroup
from icecream import ic
import numpy as np


def method_get_point_group(st: list, idx: int, Hydrogen: bool) -> str:

    if Hydrogen == False:
        idx_list = [key-1 for key, value in st[idx].names.items()
                    if value != 'H']

        # ic(idx_list)
        # ic(np.array(list(st[idx].names.values()))[idx_list])
        # ic(np.array(st[idx].coordinates)[idx_list])
        pos = np.array(st[idx].coordinates)[idx_list]
        sym = np.array(list(st[idx].names.values()))[idx_list]

    else:
        pos = np.array([a.tolist() for a in st[idx].coordinates])
        sym = np.array([a for a in st[idx].names.values()])

    pg1 = PointGroup(pos, sym)
    return pg1.get_point_group()
