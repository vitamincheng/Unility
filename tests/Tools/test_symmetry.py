import pathlib
import argparse
from pathlib import Path
import pytest
from pointgroup import PointGroup
from icecream import ic
import Tools.symmetry


def test_get_Point_Group():

    from Tools.xyzfile import ClassGeometryXYZs
    xyz = ClassGeometryXYZs(Path("tests/data/crest_conformers.xyz"))
    xyz.method_read_xyz()
    for idx in range(len(xyz)):
        pos = ([a.tolist() for a in xyz.Sts[idx].coord])
        sym = ([a for a in xyz.Sts[idx].names.values()])
        pg1 = PointGroup(pos, sym)
        Result = pg1.get_point_group()
        assert Result == 'C1'
