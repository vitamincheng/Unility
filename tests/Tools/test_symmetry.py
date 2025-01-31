from pathlib import Path
import pytest
from pointgroup import PointGroup


def test_get_Point_Group():
    # for crest_confomrers.xyz
    from censo_ext.Tools.xyzfile import GeometryXYZs
    xyz = GeometryXYZs(Path("tests/data/crest_conformers.xyz"))
    xyz.method_read_xyz()
    for idx in range(len(xyz)):
        pos = ([a.tolist() for a in xyz.Sts[idx].coord])
        sym = ([a for a in xyz.Sts[idx].names.values()])
        pg1 = PointGroup(pos, sym)
        Result = pg1.get_point_group()
        assert Result == 'C1'

    # for isomers.xyz
    xyz = GeometryXYZs(Path("tests/data/isomers.xyz"))
    xyz.method_read_xyz()

    expected_Result = ["Cs", "C2", "C1"]

    for idx in range(len(xyz)):
        pos = ([a.tolist() for a in xyz.Sts[idx].coord])
        sym = ([a for a in xyz.Sts[idx].names.values()])
        pg1 = PointGroup(pos, sym)
        Result = pg1.get_point_group()
        assert Result == expected_Result[idx]
