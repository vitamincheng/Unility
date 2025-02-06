from pathlib import Path
import pytest
from pointgroup import PointGroup

from censo_ext.Tools.symmetry import method_get_point_group


def test_method_test_point_group_crest_conformers():

    from censo_ext.Tools.xyzfile import GeometryXYZs
    xyz = GeometryXYZs(Path("tests/data/crest_conformers.xyz"))
    xyz.method_read_xyz()
    for idx in range(len(xyz)):
        Result = method_get_point_group(xyz.Sts, idx, False)
        assert Result == 'C1'


def test_method_test_point_group_crest_conformers_H():

    from censo_ext.Tools.xyzfile import GeometryXYZs
    xyz = GeometryXYZs(Path("tests/data/crest_conformers.xyz"))
    xyz.method_read_xyz()
    for idx in range(len(xyz)):
        Result = method_get_point_group(xyz.Sts, idx, True)
        assert Result == 'C1'


def test_method_test_point_group_isomers():
    from censo_ext.Tools.xyzfile import GeometryXYZs
    # for isomers.xyz
    xyz = GeometryXYZs(Path("tests/data/isomers.xyz"))
    xyz.method_read_xyz()

    expected_Result = ["Cs", "C2", "C1"]

    for idx in range(len(xyz)):
        Result = method_get_point_group(xyz.Sts, idx, False)
        assert Result == expected_Result[idx]


def test_method_test_point_group_isomers_H():
    from censo_ext.Tools.xyzfile import GeometryXYZs
    # for isomers.xyz
    xyz = GeometryXYZs(Path("tests/data/isomers.xyz"))
    xyz.method_read_xyz()

    expected_Result = ["Cs", "C2", "C1"]

    for idx in range(len(xyz)):
        Result = method_get_point_group(xyz.Sts, idx, True)
        assert Result == expected_Result[idx]
