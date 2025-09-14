from pathlib import Path
import pytest

from censo_ext.Tools.symmetry import method_get_point_group


def test_symmetry_crest_conformers_miss_args():

    from censo_ext.Tools.xyzfile import GeometryXYZs
    FileName = Path("tests/data/crest_conformers000.xyz")
    xyz = GeometryXYZs("tests/data/crest_conformers000.xyz")

    with pytest.raises(FileNotFoundError) as e:
        xyz.method_read_xyz()
    assert str(e.value) == f"  The file {FileName} is not Exist ..."


def test_symmetry_crest_conformers():

    from censo_ext.Tools.xyzfile import GeometryXYZs
    xyzFile = GeometryXYZs("tests/data/crest_conformers.xyz")
    xyzFile.method_read_xyz()
    for idx in range(len(xyzFile)):
        assert method_get_point_group(xyzFile.Sts, idx, False) == 'C1'


def test_symmetry_crest_conformers_H():

    from censo_ext.Tools.xyzfile import GeometryXYZs
    xyzFile = GeometryXYZs("tests/data/crest_conformers.xyz")
    xyzFile.method_read_xyz()
    for idx in range(len(xyzFile)):
        assert method_get_point_group(xyzFile.Sts, idx, True) == 'C1'


def test_symmetry_isomers():
    from censo_ext.Tools.xyzfile import GeometryXYZs
    # for isomers.xyz
    xyzFile = GeometryXYZs("tests/data/isomers.xyz")
    xyzFile.method_read_xyz()

    expected_Res: list[str] = ["Cs", "C2", "C1"]

    for idx in range(len(xyzFile)):
        assert method_get_point_group(
            xyzFile.Sts, idx, False) == expected_Res[idx]


def test_symmetry_isomers_H():
    from censo_ext.Tools.xyzfile import GeometryXYZs
    # for isomers.xyz
    xyzFile = GeometryXYZs("tests/data/isomers.xyz")
    xyzFile.method_read_xyz()

    expected_Res: list[str] = ["Cs", "C2", "C1"]

    for idx in range(len(xyzFile)):
        assert method_get_point_group(
            xyzFile.Sts, idx, True) == expected_Res[idx]
