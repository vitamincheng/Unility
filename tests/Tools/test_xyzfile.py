from pathlib import Path
import pytest
from censo_ext.Tools.xyzfile import GeometryXYZs
import filecmp


def test_xyzfile_read_xyz_miss_args():
    fileName: Path = Path("tests/test.xyz")
    xyzFile = GeometryXYZs(fileName)
    with pytest.raises(FileNotFoundError) as e:
        xyzFile.method_read_xyz()
    assert str(e.value) == f"  The file {fileName} is not Exist ..."


@pytest.mark.parametrize(argnames="input_Path,len_Sts_file,energy1,energy2",
                         argvalues=[(Path("tests/data/crest_conformers.xyz"), 54,
                                     -1168.36168282, -1168.35637181),
                                    (Path("tests/data/isomers.xyz"), 3,
                                     -309.820981116391, -309.815154391997)])
def test_xyzfile_read_xyz(input_Path: Path, len_Sts_file: int, energy1: float, energy2: float):
    xyzFile: GeometryXYZs = GeometryXYZs(input_Path)
    assert len(xyzFile) == 0
    xyzFile.method_read_xyz()
    assert len(xyzFile) == len_Sts_file
    assert xyzFile.Sts[0].get_comment_energy() == energy1
    assert xyzFile.Sts[-1].get_comment_energy() == energy2
    assert len(xyzFile.get_comment_energy()) == len_Sts_file
    assert xyzFile.get_comment_energy()[0] == energy1
    assert xyzFile.get_comment_energy()[-1] == energy2


@pytest.mark.parametrize(argnames="input_Path,compare_str1,compare_str2",
                         argvalues=[(Path("tests/data/crest_conformers.xyz"),
                                     " Energy =        -1168.3616828200 Eh        #Cluster:     0",
                                     " Energy =        -1168.3563718100 Eh        #Cluster:     0"),
                                    (Path("tests/data/isomers.xyz"),
                                     " Energy =        -309.8209811164 Eh        #Cluster:     1",
                                     " Energy =        -309.8151543920 Eh        #Cluster:     3")])
def test_xyzfile_rewrite_comment(input_Path: Path, compare_str1: str, compare_str2: str):
    xyzFile: GeometryXYZs = GeometryXYZs(input_Path)
    xyzFile.method_read_xyz()
    xyzFile.method_rewrite_comment()
    assert xyzFile.Sts[0].comment == compare_str1
    assert xyzFile.Sts[-1].comment == compare_str2


@pytest.mark.parametrize(argnames="input_Path,compare_str1,compare_str2",
                         argvalues=[(Path("tests/data/crest_conformers.xyz"),
                                     " Energy =        -1168.3616828200 Eh        #Cluster:     1",
                                     " Energy =        -1168.3563718100 Eh        #Cluster:     54"),
                                    (Path("tests/data/isomers.xyz"),
                                     " Energy =        -309.8209811164 Eh        #Cluster:     1",
                                     " Energy =        -309.8151543920 Eh        #Cluster:     3")])
def test_xyzfile_comment_new(input_Path: Path, compare_str1: str, compare_str2: str):
    xyzFile: GeometryXYZs = GeometryXYZs(input_Path)
    xyzFile.method_read_xyz()
    xyzFile.method_comment_new()
    assert xyzFile.Sts[0].comment == compare_str1
    assert xyzFile.Sts[-1].comment == compare_str2


@pytest.mark.parametrize(argnames="input_Path,compare_str1,compare_str2",
                         argvalues=[(Path("tests/data/crest_conformers.xyz"),
                                     " Energy =        -1168.3616828200 Eh        #Cluster:     0",
                                     " Energy =        -1168.3563718100 Eh        #Cluster:     0"),
                                    (Path("tests/data/isomers.xyz"),
                                     " Energy =        -309.8209811164 Eh        #Cluster:     1",
                                     " Energy =        -309.8151543920 Eh        #Cluster:     3")])
def test_xyzfile_comment_keep(input_Path: Path, compare_str1: str, compare_str2: str):
    xyzFile: GeometryXYZs = GeometryXYZs(input_Path)
    xyzFile.method_read_xyz()
    xyzFile.method_comment_keep()
    assert xyzFile.Sts[0].comment == compare_str1
    assert xyzFile.Sts[-1].comment == compare_str2


@pytest.mark.parametrize(argnames="input_bool,input_Path,compare_filename",
                         argvalues=[(True, Path("tests/data/crest_conformers.xyz"), "tests/compare/xyzfile-1.xyz"),
                                    (True, Path("tests/data/isomers.xyz"),
                                     "tests/compare/xyzfile-2.xyz"),
                                    (False, Path("tests/data/isomers.xyz"), "tests/compare/xyzfile-1.xyz")])
def test_xyzfile_save_xyz(input_bool: bool, input_Path: Path, compare_filename: str):
    xyzFile: GeometryXYZs = GeometryXYZs(input_Path)
    xyzFile.method_read_xyz()
    xyzFile.method_comment_new()
    outFile = Path("tests/compare/output.xyz")
    xyzFile.set_filename(outFile)
    xyzFile.method_save_xyz([])
    assert filecmp.cmp(outFile, compare_filename) == input_bool
    outFile.unlink()
