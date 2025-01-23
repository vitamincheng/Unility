import pathlib
import argparse
from pathlib import Path
import pytest
from Tools.xyzfile import ClassGeometryXYZs


def test_method_read_xyz():
    file = ClassGeometryXYZs(Path("tests/data/crest_conformers.xyz"))
    assert len(file) == 0
    file.method_read_xyz()
    assert len(file) == 54
    assert file.Sts[0].get_comment_energy() == -1168.36168282
    assert file.Sts[53].get_comment_energy() == -1168.35637181


def test_method_rewrite_comment():
    file = ClassGeometryXYZs(Path("tests/data/crest_conformers.xyz"))
    file.method_read_xyz()
    file.method_rewrite_comment()
    assert file.Sts[0].comment == " Energy =        -1168.36168282 Eh        #Cluster:     0"
    assert file.Sts[53].comment == " Energy =        -1168.35637181 Eh        #Cluster:     0"


def test_method_comment_new():
    file = ClassGeometryXYZs(Path("tests/data/crest_conformers.xyz"))
    file.method_read_xyz()
    file.method_comment_new()
    assert file.Sts[0].comment == " Energy =        -1168.36168282 Eh        #Cluster:     1"
    assert file.Sts[53].comment == " Energy =        -1168.35637181 Eh        #Cluster:     54"


def test_method_comment_keep():
    file = ClassGeometryXYZs(Path("tests/data/crest_conformers.xyz"))
    file.method_read_xyz()
    file.method_comment_keep()
    assert file.Sts[0].comment == " Energy =        -1168.36168282 Eh        #Cluster:     0"
    assert file.Sts[53].comment == " Energy =        -1168.35637181 Eh        #Cluster:     0"


def test_method_save_xyz():
    file = ClassGeometryXYZs(Path("tests/data/crest_conformers.xyz"))
    file.method_read_xyz()
    file.method_comment_new()
    outfile = Path("tests/compare/output.xyz")
    file.set_filename(outfile)
    file.method_save_xyz([])
    import filecmp
    dcmp = filecmp.cmp(outfile,
                       'tests/compare/xyzfile-1.xyz')
    assert (dcmp == True)
    import os
    os.remove(outfile)


#
#
#
#
#
# def example():
#    print("foo")
#
#
# def test_spam(capsys):
#    example()
#    captured = capsys.readouterr()
#    assert captured.out == "foo\n"
