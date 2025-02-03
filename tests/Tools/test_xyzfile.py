from pathlib import Path
import pytest
from censo_ext.Tools.xyzfile import GeometryXYZs
import os
import filecmp


def test_method_read_xyz():
    # For crest_conformers.xyz
    file = GeometryXYZs(Path("tests/data/crest_conformers.xyz"))
    assert len(file) == 0
    file.method_read_xyz()
    assert len(file) == 54
    assert file.Sts[0].get_comment_energy() == -1168.36168282
    assert file.Sts[-1].get_comment_energy() == -1168.35637181
    assert len(file.get_comment_energy()) == 54
    assert file.get_comment_energy()[0] == -1168.36168282
    assert file.get_comment_energy()[-1] == -1168.35637181
    # For isomers.xyz
    file = GeometryXYZs(Path("tests/data/isomers.xyz"))
    assert len(file) == 0
    file.method_read_xyz()
    assert len(file) == 3
    assert file.Sts[0].get_comment_energy() == -309.820981116391
    assert file.Sts[-1].get_comment_energy() == -309.815154391997


def test_method_rewrite_comment():
    # for Crest_conformers.xyz
    file = GeometryXYZs(Path("tests/data/crest_conformers.xyz"))
    file.method_read_xyz()
    file.method_rewrite_comment()
    assert file.Sts[0].comment == " Energy =        -1168.36168282 Eh        #Cluster:     0"
    assert file.Sts[-1].comment == " Energy =        -1168.35637181 Eh        #Cluster:     0"
    # for isomers.xyz
    file = GeometryXYZs(Path("tests/data/isomers.xyz"))
    file.method_read_xyz()
    file.method_rewrite_comment()
    assert file.Sts[0].comment == " Energy =        -309.820981116391 Eh        #Cluster:     1"
    assert file.Sts[-1].comment == " Energy =        -309.815154391997 Eh        #Cluster:     3"


def test_method_comment_new():
    # for crest_conformers.xyz
    file = GeometryXYZs(Path("tests/data/crest_conformers.xyz"))
    file.method_read_xyz()
    file.method_comment_new()
    assert file.Sts[0].comment == " Energy =        -1168.36168282 Eh        #Cluster:     1"
    assert file.Sts[-1].comment == " Energy =        -1168.35637181 Eh        #Cluster:     54"
    # for isomers.xyz
    file = GeometryXYZs(Path("tests/data/isomers.xyz"))
    file.method_read_xyz()
    file.method_comment_new()
    assert file.Sts[0].comment == " Energy =        -309.820981116391 Eh        #Cluster:     1"
    assert file.Sts[-1].comment == " Energy =        -309.815154391997 Eh        #Cluster:     3"


def test_method_comment_keep():
    # for crest_conformers.xyz
    file = GeometryXYZs(Path("tests/data/crest_conformers.xyz"))
    file.method_read_xyz()
    file.method_comment_keep()
    assert file.Sts[0].comment == " Energy =        -1168.36168282 Eh        #Cluster:     0"
    assert file.Sts[-1].comment == " Energy =        -1168.35637181 Eh        #Cluster:     0"
    # for isomers.xyz
    file = GeometryXYZs(Path("tests/data/isomers.xyz"))
    file.method_read_xyz()
    file.method_comment_keep()
    assert file.Sts[0].comment == " Energy =        -309.820981116391 Eh        #Cluster:     1"
    assert file.Sts[-1].comment == " Energy =        -309.815154391997 Eh        #Cluster:     3"


def test_method_save_xyz():
    # for crest_conformers.xyz
    file = GeometryXYZs(Path("tests/data/crest_conformers.xyz"))
    file.method_read_xyz()
    file.method_comment_new()
    outfile = Path("tests/compare/output.xyz")
    file.set_filename(outfile)
    file.method_save_xyz([])
    dcmp = filecmp.cmp(outfile,
                       'tests/compare/xyzfile-1.xyz')
    assert (dcmp == True)
    os.remove(outfile)
    # for isomers.xyz
    file = GeometryXYZs(Path("tests/data/isomers.xyz"))
    file.method_read_xyz()
    file.method_comment_new()
    outfile = Path("tests/compare/output.xyz")
    file.set_filename(outfile)
    file.method_save_xyz([])
    dcmp = filecmp.cmp(outfile,
                       'tests/compare/xyzfile-2.xyz')
    assert (dcmp == True)
    os.remove(outfile)
