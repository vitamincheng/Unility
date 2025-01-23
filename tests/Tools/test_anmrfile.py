import pathlib
import argparse
from pathlib import Path
import pytest
from Tools.anmrfile import Anmrrc
from Tools.anmrfile import ClassAnmr


def test_method_anmrrc():
    file = ClassAnmr(Path("tests/data/04.Hydrogen"))
    file.method_read_anmrrc()
    import sys
    original_stdout = sys.stdout
    filename = "tests/compare/.anmrrctest"
    with open(filename, "w") as f:
        sys.stdout = f
        file.method_print_anmrrc()
    sys.stdout = original_stdout

    import filecmp
    source = "tests/data/04.Hydrogen/.anmrrc"
    dcmp = filecmp.cmp(filename, source)
    assert (dcmp == True)
    import os
    os.remove(filename)
    path = file.get_Directory()
    assert path == Path("tests/data/04.Hydrogen")
    active = file.get_Anmr_Active()
    assert active == ['H']


def test_get_average_orcaSJ_Exist():
    file = ClassAnmr(Path("tests/data/04.Hydrogen"))
    exist = file.get_average_orcaSJ_Exist()
    assert exist == True


def test_method_read_enso():
    file = ClassAnmr(Path("tests/data/04.Hydrogen"))
    file.method_read_enso()

    import sys
    original_stdout = sys.stdout
    filename = "tests/compare/ensotest"
    with open(filename, "w") as f:
        sys.stdout = f
        file.method_print_enso()
    sys.stdout = original_stdout

    import filecmp
    source = "tests/data/04.Hydrogen/anmr_enso"
    dcmp = filecmp.cmp(filename, source)
    assert (dcmp == True)
    import os
    os.remove(filename)

# def test_method_read_nucinfo():
#    file = ClassAnmr(Path("tests/04.Hydrogen"))
#    file.method_read_nucinfo()
#    import sys
#    original_stdout = sys.stdout
#    filename = "tests/04.Hydrogen/.nucinfo"
#    with open(filename, "w") as f:
#        sys.stdout = f
#        file.method_print_nucinfo()
#    sys.stdout = original_stdout
#
#    import filecmp
#    source = "tests/04.Hydrogen/anmr_nucinfo"
#    dcmp = filecmp.cmp(filename, source)
#    assert (dcmp == True)
#    import os
#    os.remove(filename)

    # file = ClassGeometryXYZs(Path("tests/crest_conformers.xyz"))
    # assert len(file) == 0
    # file.method_read_xyz()
    # assert len(file) == 54
    # assert file.Sts[0].get_comment_energy() == -1168.36168282
    # assert file.Sts[53].get_comment_energy() == -1168.35637181
