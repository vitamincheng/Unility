import pathlib
import argparse
from pathlib import Path
import pytest
from censo_ext.Tools.anmrfile import Anmr
from censo_ext.Tools.anmrfile import OrcaSJ
import sys
import os
import filecmp


def test_anmrfile_miss_args():

    Missing: OrcaSJ = OrcaSJ()

    with pytest.raises(SystemExit) as e:
        Missing.method_read_orcaJ()
    assert e.type == SystemExit
    assert e.value.code == 1  # for argparse error

    with pytest.raises(SystemExit) as e:
        Missing.method_read_orcaS()
    assert e.type == SystemExit
    assert e.value.code == 1  # for argparse error


def test_anmrfile_read_OrcaSJ():
    Hydrogen: OrcaSJ = OrcaSJ()
    assert Hydrogen.method_read_orcaJ(
        Path("tests/data/EthylAcetate/03.Censo/CONF1/NMR/orcaJ.out"))
    assert Hydrogen.method_read_orcaS(
        Path("tests/data/EthylAcetate/03.Censo/CONF1/NMR/orcaS.out"))


def test_anmrfile_anmrrc():
    # For Hydrogen
    file = Anmr(Path("tests/data/04.Hydrogen"))
    file.method_read_anmrrc()
    original_stdout = sys.stdout
    filename = "tests/compare/.anmrrc"
    with open(filename, "w") as f:
        sys.stdout = f
        file.method_print_anmrrc()
    sys.stdout = original_stdout

    source = "tests/data/04.Hydrogen/.anmrrc"
    dcmp = filecmp.cmp(filename, source)
    assert (dcmp == True)
    os.remove(filename)
    path = file.get_Directory()
    assert path == Path("tests/data/04.Hydrogen")
    active = file.get_Anmr_Active()
    assert active == ['H']

    # For Carbon
    file = Anmr(Path("tests/data/07.Carbon"))
    file.method_read_anmrrc()
    original_stdout = sys.stdout
    filename = "tests/compare/.anmrrctest"
    with open(filename, "w") as f:
        sys.stdout = f
        file.method_print_anmrrc()
    sys.stdout = original_stdout

    source = "tests/data/07.Carbon/.anmrrc"
    dcmp = filecmp.cmp(filename, source)
    assert (dcmp == True)
    os.remove(filename)
    path = file.get_Directory()
    assert path == Path("tests/data/07.Carbon")
    active = file.get_Anmr_Active()
    assert active == ['C']


def test_anmrfile_get_avg_orcaSJ_Exist():
    # For Hydrogen
    file = Anmr(Path("tests/data/04.Hydrogen"))
    exist = file.get_avg_orcaSJ_Exist()
    assert exist == True
    # for Carbon
    file = Anmr(Path("tests/data/07.Carbon"))
    exist = file.get_avg_orcaSJ_Exist()
    assert exist == True


def test_anmrfile_read_enso_miss_args():
    # For Hydrogen
    file = Anmr(Path("tests/data/04.Hydrogen"))
    with pytest.raises(SystemExit) as e:
        file.method_read_enso(Path("anmr_enso_error"))
    assert e.type == SystemExit
    assert e.value.code == 1


def test_anmrfile_read_enso():
    # For Hydrogen
    file = Anmr(Path("tests/data/04.Hydrogen"))
    file.method_read_enso()

    original_stdout = sys.stdout
    filename = "tests/compare/ensotest"
    with open(filename, "w") as f:
        sys.stdout = f
        file.method_print_enso()
    sys.stdout = original_stdout

    source = "tests/data/04.Hydrogen/anmr_enso"
    assert filecmp.cmp(filename, source)
    os.remove(filename)

    # for Carbon
    file = Anmr(Path("tests/data/07.Carbon"))
    file.method_read_enso()

    original_stdout = sys.stdout
    filename = "tests/compare/ensotest"
    with open(filename, "w") as f:
        sys.stdout = f
        file.method_print_enso()
    sys.stdout = original_stdout

    source = "tests/data/07.Carbon/anmr_enso"
    dcmp = filecmp.cmp(filename, source)
    assert (dcmp == True)
    os.remove(filename)


def test_anmrfile_Censo():
    from censo_ext.Tools.anmrfile import CensoDat
    infile: CensoDat = CensoDat(Path("tests/data/04.Hydrogen/anmr.dat"))
    assert len(infile) == 77868
    assert infile.get_fileName() == Path("tests/data/04.Hydrogen/anmr.dat")
    infile.method_normalize_dat()
    outfile = Path("tests/compare/out.dat")
    source = Path("tests/compare/anmr_normal.dat")
    infile.set_fileName(outfile)
    infile.method_save_dat()

    dcmp = filecmp.cmp(outfile, source)
    assert (dcmp == True)
    os.remove(outfile)


def test_anmrfile_read_anmrSJ_H():
    file = Anmr(Path("tests/data/04.Hydrogen"))
    file.method_read_anmrSJ(Path("anmrh0.out"))
    assert len(file.anmrS) == 33
    assert file.anmrS[0] == [1, 4, 1, 2.317]
    assert file.anmrS[-1] == [33, 71, 3, 1.104]
    assert len(file.anmrJ[0]) == 33
    assert file.anmrJ[0][1] == -13.04873
    assert file.anmrJ[-1][-2] == -0.00107


def test_anmrfile_read_anmrSJ_C():
    file = Anmr(Path("tests/data/07.Carbon"))
    file.method_read_anmrSJ(Path("anmrc.out"))
    assert len(file.anmrS) == 28
    assert file.anmrS[0] == [1, 1, 1, 125.371]
    assert file.anmrS[-1] == [28, 61, 1, 9.315]
    assert len(file.anmrJ[0]) == 28
    assert file.anmrJ[0][1] == 0.00000
    assert file.anmrJ[-1][-2] == -0.00000
