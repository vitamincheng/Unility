import pathlib
import argparse
from pathlib import Path
import pytest
from censo_ext.Tools.anmrfile import Anmrrc
from censo_ext.Tools.anmrfile import ClassAnmr


def test_method_anmrrc():
    # For Hydrogen
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

    # For Carbon
    file = ClassAnmr(Path("tests/data/07.Carbon"))
    file.method_read_anmrrc()
    import sys
    original_stdout = sys.stdout
    filename = "tests/compare/.anmrrctest"
    with open(filename, "w") as f:
        sys.stdout = f
        file.method_print_anmrrc()
    sys.stdout = original_stdout

    import filecmp
    source = "tests/data/07.Carbon/.anmrrc"
    dcmp = filecmp.cmp(filename, source)
    assert (dcmp == True)
    import os
    os.remove(filename)
    path = file.get_Directory()
    assert path == Path("tests/data/07.Carbon")
    active = file.get_Anmr_Active()
    assert active == ['C']


def test_get_average_orcaSJ_Exist():
    # For Hydrogen
    file = ClassAnmr(Path("tests/data/04.Hydrogen"))
    exist = file.get_average_orcaSJ_Exist()
    assert exist == True
    # for Carbon
    file = ClassAnmr(Path("tests/data/07.Carbon"))
    exist = file.get_average_orcaSJ_Exist()
    assert exist == True


def test_method_read_enso():
    # For Hydrogen
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

    # for Carbon
    file = ClassAnmr(Path("tests/data/07.Carbon"))
    file.method_read_enso()

    import sys
    original_stdout = sys.stdout
    filename = "tests/compare/ensotest"
    with open(filename, "w") as f:
        sys.stdout = f
        file.method_print_enso()
    sys.stdout = original_stdout

    import filecmp
    source = "tests/data/07.Carbon/anmr_enso"
    dcmp = filecmp.cmp(filename, source)
    assert (dcmp == True)
    import os
    os.remove(filename)
