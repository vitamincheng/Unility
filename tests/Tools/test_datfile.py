from pathlib import Path
import os
import pytest
import filecmp
from censo_ext.Tools.datfile import CensoDat


def test_anmrfile_miss_args() -> None:
    Missing: CensoDat
    with pytest.raises(FileNotFoundError) as e:
        Missing = CensoDat()                    # noqa: F841
    assert str(e.value) == "anmr.dat The file is not Exist ..."


def test_anmrfile_Censo() -> None:
    infile: CensoDat = CensoDat(
        Path("tests/data/34.Ergocalciferol/04.Hydrogen/anmr.dat"))
    assert len(infile) == 77868
    assert infile.get_fileName() == Path(
        "tests/data/34.Ergocalciferol/04.Hydrogen/anmr.dat")
    infile.method_normalize_dat()
    outfile = Path("tests/compare/out.dat")
    source = Path("tests/compare/anmr_normal.dat")
    infile.set_fileName(outfile)
    infile.method_save_dat()

    assert filecmp.cmp(outfile, source)
    os.remove(outfile)
