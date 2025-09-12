from pathlib import Path
import os
import pytest
import filecmp
from censo_ext.Tools.datfile import CensoDat


def test_anmrfile_miss_args() -> None:
    Missing: CensoDat
    fileName: Path = Path("anmr.dat")  # default FileName
    with pytest.raises(FileNotFoundError) as e:
        Missing = CensoDat()                    # noqa: F841
    assert str(e.value) == f"  The file {fileName} is not Exist ..."


def test_anmrfile_Censo() -> None:
    anmrFile: Path = Path("tests/data/34.Ergocalciferol/04.Hydrogen/anmr.dat")
    infile: CensoDat = CensoDat(anmrFile)
    assert len(infile) == 77868
    assert infile.get_fileName() == Path(anmrFile)
    infile.method_normalize_dat()
    outfile = Path("tests/compare/out.dat")
    source = Path("tests/compare/anmr_normal.dat")
    infile.set_fileName(outfile)
    infile.method_save_dat()

    assert filecmp.cmp(outfile, source)
    os.remove(outfile)
