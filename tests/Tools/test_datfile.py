from pathlib import Path
import os
import numpy as np
import pytest
import filecmp
from censo_ext.Tools.datfile import CensoDat, Peaks_npz, unit_conversion


def test_datfile_miss_args() -> None:
    with pytest.raises(SystemExit) as e:
        Missing = CensoDat()                # noqa: F841
    assert e.type is SystemExit
    assert e.value.code == 0


def test_datfile_Censo() -> None:
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


def test_unit_conversion() -> None:
    uc = unit_conversion(np.array([1.0, 2.0, 3.0, 5.0]))                # noqa: F841
    assert uc.index(5.0) == 3
    assert uc.ppm(1) == 2.0
    result = uc.ppm_scale()
    expected = np.array([1.0, 2.0, 3.0, 5.0])
    np.testing.assert_array_equal(result, expected)
    assert uc.ppm_limits() == (1.0, 5.0)


def test_Peaks_npz():
    uc = unit_conversion(np.array([1.0, 2.0, 3.0, 5.0]))                # noqa: F841
    npz = Peaks_npz(uc)
    npz.method_load_Data(
        [(1, 0.5, 0.4, 400), (2, 2.5, 2.4, 800), (5, 1.8, 1.6, 100)])
    assert len(npz) == 3
    assert npz.method_ppm2cID(2.44) == 2
    npz.method_merge_cID([2, 5])
    assert len(npz) == 2
