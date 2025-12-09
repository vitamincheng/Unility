
import pytest
from censo_ext.Tools.factor import method_factor_analysis, method_factor_opt
from pathlib import Path

from censo_ext.Tools.utility import AtomID


def test_factor_analysis_miss_args() -> None:
    fileName: Path = Path("tests/data/crest_conformers000.xyz")

    with pytest.raises(SystemExit) as e:
        method_factor_analysis(inFile=fileName, _factor=0.5)
    assert e.type is SystemExit
    assert e.value.code == 0


def test_factor_analysis() -> None:
    fileName: Path = Path("tests/data/crest_conformers.xyz")
    idx1_minor_factor, Table_S = method_factor_analysis(
        inFile=fileName, _factor=0.50)
    assert idx1_minor_factor == [1, 2, 7, 8, 15, 19, 21,
                                 23, 24, 27, 30, 33, 35, 37, 40, 44, 45, 50, 52, 55]
    assert len(Table_S) == 29
    assert Table_S[AtomID(1)] == pytest.approx(0.41868188431231546)
    assert Table_S[AtomID(61)] == pytest.approx(3.8871420948650774)

    fileName: Path = Path("tests/data/crest_conformers1.xyz")
    idx1_minor_factor, Table_S = method_factor_analysis(
        inFile=fileName, _factor=0.50)
    assert idx1_minor_factor == []
    assert len(Table_S) == 29
    assert Table_S[AtomID(1)] == 0.0
    assert Table_S[AtomID(61)] == 0.0


def test_factor_opt() -> None:
    fileName: Path = Path("tests/data/crest_conformers.xyz")
    idx1_minor_factor, Table_S = method_factor_analysis(
        inFile=fileName, _factor=0.50)

    a0, a1, a2 = method_factor_opt(
        inFile=fileName, _lowFactor=idx1_minor_factor, table_std=Table_S)  # type: ignore
    assert a0
    assert a1 == [52, 55]
    assert a2 == pytest.approx(0.9863091021443952)
