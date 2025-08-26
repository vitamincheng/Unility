#!/usr/bin/env python
import pytest
from pathlib import Path
import argparse
import censo_ext.FactorCompare as FactorCompare
import platform
inFile: Path = Path("tests/data/crest_conformers.xyz")
_system = platform.system()


def test_FactorCompare_miss_args():
    x: dict = {}
    with pytest.raises(SystemExit) as e:
        FactorCompare.main(argparse.Namespace(**x))
    assert e.type is SystemExit
    assert e.value.code == 2  # for argparse error


@pytest.mark.skipif(_system == "Darwin", reason="Crest only work under Linux")
def test_FactorCompare():
    x: dict = {"file": [inFile, inFile]}
    Res = FactorCompare.main(argparse.Namespace(**x))

    if _system == "Darwin":
        compare = ""
    elif _system == "Linux":
        compare = ""  # noqa: F841
    else:
        pytest.raises(
            ValueError, match="OS system only can run under Darwin or Linux")
    assert not Res
