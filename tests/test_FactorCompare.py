#!/usr/bin/env python
import os
import pytest
from pathlib import Path
import argparse
import censo_ext.FactorCompare as FactorCompare
import filecmp
import platform
in_file = f"tests/data/crest_conformers.xyz"
_system = platform.system()


def test_FactorCompare_miss_args():
    x: dict = {}
    args = argparse.Namespace(**x)
    with pytest.raises(SystemExit) as e:
        FactorCompare.main(args)
    assert e.type == SystemExit
    assert e.value.code == 2  # for argparse error


@pytest.mark.skipif(_system == "Darwin", reason="Crest only work under Linux")
def test_FactorCompare():
    x: dict = {"file": [in_file, in_file]}
    args = argparse.Namespace(**x)
    Result = FactorCompare.main(args)

    if _system == "Darwin":
        compare = ""
    elif _system == "Linux":
        compare = ""
    else:
        pytest.raises(
            ValueError, match="OS system only can run under Darwin or Linux")
    assert Result == None
