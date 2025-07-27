#!/usr/bin/env python3
import os
import pytest
from pathlib import Path
import argparse
import censo_ext.cregen as cregen
import filecmp
import platform

_system = platform.system()


def test_cregen_miss_args():
    x: dict = {}
    args = argparse.Namespace(**x)
    with pytest.raises(SystemExit) as e:
        cregen.main(args)
    assert e.type == SystemExit
    assert e.value.code == 2  # for argparse error


@pytest.mark.skipif(_system == "Darwin", reason="crest only work under linux")
def test_cregen_crest():
    x: dict = {"file": "tests/data/crest_conformers2.xyz", "rthr": 0.175, "bthr": 0.03,
               "ethr": 0.15, "ewin": 4, "out": "cluster.xyz"}
    args = argparse.Namespace(**x)
    cregen.main(args)

    compare: Path = Path("")
    if _system == "Linux":
        compare = Path("tests/compare/cregen_cluster.xyz")
    elif _system == "Darwin":  # No crest under Darwin system
        compare = Path("tests/compare/cregen_cluster_Darwin.xyz")
    else:
        pytest.raises(
            ValueError, match="OS system only can run under Darwin or Linux")
    assert filecmp.cmp(args.out, compare) == True
    os.remove(args.out)
    os.remove("isomers.out")
