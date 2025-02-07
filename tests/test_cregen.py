#!/usr/bin/env python3
import os
import string
import pytest
from pathlib import Path
import argparse
import censo_ext.cregen as cregen
import filecmp
import platform


def test_cregen_miss_args():
    x: dict = {}
    args = argparse.Namespace(**x)
    with pytest.raises(SystemExit) as e:
        cregen.main(args)
    assert e.type == SystemExit
    assert e.value.code == 2  # for argparse error


def test_cregen_miss_crest():
    x: dict = {"file": "tests/data/crest_conformers2.xyz", "rthr": 0.175, "bthr": 0.03,
               "ethr": 0.15, "ewin": 4, "out": "cluster.xyz"}
    args = argparse.Namespace(**x)

    if platform.system() != "Darwin":
        with pytest.raises(SystemExit) as e:
            cregen.main(args)
        assert e.type == SystemExit
        assert e.value.code == 2  # for argparse error
