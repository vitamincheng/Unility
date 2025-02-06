#!/usr/bin/env python3
import os
import string
import pytest
from pathlib import Path
import argparse
import censo_ext.molclus_orca as orca
import filecmp
import platform


def test_orca_miss_args():
    x: dict = {}
    args = argparse.Namespace(**x)
    with pytest.raises(SystemExit) as e:
        orca.main(args)
    assert e.type == SystemExit
    assert e.value.code == 2  # for argparse error
