#!/usr/bin/env python3
import os
import string
import pytest
from pathlib import Path
import argparse
import censo_ext.datNormalized as datNormalized
import filecmp
import platform


def test_datNormalized_miss_args():
    x: dict = {}
    args = argparse.Namespace(**x)
    with pytest.raises(SystemExit) as e:
        datNormalized.main(args)
    assert e.type == SystemExit
    assert e.value.code == 2  # for argparse error
