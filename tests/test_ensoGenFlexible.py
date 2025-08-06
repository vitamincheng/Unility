#!/usr/bin/env python
import pytest
import argparse
import censo_ext.ensoGenFlexible as ensoGenFlexible


def test_ensoGenFlexible_miss_args():
    x: dict = {}
    args = argparse.Namespace(**x)
    with pytest.raises(SystemExit) as e:
        ensoGenFlexible.main(args)
    assert e.type is SystemExit
    assert e.value.code == 2  # for argparse error
