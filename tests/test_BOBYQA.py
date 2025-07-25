
#!/usr/bin/env python3
import os
import pytest
from pathlib import Path
import argparse
import censo_ext.BOBYQA as BOBYQA
import censo_ext.anmr as anmr
import filecmp
import platform

_system = platform.system()


def test_cregen_miss_args():
    x: dict = {}
    args = argparse.Namespace(**x)
    with pytest.raises(SystemExit) as e:
        BOBYQA.main(args)
    assert e.type == SystemExit
    assert e.value.code == 2  # for argparse error


# @pytest.mark.skipif(_system == "Darwin", reason="crest only work under linux")
def test_BOBYQA():
    x: dict = {"auto": True, "average": False, "dir": "tests/data/31.Cyclohexanone/03.Censo_for_Hydrogen_(revTPSS)",
               "bobyqa": False, "mf": 500, "lw": None, "thr": None, "json": None, "thrab": 0.025,
               "tb": 4, "mss": 9, "cutoff": 0.001, "show": False, "start": None, "end": None, "out": "output.dat"}
    args = argparse.Namespace(**x)
    anmr.main(args)

    x2: dict = {
        "dir": "tests/data/31.Cyclohexanone/03.Censo_for_Hydrogen_(revTPSS)", "ref": "1r_h.dat", "limit": 0.20}
    args2 = argparse.Namespace(**x2)
    BOBYQA.main(args2)
    # Create orcaS_BOBYQA file
    BOBYQA.main(args2)
    # run single peaks or multipeaks

    # compare = ""
    # if _system == "Linux":
    #    compare = f"tests/compare/cregen_cluster.xyz"
    # elif _system == "Darwin":  # No crest under Darwin system
    #    compare = f"tests/compare/cregen_cluster_Darwin.xyz"
    # else:
    #    pytest.raises(
    #        ValueError, match="OS system only can run under Darwin or Linux")
    # assert filecmp.cmp(args.out, compare) == True
    # os.remove(args.out)
    # os.remove("isomers.out")

# BOBYQA.py -d tests/data/31.Cyclohexanone/03.Censo_For_Hydorgen(revTPSS) -r 1r_h.dat
