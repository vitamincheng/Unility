
#!/usr/bin/env python
import pytest
from pathlib import Path
import argparse
import censo_ext.Pending_dat2BOBYQA as Pending_dat2BOBYQA
import censo_ext.BOBYQA as BOBYQA
import censo_ext.anmr as anmr

DirName: Path = Path(
    "tests/data/31.Cyclohexanone/03.Censo_for_Hydrogen_(revTPSS)")
RefDat: Path = Path("1r_h.dat")


def anmr_init():
    # Create Average/NMR/orcaS.out and Average/NMR/orcaS_BOBYQA.out file
    args_x: dict = {"auto": True, "average": False, "dir": DirName, "bobyqa": False, "mf": 500,
                    "lw": None, "thr": None, "json": None, "thrab": 0.025, "tb": 4, "mss": 9,
                    "cutoff": 0.001, "show": False, "start": None, "end": None, "out": "output.dat"}
    args = argparse.Namespace(**args_x)
    anmr.main(args)


def BOBYQA_init():
    args_y: dict = {"dir": DirName, "ref": RefDat, "limit": 0.20, "prog": None}
    args2 = argparse.Namespace(**args_y)

    with pytest.raises(SystemExit) as e:
        BOBYQA.main(args2)
    assert e.type is SystemExit
    assert e.value.code == 0  # for normal exit(0)


def BOBYQA_final_remove_files():
    from censo_ext.Tools.utility import delete_all_files
    delete_all_files(DirName / "output.dat",
                     DirName / "peaks.json",
                     DirName / "Average/NMR/orcaA.out",
                     DirName / "Average/NMR/orcaJ.out",
                     DirName / "Average/NMR/orcaS.out",
                     DirName / "Average/NMR/orcaS-BOBYQA.out",)


def test_dat2BOBYQA_missing_args():
    x: dict = {}
    args = argparse.Namespace(**x)
    with pytest.raises(SystemExit) as e:
        Pending_dat2BOBYQA.main(args)
    assert e.type is SystemExit
    assert e.value.code == 2  # for argparse error


def test_dat2BOBYQA_Full():
    anmr_init()
    BOBYQA_init()
    x: dict = {"dir": DirName, "file": DirName /
               RefDat, "start": -5.0, "end": 15.0, "dpi": 100, "thr": 10}
    Pending_dat2BOBYQA.main(argparse.Namespace(**x))
    BOBYQA_final_remove_files()
