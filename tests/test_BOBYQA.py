#!/usr/bin/env python3
import pytest
from pathlib import Path
import argparse
import censo_ext.BOBYQA as BOBYQA
import censo_ext.anmr as anmr
import filecmp
import platform

_system = platform.system()
DirName: Path = Path(
    "tests/data/31.Cyclohexanone/03.Censo_for_Hydrogen_(revTPSS)")
RefDat: Path = Path("1r_h.dat")


def anmr_BOBYQA_init():
    # Create Average/NMR/orcaS.out and Average/NMR/orcaS_BOBYQA.out file
    args_x: dict = {"auto": True, "average": False, "dir": DirName,
                    "bobyqa": False, "mf": 500, "lw": None, "thr": None, "json": None, "thrab": 0.025,
                    "tb": 4, "mss": 9, "cutoff": 0.001, "show": False, "start": None, "end": None, "out": "output.dat"}
    args = argparse.Namespace(**args_x)
    anmr.main(args)

    args_y: dict = {"dir": DirName, "ref": RefDat, "limit": 0.20, "prog": None}
    args2 = argparse.Namespace(**args_y)
    assert BOBYQA.main(args2) == (False, True)


def BOBYQA_final_remove_files():
    from censo_ext.Tools.utility import delete_all_files
    delete_all_files(DirName / "output.dat",
                     DirName / "peaks.json",
                     DirName / "Average/NMR/orcaA.out",
                     DirName / "Average/NMR/orcaJ.out",
                     DirName / "Average/NMR/orcaS.out",
                     DirName / "Average/NMR/orcaS-BOBYQA.out",)


def test_BOBYQA_cregen_miss_args():
    args_x: dict = {}
    args = argparse.Namespace(**args_x)
    with pytest.raises(SystemExit) as e:
        BOBYQA.main(args)
    assert e.type is SystemExit
    assert e.value.code == 2  # for argparse error


def test_BOBYQA_block_file():
    # run block orcaS_BOBYQA file
    anmr_BOBYQA_init()
    args_y: dict = {"dir": DirName, "ref": RefDat, "limit": 0.20, "prog": None}
    args2 = argparse.Namespace(**args_y)
    assert BOBYQA.main(args2) == (True, False)


def test_BOBYQA_block_file_final_remove_files():
    BOBYQA_final_remove_files()


def test_BOBYQA_single():
    anmr_BOBYQA_init()
    import shutil
    shutil.copyfile(DirName / Path("orcaS-BOBYQA.out"),
                    DirName / Path("Average/NMR/orcaS-BOBYQA.out"))
    args_y: dict = {"dir": DirName, "ref": RefDat, "limit": 0.20, "prog": None}
    args2 = argparse.Namespace(**args_y)
    assert BOBYQA.main(args2) == (True, False)
    assert filecmp.cmp(DirName / Path("output.dat"),
                       DirName / Path("orcaS-BOBYQA-anmrpy.dat"))
    BOBYQA_final_remove_files()


@pytest.mark.skipif(_system == "Darwin", reason="anmr only work under linux CLI")
def test_BOBYQA_single_external_prog(monkeypatch):
    anmr_BOBYQA_init()
    import shutil
    shutil.copyfile(DirName / Path("orcaS-BOBYQA.out"),
                    DirName / Path("Average/NMR/orcaS-BOBYQA.out"))

    args_y: dict = {"dir": DirName, "ref": RefDat, "limit": 0.20, "prog": True}
    from io import StringIO
    monkeypatch.setattr('sys.stdin', StringIO('Y\n'))
    args2 = argparse.Namespace(**args_y)
    assert BOBYQA.main(args2) == (True, False)
    assert filecmp.cmp(DirName / Path("anmr.dat"),
                       DirName / Path("orcaS-BOBYQA-anmr.dat"))
    BOBYQA_final_remove_files()
    from censo_ext.Tools.utility import delete_all_files
    delete_all_files(DirName / "anmr.dat",
                     DirName / "anmr.out",
                     DirName / "tmpanmr.1",
                     DirName / "tmpanmr_frag.av",
                     DirName / "tmpanmr_full.av")


def test_BOBYQA_group():
    anmr_BOBYQA_init()
    import shutil
    shutil.copyfile(DirName / Path("orcaS-BOBYQA-group.out"),
                    DirName / Path("Average/NMR/orcaS-BOBYQA.out"))
    args_y: dict = {"dir": DirName, "ref": RefDat, "limit": 0.20, "prog": None}
    args2 = argparse.Namespace(**args_y)
    assert BOBYQA.main(args2) == (True, False)
    assert filecmp.cmp(DirName / Path("output.dat"),
                       DirName / Path("orcaS-BOBYQA-group-anmrpy.dat"))
    BOBYQA_final_remove_files()


@pytest.mark.skipif(_system == "Darwin", reason="anmr only work under linux CLI")
def test_BOBYQA_group_external_prog():
    raise NotImplementedError("pending")


if __name__ == "__main__":
    import cProfile
    cProfile.run("test_BOBYQA_init()", sort="cumtime")
    cProfile.run("test_BOBYQA_single()", sort="cumtime")

#   on Unility directory
#   pytest -v tests/test_BOBYQA.py::test_BOBYQA_block_file
#   pytest -v tests/test_BOBYQA.py::test_BOBYQA_single_external_prog
#   pytest -v tests/test_BOBYQA.py::test_BOBYQA_group_external_prog
#   pytest -v --pdb tests/test_BOBYQA.py::test_BOBYQA_single_external_prog   (--pdb is debug)
#
