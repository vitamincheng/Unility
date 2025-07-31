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


def test_BOBYQA_cregen_miss_args():
    args_x: dict = {}
    args = argparse.Namespace(**args_x)
    with pytest.raises(SystemExit) as e:
        BOBYQA.main(args)
    assert e.type == SystemExit
    assert e.value.code == 2  # for argparse error


def BOBYQA_init():
    # Create orcaS_BOBYQA file
    args_x: dict = {"auto": True, "average": False, "dir": DirName,
                    "bobyqa": False, "mf": 500, "lw": None, "thr": None, "json": None, "thrab": 0.025,
                    "tb": 4, "mss": 9, "cutoff": 0.001, "show": False, "start": None, "end": None, "out": "output.dat"}
    args = argparse.Namespace(**args_x)
    anmr.main(args)

    args_y: dict = {"dir": DirName, "ref": RefDat, "limit": 0.20, "prog": None}
    args2 = argparse.Namespace(**args_y)

    Res_init: tuple[bool, bool] = BOBYQA.main(args2)
    assert Res_init == (False, True)


def test_BOBYQA_block_file():
    # run block orcaS_BOBYQA file
    BOBYQA_init()
    args_y: dict = {"dir": DirName, "ref": RefDat, "limit": 0.20, "prog": None}
    args2 = argparse.Namespace(**args_y)
    Res: tuple[bool, bool] = BOBYQA.main(args2)
    assert Res == (True, False)
    BOBYQA_final_remove_files()


def test_BOBYQA_single():
    BOBYQA_init()
    import shutil
    shutil.copyfile(DirName / Path("orcaS-BOBYQA.out"),
                    DirName / Path("Average/NMR/orcaS-BOBYQA.out"))
    args_y: dict = {"dir": DirName, "ref": RefDat, "limit": 0.20, "prog": None}
    args2 = argparse.Namespace(**args_y)
    Res: tuple[bool, bool] = BOBYQA.main(args2)
    assert Res == (True, False)
    assert filecmp.cmp(DirName / Path("output.dat"),
                       DirName / Path("output-BOBYQA-anmrpy.dat")) == True
    BOBYQA_final_remove_files()


@pytest.mark.skipif(_system == "Darwin", reason="crest only work under linux")
def test_BOBYQA_single_external_prog():
    BOBYQA_init()
    import shutil
    shutil.copyfile(DirName / Path("orcaS-BOBYQA.out"),
                    DirName / Path("Average/NMR/orcaS-BOBYQA.out"))

    args_y: dict = {"dir": DirName, "ref": RefDat, "limit": 0.20, "prog": True}
    args2 = argparse.Namespace(**args_y)
    Res: tuple[bool, bool] = BOBYQA.main(args2)
    assert Res == (True, False)
    assert filecmp.cmp(DirName / Path("anmr.dat"),
                       DirName / Path("output-BOBYQA-anmr.dat")) == True
    BOBYQA_final_remove_files()


def test_BOBYQA_group():
    raise NotImplementedError("pending")
    BOBYQA_init()
    import shutil
    shutil.copyfile(DirName / Path("orcaS-BOBYQA-group.out"),
                    DirName / Path("Average/NMR/orcaS-BOBYQA.out"))
    args_y: dict = {"dir": DirName, "ref": RefDat, "limit": 0.20, "prog": None}
    args2 = argparse.Namespace(**args_y)
    Res: tuple[bool, bool] = BOBYQA.main(args2)
    assert Res == (True, False)
    assert filecmp.cmp(DirName / Path("output.dat"),
                       DirName / Path("output-BOBYQA-anmrpy.dat")) == True
    BOBYQA_final_remove_files()


def BOBYQA_final_remove_files():
    from censo_ext.Tools.utility import delete_all_files
    delete_all_files(DirName / "output.dat",
                     DirName / "peaks.json",
                     DirName / "Average/NMR/orcaA.out",
                     DirName / "Average/NMR/orcaJ.out",
                     DirName / "Average/NMR/orcaS.out",
                     DirName / "Average/NMR/orcaS-BOBYQA.out",)


if __name__ == "__main__":
    import cProfile
    cProfile.run("test_BOBYQA_init()", sort="cumtime")
    cProfile.run("test_BOBYQA_single()", sort="cumtime")
