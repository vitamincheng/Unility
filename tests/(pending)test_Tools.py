#!/usr/bin/env python3
if __name__ == "__main__":

    import argparse
    from icecream import ic
    from pathlib import Path
    content: dict[int, str] = {1: "anmrfile.py",
                               2: "calculate_rmsd.py",
                               3: "ml4nmr.py",
                               4: "Parameter.py",
                               5: "qm.py",
                               6: "spectra.py",
                               7: "symmetry.py",
                               8: "topo.py",
                               9: "unility",
                               10: "xyzfile.py",
                               99: "quit"}

    ic(content)
    nCH: int = int(input())
    while nCH != 99:

        print("="*80)
        print(content[nCH])
        print("="*80)

        ################################################
        # anmrfile.py
        ################################################
        if nCH == 1:

            from censo_ext.Tools.anmrfile import *
            anmr = Anmr(Path("../tests/04.Hydrogen"))
            anmr.method_read_anmrrc()
            anmr.method_print_anmrrc()
            anmr.method_read_nucinfo()
            anmr.method_print_nucinfo()
            anmr.method_read_enso()
            anmr.method_print_enso()
            anmr.method_read_folder_orcaSJ()
            print(""*3)
            # anmr.method_update_equivalent_orcaSJ()
            # print(""*3)
            # anmr.method_average_orcaSJ()
            # anmr.method_save_average_CONF()
            # something wrong in orcaJcoups not use it!!!
            ##

        #########################################
        #  calculate_rmsd.py
        #########################################

        if nCH == 2:
            # Example for main
            # main(["crest_conformers.xyz","crest_conformers2.xyz","-nh","-q"])
            # main_module_xyz(["crest_conformers.xyz","crest_conformers2.xyz","-nh","-q"])

            # Example for main_xyz
            from censo_ext.Tools.calculate_rmsd import *

            xyzfile = GeometryXYZs(
                Path("../tests/crest_conformers.xyz"))
            xyzfile.method_read_xyz()

            x: dict = {"remove_idx": None, "add_idx": None, "bond_broken": None,
                       "ignore_Hydrogen": False, "debug": False}
            idx_atom1 = cal_rmsd_xyz(
                xyzfile, 1, 2, args=argparse.Namespace(**x))
            ic(idx_atom1)

            x = {"remove_idx": None, "add_idx": None, "bond_broken": [
                52, 55], "ignore_Hydrogen": True, "debug": False}
            idx_atom1 = cal_rmsd_xyz(
                xyzfile, 1, 2, args=argparse.Namespace(**x))
            ic(idx_atom1)

        ###############################################################################
        # ml4nmr.py
        ###############################################################################

        if nCH == 3:

            from censo_ext.Tools.ml4nmr import *

            mol, neighbors = read_mol_neighbors(
                Path("../tests/crest_conformers.xyz"))
            ic(mol, neighbors)
            mol, neighbors, bond_order = read_mol_neighbors_bond_order(
                Path("../tests/crest_conformers.xyz"))
            ic([atom for atom in mol])
            ic([atom.number for atom in mol])  # type: ignore
            ic(mol, neighbors, bond_order)

        ###############################################################################
        # Parameter.py
        ###############################################################################

        if nCH == 4:

            from censo_ext.Tools.Parameter import ELEMENT_NAMES, ELEMENT_WEIGHTS, NAMES_ELEMENT

            ic(ELEMENT_NAMES)
            ic(ELEMENT_WEIGHTS)
            ic(NAMES_ELEMENT)

        ###############################################################################
        # qm.py
        ###############################################################################

        if nCH == 5:

            from censo_ext.Tools.qm import *
            x = {"out": "output.dat", "start": -
                 0.5, "end": 10.5, "lw": 1, "mf": 500.0, "cutoff": 0.001, "debug": False}
            v = [964, 2775.76, 2768.20, 928]

            J = np.array([[0.0,   0.0,   0.0,   0.0],
                          [0.0,   0.0, 16.97,   0.0],
                          [0.0, 16.97,   0.0,   7.0],
                          [0.0,   0.0,   7.0,   0.0]])

            # J = np.array([[   0,   6.8, -1.7,  1.2 ],
            #              [ 6.8,     0,16.97, -1.97],
            #              [-1.7, 16.97,    0,  6.43],
            #              [ 1.2, -1.97, 6.43,     0]])

            # R_peak = qm_full(v=v, J=J, nIntergals=1, args=argparse.Namespace(**x))

            R_peak = qm_partial(v=v, J=J, nIntergals=1, idx0_nspins=1,
                                args=argparse.Namespace(**x))
            ic(R_peak)
            print_plot(inpeaklist=R_peak, dpi=10000, nIntergals=2,
                       Active_range=10, args=argparse.Namespace(**x), hidden=False)

        ###############################################################################
        # spectra.py
        ###############################################################################

        if nCH == 6:

            from censo_ext.Tools.spectra import *
            a: np.ndarray = np.array(
                [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 22, 39, 78, 12, 4, 4, 9, -10, -12])
            ic(numpy_threshold_10(a))
            ic(numpy_threshold_3(a))
            ic(numpy_threshold_mean_3(a))
            ic(find_nearest(a, 24))

        ###############################################################################
        # symmetry.py
        ###############################################################################

        if nCH == 7:

            from censo_ext.Tools.symmetry import *
            from censo_ext.Tools.xyzfile import GeometryXYZs
            xyz = GeometryXYZs(Path("../tests/crest_conformers.xyz"))
            xyz.method_read_xyz()
            for idx in range(len(xyz)):
                pos = ([a.tolist() for a in xyz.Sts[idx].coord])
                sym = ([a for a in xyz.Sts[idx].names.values()])
                pg1 = PointGroup(pos, sym)
                ic(pg1.get_point_group())

        ###############################################################################
        # topo.py
        ########################################################

        if nCH == 8:

            from censo_ext.Tools.topo import *
            print(" Bonding : ")
            x = {"file": "../tests/crest_conformers.xyz",
                 "bonding": 51, "print": True, "debug": False}
            Sts_topo: Topo = Topo(x["file"])
            Sts_topo.method_bonding(argparse.Namespace(**x))
            print("")

            print(" Broken_bond_H : ")
            x = {"file": "../tests/crest_conformers.xyz",
                 "bond_broken": [40, 44], "print": True, "debug": False}
            Sts_topo.method_broken_bond_H(argparse.Namespace(**x))

            print(" Broken_bond : ")
            x = {"file": "../tests/crest_conformers.xyz",
                 "bond_broken": [52, 55], "print": True, "debug": False}
            Sts_topo.method_broken_bond(argparse.Namespace(**x))

            mol, neighbors, circle_Mols, residual_Mols = Sts_topo.topology()
            ic(mol)
            ic(neighbors)
            ic(circle_Mols)
            ic(residual_Mols)
            ic(Sts_topo.get_cn())

        ###############################################################################
        # unility.py
        ###############################################################################

        if nCH == 9:
            infile = Path("../tests/crest_conformers.xyz")
            from censo_ext.Tools.utility import IsExist, IsExist_return_bool, program_IsExist
            ic(IsExist(infile))
            ic(IsExist_return_bool(infile))
            ic(program_IsExist("xtb"))

        ###############################################################################
        # xyzfile.py
        #####################################################

        if nCH == 10:

            from censo_ext.Tools.xyzfile import *
            infile = GeometryXYZs()
            infile.set_filename(Path("../tests/crest_conformers.xyz"))
            infile.method_read_xyz()
            infile.method_print([])
            infile.Sts[0].method_update_masses()
            print(infile.Sts[0].method_computeCOM())
            print(infile.Sts[0].method_computeInertia())

            infile.method_idx_molecules_xyzs(idx1=1)

        ic(content)
        nCH = int(input())
