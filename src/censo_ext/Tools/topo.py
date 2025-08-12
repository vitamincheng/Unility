#!/usr/bin/env python
#                                   [07.11.2023] vitamin.cheng@gmail.com
# Need ml4nmr.py
#  ASE library / graph-theory library

import argparse
import numpy as np
import numpy.typing as npt
import censo_ext.Tools.ml4nmr as ml4nmr
from graph import Graph
from ase.atoms import Atoms
from pathlib import Path


class Topo():
    """
    The Topo class is designed to analyze molecular structures using graph theory and the ASE library.
    It provides functionality to read molecular data, compute coordination numbers, identify terminal atoms,
    determine bonding information, and classify molecular topology into circular and residual molecules.
    """

    def __init__(self, file: Path) -> None:
        """ 
        Initializes a new instance of the Topo class with the provided file path.

        Args:
            __fileName(Path): The path to the input file containing molecular data.
            __mol(Atoms|list[Atoms]): The molecular structure(s) read from the input file.
            __neighbors(dict[int,npt.NDArray[np.int64]]): A dictionary of atom indices and their respective neighbors.
            idx1_Hydrogen_atom(list[int]): A list of indices for hydrogen atoms in the molecule.
        """
        self.__fileName: Path = Path(file)
        self.__mol: Atoms | list[Atoms]
        self.__neighbors: dict[int, npt.NDArray[np.int64]]
        self.__mol, self.__neighbors = ml4nmr.read_mol_neighbors(
            self.__fileName)
        self.idx1_Hydrogen_atom: list[int] = [idx+1 for idx,
                           i in enumerate(self.__mol) if i.symbol == "H"]  # type: ignore # nopep8

    def get_cn(self) -> dict[int, int]:
        """ 
        Computes and returns the coordination number (CN) for each atom in the molecular structure.

        Returns: 
            dict[int,int]: A dict of atom indice and coordination number for each atom.  
        """
        cn: dict[int, int] = {}
        for key, value in self.__neighbors.items():
            cn[key] = len(value)
        return cn

    def method_broken_bond_H(self, args: argparse.Namespace) -> list[int]:
        """ 
        Identifies terminal atoms involved in a broken bond, including hydrogen atoms.

        Args:
            args[argparse.Namespace]: Command-line arguments containing information about the broken bond and whether to print results.
            args.bond_broken[int,int]: atom's index of broken-bond[include 1's atom , not include 2's atom]
            args.print[bool]: print the final data on screen   

        Returns:
            list[int]: A list of atom indices involved in the broken bond (including H atoms).
        """
        Res: list[int] = self.method_broken_bond(args)
        neighbors: dict[int, npt.NDArray] = self.__neighbors
        idx1_H_atoms: list[int] = self.idx1_Hydrogen_atom
        NeighborsAtoms_H_atoms: dict[int, int] = {}  # {H:C}
        for idx in idx1_H_atoms:
            NeighborsAtoms_H_atoms[idx] = neighbors[idx][0]
        addition: list[int] = []
        for idx in Res:
            for key, value in NeighborsAtoms_H_atoms.items():
                if idx == value:
                    addition.append(key)
        Res = Res + addition
        Res.sort()
        if args.print:
            print(" Terminal_Atoms_int (include H) : ", Res)
        return Res

    def method_broken_bond(self, args: argparse.Namespace) -> list[int]:
        """
        Identifies terminal atoms involved in a broken bond, excluding hydrogen atoms.

        Args:
            args[argparse.Namespace]: Command-line arguments containing information about the broken bond and whether to print results.
            args.bond_broken[int,int] : atom's index of broken-bond [include 1's atom, not include 2's atom]
            args.print[bool] : print the final data on screen   

        Returns:
            list[int]: A list of atom indices involved in the broken bond (excluding H atoms).
        """
        idx_p, idx_q = args.bond_broken

        neighbors: dict[int, npt.NDArray] = self.__neighbors
        idx1_Hydrogen_atoms: list[int] = self.idx1_Hydrogen_atom
        idx1_Hydrogen_atoms.append(idx_q)
        NeighborsAtoms_not_H: dict[int, npt.NDArray] = {}
        for idx in neighbors.keys():
            import numpy as np
            NeighborsAtoms_not_H[idx] = np.array(
                [i for i in neighbors[idx] if int(i) not in idx1_Hydrogen_atoms])
        del NeighborsAtoms_not_H[idx_q]
        Terminal_Atoms: list[int] = [idx_p]
        Complete_Atoms: bool = False
        while (not Complete_Atoms):
            Is_Terminal_Atoms: bool = True
            for idx in Terminal_Atoms:
                for idy in NeighborsAtoms_not_H[idx]:
                    if idy in Terminal_Atoms and Is_Terminal_Atoms:
                        Is_Terminal_Atoms = True
                    else:
                        Terminal_Atoms.append(int(idy))
                        Terminal_Atoms = list(set(Terminal_Atoms))
                        Complete_Atoms, Is_Terminal_Atoms = False, False
            if Is_Terminal_Atoms:
                Complete_Atoms = True
        if args.print:
            print(" Terminal_Atoms (not H) : ", Terminal_Atoms)
        return Terminal_Atoms

    def method_bonding(self, args: argparse.Namespace) -> list[int]:
        """ 
        Retrieves the bonding partners for a specified atom, excluding hydrogen atoms.

        Args:
            args.bonding[int]: atom's index
            args.print[bool]: print the List of bonding

        Returns:
            list[int]: A list of atom indices bonded to the specified atom (excluding H atoms).
        """
        idx_p: int = args.bonding
        # ic(args.file)
        neighbors: dict[int, npt.NDArray] = self.__neighbors
        idx1_Hydrogen_atoms: list[int] = self.idx1_Hydrogen_atom
        Neighbors_Atoms: list[int] = neighbors[idx_p].tolist()
        Neighbors_Atoms = [
            i for i in Neighbors_Atoms if i not in idx1_Hydrogen_atoms]
        Neighbors_Atoms.sort()
        if args.print:
            print(" Bonding : ", idx_p, " @ ", Neighbors_Atoms)
        return Neighbors_Atoms

    def topology(self) -> tuple[ml4nmr.Atoms | list[ml4nmr.Atoms], dict[int, npt.NDArray], list[list[int]], list[list[np.int64]]]:
        """ 
        Analyzes the molecular structure to classify it into circular and residual molecules.

        Returns:
            tuple: A tuple containing the following elements:
                mol(ml4nmr.Atoms|list[ml4nmr.Atoms]): The original molecular structure(s).
                neighbors(dict[int,npt.NDArray]): Updated dictionary of atom indices and their respective neighbors (excluding H atoms).
                circle_Mols(list[list[int]]): A list of circular molecules identified in the structure.
                residual_Mols (list[list[np.int64]]): A list of residual molecules identified in the structure.
        """

        mol: ml4nmr.Atoms | list[ml4nmr.Atoms] = self.__mol
        neighbors: dict[int, npt.NDArray] = self.__neighbors.copy()
        # neighbors is removed all H-atoms
        idx1_Hydorgen_atoms: list[int] = self.idx1_Hydrogen_atom
        for key, value in neighbors.copy().items():
            if key in idx1_Hydorgen_atoms:
                del neighbors[key]
        for key, value in neighbors.copy().items():
            neighbors[key] = np.array(
                [a for a in value if a not in idx1_Hydorgen_atoms])

        # Tranfer neighbors to Graph
        graph_in: list[tuple[int, npt.NDArray]] = list()
        for key, value in neighbors.items():
            for x in value:
                graph_in.append((key, x))
        g = Graph(from_list=graph_in)

        # Get the node of bonding numbers 3 to 6
        circle_Atoms: list[int] = list()
        # for a in [3, 4, 5, 6]:
        for a in [3,]:
            circle_Atoms.extend(g.nodes(in_degree=a))  # type: ignore
        circle_Atoms.sort()

        # use Graph Theory to collect the circle sturcutres and not repeated
        circle_Mols: list[list[int]] = list()
        for atom in circle_Atoms:
            for neighbors_atoms in neighbors[atom]:
                start: int = atom
                end = neighbors_atoms
                g = Graph(from_list=graph_in)
                g.del_edge(start, end)
                # circle_Mols.append
                for x in g.all_paths(start, end):
                    the_same: bool = False
                    for circle_Mol in circle_Mols:
                        if set(x) == set(circle_Mol):
                            the_same = True
                    if not the_same:
                        circle_Mols.append(x)

        # Remove the repeated the same Atoms by use the set function (the same of the length)
        for x in circle_Mols.copy():
            if len(set(x)) != len(x):
                circle_Mols.remove(x)

        # 2D list to flatten to 1D list
        flat_circle_Mols: list[int | list[int]] = []
        for row in circle_Mols:
            flat_circle_Mols += row

        # Get residual atoms of circule molecules by use difference set
        residual_atoms: list[int] = list(
            set(neighbors.keys()).difference(set(flat_circle_Mols)))
        residual_atoms.sort()

        # g_straight is the Graph and delete the edge of every circle_Mol
        g = Graph(from_list=graph_in)
        for circle_Mol in circle_Mols:
            len_circle_Mol: int = len(circle_Mol)
            for n in range(len_circle_Mol):
                start, end = circle_Mol[n-1], circle_Mol[n]
                g.del_edge(start, end)
                g.del_edge(end, start)
        g_straight: Graph = g.copy()

        # residual_Mols is use graph : is_connected to find the connect node and append
        residual_Mols: list[list[np.int64]] = []
        for atom in residual_atoms:
            Molecules: list[np.int64] = []
            for node in list(g_straight.nodes()):  # type: ignore
                if g_straight.is_connected(atom, node):
                    Molecules.append(node)
            the_same: bool = False
            for residual_Mol in residual_Mols:
                if set(Molecules) == set(residual_Mol):
                    the_same = True
            if not the_same:
                residual_Mols.append(Molecules)
        return mol, neighbors, circle_Mols, residual_Mols
