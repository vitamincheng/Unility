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
            file (Path): The path to the input file containing molecular data.

        Attributes:
            __fileName (Path): The path to the input file containing molecular data.
            __mol (Atoms | list[Atoms]): The molecular structure(s) read from the input file.
            __neighbors (dict[int, npt.NDArray[np.int64]]): A dictionary of atom indices and their respective neighbors.
            idx1_Hydrogen_atom (list[int]): A list of indices for hydrogen atoms in the molecule (1-indexed).
        """

        self.__fileName: Path = Path(file)
        self.__mol: Atoms | list[Atoms]
        self.__neighbors: dict[int, npt.NDArray[np.int64]]
        self.__mol, self.__neighbors = ml4nmr.read_mol_neighbors(
            self.__fileName)
        self.idx1_Hydrogen_atom: list[int] = [idx+1 for idx,
                           i in enumerate(self.__mol) if i.symbol == "H"]  # type: ignore # nopep8

    def get_cn(self) -> dict[int, int]:
        """Computes and returns the coordination number (CN) for each atom in the molecular structure.

        The coordination number represents the number of nearest neighbors for each atom.

        Returns:
            dict[int,int]: A dictionary mapping atom indices to their respective coordination numbers.
                - Key: Atom index (int)
                - Value: Coordination number (int), representing the count of nearest neighbors

        Example:
            >>> topology.get_cn()
            {0: 4, 1: 2, 2: 4, 3: 2}

        Note:
            This method relies on the internal `self.__neighbors` attribute which should be 
            pre-computed and contain the neighbor information for each atom.
        """

        idx_cn: dict[int, int] = {}
        for key, value in self.__neighbors.items():
            idx_cn[key] = len(value)
        return idx_cn

    def method_broken_bond_H(self, args: argparse.Namespace) -> list[int]:
        """ 
        Identifies terminal atoms involved in a broken bond, including hydrogen atoms.

        Args:
            args (argparse.Namespace): Command-line arguments containing information about 
                the broken bond and whether to print results.
            args.bond_broken (tuple[int, int]): atom's index of broken-bond[include 1's atom , not include 2's atom]
            args.print (bool): print the final data on screen   

        Returns:
            list[int]: A list of atom indices involved in the broken bond (including H atoms).

        Example:
            >>> args = argparse.Namespace(bond_broken=(5, 10), print=True)
            >>> result = method_broken_bond_H(args)
            >>> print(result)
            [5, 6, 7, 10, 11]
        """

        Res: list[int] = self.method_broken_bond(args)
        idx_neighbors: dict[int, npt.NDArray] = self.__neighbors
        idx1_H_atoms: list[int] = self.idx1_Hydrogen_atom
        NeighborsAtoms_H_atoms: dict[int, int] = {}  # {H:C}
        for idx in idx1_H_atoms:
            NeighborsAtoms_H_atoms[idx] = idx_neighbors[idx][0]
        addition: list[int] = []
        for idx in Res:
            for key, value in NeighborsAtoms_H_atoms.items():
                if idx == value:
                    addition.append(key)
        Res = Res + addition
        Res.sort()
        if args.print:
            print(f" Terminal_Atoms_int (include H) : {Res}")
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

        idx_neighbors: dict[int, npt.NDArray] = self.__neighbors
        idx1_Hydrogen_atoms: list[int] = self.idx1_Hydrogen_atom
        idx1_Hydrogen_atoms.append(idx_q)
        NeighborsAtoms_not_H: dict[int, npt.NDArray] = {}
        for idx in idx_neighbors.keys():
            import numpy as np
            NeighborsAtoms_not_H[idx] = np.array(
                [i for i in idx_neighbors[idx] if int(i) not in idx1_Hydrogen_atoms])
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
            print(f" Terminal_Atoms (not H) : {Terminal_Atoms}")
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
        idx_neighbors: dict[int, npt.NDArray] = self.__neighbors
        idx1_Hydrogen_atoms: list[int] = self.idx1_Hydrogen_atom
        Neighbors_Atoms: list[int] = idx_neighbors[idx_p].tolist()
        Neighbors_Atoms = [
            i for i in Neighbors_Atoms if i not in idx1_Hydrogen_atoms]
        Neighbors_Atoms.sort()
        if args.print:
            print(f" Bonding : {idx_p} @ Neighbors_Atoms")
        return Neighbors_Atoms

    def topology(self) -> tuple[ml4nmr.Atoms | list[ml4nmr.Atoms], dict[int, npt.NDArray], list[list[int]], list[list[np.int64]]]:
        """Analyzes the molecular structure to classify it into circular and residual molecules.

        This method identifies circular (ring) structures and residual (non-ring) fragments
        within a molecular system. It processes the molecular graph by removing hydrogen atoms
        and then applies graph theory algorithms to detect ring systems and connected components
        of non-ring atoms.

        The algorithm:
        1. Removes hydrogen atoms from neighbor lists
        2. Constructs a molecular graph from remaining bonds
        3. Identifies potential ring atoms (degree 3-6)
        4. Finds complete ring structures by removing edges and finding paths
        5. Classifies remaining atoms as residual fragments

        Returns:
            tuple: A tuple containing:
                mol (ml4nmr.Atoms | list[ml4nmr.Atoms]): The original molecular structure(s).
                neighbors (dict[int, npt.NDArray]): Updated dictionary of atom indices and their respective neighbors (excluding H atoms).
                circle_Mols (list[list[int]]): A list of circular molecules identified in the structure.
                residual_Mols (list[list[np.int64]]): A list of residual molecules identified in the structure.

        Note:
            - Ring detection is performed by identifying atoms with degrees 3-6
            - Duplicate ring structures are automatically removed
            - Residual fragments are connected components of non-ring atoms
        """

        mol: ml4nmr.Atoms | list[ml4nmr.Atoms] = self.__mol
        idx_neighbors: dict[int, npt.NDArray] = self.__neighbors.copy()
        # neighbors is removed all H-atoms
        idx1_Hydorgen_atoms: list[int] = self.idx1_Hydrogen_atom
        for key, value in idx_neighbors.copy().items():
            if key in idx1_Hydorgen_atoms:
                del idx_neighbors[key]
        for key, value in idx_neighbors.copy().items():
            idx_neighbors[key] = np.array(
                [a for a in value if a not in idx1_Hydorgen_atoms])

        # Tranfer neighbors to Graph
        graph_in: list[tuple[int, npt.NDArray]] = list()
        for key, value in idx_neighbors.items():
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
            for neighbors_atoms in idx_neighbors[atom]:
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
            set(idx_neighbors.keys()).difference(set(flat_circle_Mols)))
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
        return mol, idx_neighbors, circle_Mols, residual_Mols
