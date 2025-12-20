#!/usr/bin/env python
import numpy as np
import numpy.typing as npt
import censo_ext.Tools.ml4nmr as ml4nmr
from graph import Graph
from ase.atoms import Atoms
from censo_ext.Tools.utility import AtomID
from censo_ext.Tools.xyzfile import GeometryXYZs


class Topo():
    """
    The Topo class is designed to analyze molecular structures using graph theory and the ASE library.
    It provides functionality to read molecular data, compute coordination numbers, identify terminal atoms,
    determine bonding information, and classify molecular topology into circular and residual molecules.
    """

    def __init__(self, xyzFile: GeometryXYZs) -> None:
        """ 
        Initializes a new instance of the Topo class with the provided file path.

        Args:
            file(Path): The path to the input file containing molecular data.

        Attributes:
            __fileName (Path): The path to the input file containing molecular data.
            __mol (Atoms | list[Atoms]): The molecular structure(s) read from the input file.
            __neighbors (dict[int, npt.NDArray[np.int64]]): A dictionary of atom indices and their respective neighbors.
            idx1_Hydrogen_atom (list[int]): A list of indices for hydrogen atoms in the molecule (1-indexed).
        """

        self.__mol: Atoms | list[Atoms]
        self.__neighbors: dict[AtomID, npt.NDArray[np.int64]]
        self.__xyzFile: GeometryXYZs = xyzFile
        self.__mol, self.__neighbors = ml4nmr.read_mol_neighbors(
            self.__xyzFile)
        self._H_atomIDs: list[AtomID] = [AtomID(idx1) for idx1,
                           i in enumerate(self.__mol, 1) if i.symbol == "H"]  # type: ignore # nopep8

    def get_cn(self) -> dict[AtomID, int]:
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

        idx_cn: dict[AtomID, int] = {}
        for key, value in self.__neighbors.items():
            idx_cn[key] = len(value)
        return idx_cn

    def method_broken_bond_H(self, _bond_broken: tuple[int, int], _print: bool) -> list[AtomID]:
        """ 
        Identifies terminal atoms involved in a broken bond, including hydrogen atoms.

        Args:
            args(argparse.Namespace): Command-line arguments containing information about 
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

        Res: list[AtomID] = self.method_broken_bond(
            _bond_broken=_bond_broken, _print=_print)
        Neighbors_H_atom: dict[AtomID, AtomID] = {}  # {H:C}
        for idx in self._H_atomIDs:
            Neighbors_H_atom[idx] = AtomID(self.__neighbors[idx][0])
        addition: list[AtomID] = []
        for idx in Res:
            for key, value in Neighbors_H_atom.items():
                if idx == value:
                    addition.append(key)
        Res = Res + addition
        Res.sort()
        if _print:
            print(f" Terminal_Atoms_int (include H) : {Res}")
        return Res

    def method_broken_bond(self, _bond_broken: tuple[int, int], _print: bool) -> list[AtomID]:
        """
        Identifies terminal atoms involved in a broken bond, excluding hydrogen atoms.

        Args:
            args[argparse.Namespace]: Command-line arguments containing information about the broken bond and whether to print results.
            args.bond_broken[int,int] : atom's index of broken-bond [include 1's atom, not include 2's atom]
            args.print[bool] : print the final data on screen   

        Returns:
            list[int]: A list of atom indices involved in the broken bond (excluding H atoms).
        """
        idx1_p, idx1_q = _bond_broken
        # neighbors: dict[AtomID, npt.NDArray[np.int64]] = self.__neighbors
        H_atoms: list[AtomID] = self._H_atomIDs.copy()
        H_atoms.append(AtomID(idx1_q))
        Neighbors_not_H: dict[AtomID, npt.NDArray[np.int64]] = {}
        for idx in self.__neighbors.keys():
            Neighbors_not_H[idx] = np.array(
                [x for x in self.__neighbors[idx] if int(x) not in H_atoms])
        del Neighbors_not_H[AtomID(idx1_q)]
        Terminal_Atoms: list[AtomID] = [AtomID(idx1_p)]
        Complete_Atoms: bool = False
        while (not Complete_Atoms):
            Is_Terminal_Atoms: bool = True
            for idx in Terminal_Atoms:
                for idy in Neighbors_not_H[idx]:
                    if idy in Terminal_Atoms and Is_Terminal_Atoms:
                        Is_Terminal_Atoms = True
                    else:
                        Terminal_Atoms.append(AtomID(int(idy)))
                        Terminal_Atoms = list(set(Terminal_Atoms))
                        Complete_Atoms, Is_Terminal_Atoms = False, False
            if Is_Terminal_Atoms:
                Complete_Atoms = True
        if _print:
            print(f" Terminal_Atoms (not H) : {Terminal_Atoms}")
        return Terminal_Atoms

    def method_bonding(self, _bonding: AtomID, _print: bool) -> list[AtomID]:
        # def method_bonding(self, args: argparse.Namespace) -> list[AtomID]:
        """ 
        Retrieves the bonding partners for a specified atom, excluding hydrogen atoms.

        Args:
            args.bonding[int]: atom's index
            args.print[bool]: print the List of bonding

        Returns:
            list[int]: A list of atom indices bonded to the specified atom (excluding H atoms).
        """
        idx1_p: AtomID = _bonding
        _Bonding_AtomIDs: list[AtomID] = self.__neighbors[idx1_p].tolist()
        _Bonding_AtomIDs = [
            x for x in _Bonding_AtomIDs if x not in self._H_atomIDs]
        _Bonding_AtomIDs.sort()
        if _print:
            print(f" Bonding : {idx1_p} @ Neighbors_Atoms")
        return _Bonding_AtomIDs

    def topology(self) -> tuple[dict[AtomID, npt.NDArray[np.int64]], list[list[AtomID]], list[set[int]], dict]:
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

        neighbors: dict[AtomID, npt.NDArray[np.int64]
                        ] = self.__neighbors.copy()

        # neighbors is removed all H-atoms
        for key, value in neighbors.copy().items():
            if key in self._H_atomIDs:
                del neighbors[key]
        for key, value in neighbors.copy().items():
            neighbors[key] = np.array(
                [x for x in value if x not in self._H_atomIDs])

        # Transfer neighbors to Graph
        graph_in: list[tuple[AtomID, AtomID]] = list()
        for key, value in neighbors.items():
            for x in value:
                graph_in.append((key, AtomID(int(x))))
        g = Graph(from_list=graph_in)

        # Get the node of bonding numbers 3 to 6
        circle_AtomIDs: list[AtomID] = list()
        # for a in [3, 4, 5, 6]:
        for a in [3,]:
            circle_AtomIDs.extend(g.nodes(in_degree=a))  # type: ignore
        circle_AtomIDs.sort()

        # use Graph Theory to collect the circle sturcutres and not repeated
        circle_Mols: list[list[AtomID]] = list()
        for atom in circle_AtomIDs:
            for neighbors_atoms in neighbors[atom]:
                start: int = atom
                end = neighbors_atoms
                g = Graph(from_list=graph_in)
                g.del_edge(start, end)
                for x in g.all_paths(start, end):
                    the_same: bool = False
                    for circle_Mol in circle_Mols:
                        if set(x) == set(circle_Mol):
                            the_same = True
                    if not the_same:
                        circle_Mols.append([AtomID(int(a)) for a in x])

        # Remove the repeated the same Atoms by use the set function (the same of the length)
        for x in circle_Mols.copy():
            if len(set(x)) != len(x):
                circle_Mols.remove(x)

        # 2D list to flatten to 1D list
        flat_circle_Mols: list[int | list[AtomID]] = []
        for row in circle_Mols:
            flat_circle_Mols += row

        # Get residual atoms of circule molecules by use difference set
        residual_atoms: list[int] = list(
            set(neighbors.keys()).difference(set(flat_circle_Mols)))
        residual_atoms.sort()

        # g_straight is the Graph and delete the edge of every circle_Mol
        g = Graph(from_list=graph_in)
        for circle_Mol in circle_Mols:
            for n in range(len(circle_Mol)):
                start, end = circle_Mol[n-1], circle_Mol[n]
                g.del_edge(start, end)
                g.del_edge(end, start)
        g_straight: Graph = g.copy()

        # residual_Mols is use graph : is_connected to find the connect node and append
        g_components = g_straight.components()
        residual_Mols: list[set[int]] = []
        for g_component in g_components:
            if len(g_component) != 1:
                residual_Mols.append({int(a) for a in g_component})
        residual_Mols_all_pairs: dict = g_straight.all_pairs_shortest_paths()

        return neighbors, circle_Mols, residual_Mols, residual_Mols_all_pairs

    def topology_components(self) -> list[set[int]]:
        '''
        Get the index of every components, but not include one independence atom.
        '''
        # Transfer neighbors to Graph
        graph_in: list[tuple[AtomID, AtomID]] = list()
        for key, value in self.__neighbors.items():
            for x in value:
                graph_in.append((key, AtomID(int(x))))
        g = Graph(from_list=graph_in)
        return g.components()
