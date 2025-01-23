#!/usr/bin/env python3
#                                   [07.11.2023] vitamin.cheng@gmail.com
# Need ml4nmr.py
#  ASE library / graph-theory library

from sys import argv as sysargv
import argparse
from turtle import circle
import Tools.ml4nmr as ml4nmr
from icecream import ic
from graph import Graph
from pathlib import Path


class topo():

    def __init__(self, fileName: Path) -> None:
        self.__fileName: Path = Path(fileName)
        self.__mol, self.__neighbors = ml4nmr.read_mol_neighbors(
            self.__fileName)
        self.idx_H_atom: list = [idx+1 for idx,
                           i in enumerate(self.__mol) if i.symbol == "H"]  # type: ignore # nopep8

    def get_cn(self) -> dict[int, int]:
        """
        Get the coordination number (CN) as number of neighbors (input is neighbors from read_mol)
        """
        # ic(self.__neighbors)
        cn = {}
        for key, value in self.__neighbors.items():
            cn[key] = len(value)
        return cn

    def Broken_bond_H(self, args: argparse.Namespace) -> list[int]:
        """_summary_
        Atom's index of your assigned bond (include 1's atom and include H atom)
        Args:
            args [argparse.Namespace()]
            args.bond_broken [int,int] : atom's index of broken-bond
                                         [include 1's atom , not include 2's atom]
            args.print [bool] : print the final data on screen   
        Returns:
            [list,int] : the list of atom's index in your assigned bond (include H atom)
        """
        Result: list[int] = self.Broken_bond(args)
        mol, neighbors = self.__mol, self.__neighbors
        idx_H_atom = self.idx_H_atom
        NeighborsAtoms_H_atom: dict[int, int] = {}  # {H:C}
        for idx in idx_H_atom:
            NeighborsAtoms_H_atom[idx] = neighbors[idx][0]
        addition: list = []
        for idx in Result:
            for key, value in NeighborsAtoms_H_atom.items():
                if idx == value:
                    addition.append(key)
        Result = Result + addition
        Result.sort()
        if args.print:
            print(" Terminal_Atoms_int (include H) : ", Result)
        return Result

    def Broken_bond(self, args: argparse.Namespace) -> list[int]:
        """_summary_
        Atom's index of your assigned bond (include 1's atom and not include H atom)
        Args:
            args [argparse.Namespace()]
            args.bond_broken [int,int] : atom's index of broken-bond
                                        [include 1's atom, not include 2's atom]
            args.print [bool] : print the final data on screen   
        Returns:
            [list,int] : the list of atom's index in your assigned bond (not include H atom)
        """
        idx_p, idx_q = args.bond_broken[0], args.bond_broken[1]
        mol, neighbors = self.__mol, self.__neighbors
        idx_H_atom: list[int] = self.idx_H_atom
        idx_H_atom.append(idx_q)
        NeighborsAtoms_not_H: dict = {}
        for idx in neighbors.keys():
            import numpy as np
            NeighborsAtoms_not_H[idx] = np.array(
                [i for i in neighbors[idx] if int(i) not in idx_H_atom])
        del NeighborsAtoms_not_H[idx_q]
        Terminal_Atoms: list[int] = [idx_p]
        Complete_Atoms: bool = False
        while (not Complete_Atoms):
            Is_Terminal_Atoms: bool = True
            for idx in Terminal_Atoms:
                for idy in NeighborsAtoms_not_H[idx]:
                    if idy in Terminal_Atoms and (Is_Terminal_Atoms == True):
                        Is_Terminal_Atoms = True
                    else:
                        Terminal_Atoms.append(int(idy))
                        Terminal_Atoms = list(set(Terminal_Atoms))
                        Complete_Atoms, Is_Terminal_Atoms = False, False
            if Is_Terminal_Atoms == True:
                Complete_Atoms = True
        if args.print:
            print(" Terminal_Atoms (not H) : ", Terminal_Atoms)
        return Terminal_Atoms

    def Bonding(self, args: argparse.Namespace) -> list[int]:
        """_summary_
            Return to list of atom's index of bonding in your speific atom
        Args:
            args [argparse.Namespace()]
            args.bonding [int] : atom's index
            args.print [bool] : print the List of bonding
            args.debug [bool] : Debug mode and show the detail of cal.   
        Returns:
            [list,int] : atom's index of bonding
        """
        from icecream import ic
        idx_p: int = args.bonding
        # ic(args.file)
        mol, neighbors = self.__mol, self.__neighbors
        idx_H_atom: list = self.idx_H_atom
        Neighbors_Atoms: list[int] = neighbors[idx_p].tolist()
        Neighbors_Atoms = [i for i in Neighbors_Atoms if i not in idx_H_atom]
        Neighbors_Atoms.sort()
        if args.print:
            print(" Bonding : ", idx_p, " @ ", Neighbors_Atoms)
        return Neighbors_Atoms

    def topology(self):
        '''
        '''
        mol, neighbors = self.__mol, self.__neighbors.copy()
        # neighbors is removed all H-atoms
        idx_H_atoms: list[int] = self.idx_H_atom
        # ic(mol)
        # ic(idx_H_atoms)
        # ic(self.__neighbors)
        # ic(neighbors)
        for key, value in neighbors.copy().items():
            if key in idx_H_atoms:
                del neighbors[key]
        import numpy as np
        for key, value in neighbors.copy().items():
            neighbors[key] = np.array(
                [a for a in value if a not in idx_H_atoms])
        # ic(idx_H_atoms)
        # ic(mol, neighbors)
        # Tranfer neighbors to Graph
        graph_in: list = list()
        for key, value in neighbors.items():
            for x in value:
                graph_in.append((key, x))
        g = Graph(from_list=graph_in)
        # ic(graph_in)
        # Get the node of bonding numbers 3 to 6
        circle_Atoms: list[int] = list()
        # for a in [3, 4, 5, 6]:
        for a in [3,]:
            # ic(g.nodes(in_degree=a))  # type: ignore
            circle_Atoms.extend(g.nodes(in_degree=a))  # type: ignore
        circle_Atoms.sort()
        # ic(circle_Atoms)
        # use Graph Theory to collect the circle sturcutres and not repeated
        circle_Mols: list[list[int]] = list()
        for atom in circle_Atoms:
            for neighbors_atoms in neighbors[atom]:
                # ic(atom, neighbors_atoms)
                start: int = atom
                end = neighbors_atoms
                g = Graph(from_list=graph_in)
                g.del_edge(start, end)
                # circle_Mols.append
                for x in g.all_paths(start, end):
                    # ic(start, end, mol)
                    the_same: bool = False
                    for circle_Mol in circle_Mols:
                        if set(x) == set(circle_Mol):
                            the_same = True
                    if the_same == False:
                        circle_Mols.append(x)
        # ic(circle_Mols)
        # Remove the repeated the same Atoms by use the set function (the same of the length)
        for x in circle_Mols.copy():
            if len(set(x)) != len(x):
                circle_Mols.remove(x)
        # ic(circle_Mols)
        # 2D list to flatten to 1D list
        flat_circle_Mols: list = []
        for row in circle_Mols:
            flat_circle_Mols += row
        # Get residual atoms of circule molecules by use difference set
        residual_atoms = list(
            set(neighbors.keys()).difference(set(flat_circle_Mols)))
        residual_atoms.sort()
        # ic(flat_circle_Mols)
        # ic(circle_Mols)
        # ic(residual_atoms)
        # g_straight is the Graph and delete the edge of every circle_Mol
        g = Graph(from_list=graph_in)
        for circle_Mol in circle_Mols:
            len_circle_Mol: int = len(circle_Mol)
            for n in range(len_circle_Mol):
                start, end = circle_Mol[n-1], circle_Mol[n]
                g.del_edge(start, end)
                g.del_edge(end, start)
        g_straight: Graph = g.copy()
        # ic(g_straight.nodes())
        # residual_Mols is use graph : is_connected to find the connect node and append
        residual_Mols: list = []
        for atom in residual_atoms:
            Mols: list = []
            for node in list(g_straight.nodes()):  # type: ignore
                if g_straight.is_connected(atom, node) == True:
                    Mols.append(node)
            # ic(Mols)
            the_same: bool = False
            for residual_Mol in residual_Mols:
                if set(Mols) == set(residual_Mol):
                    the_same = True
            if the_same == False:
                residual_Mols.append(Mols)
        # ic(residual_Mols)
        return mol, neighbors, circle_Mols, residual_Mols
