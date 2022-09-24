### structures.py
### Author Paul A. Wambo

### Describes the data structures used to build the peptide graph

from collections import namedtuple
from itertools import combinations, product
from typing import List, NamedTuple, Set

import networkx as nx
import numpy as np
from Bio.PDB.Residue import Residue
from scipy.spatial import distance

from .aas import get_acidic_aas, get_all_aas, get_basic_aas, get_neutral_aas

all_side_chains: Set[str] = get_all_aas()
neutral_side_chains: Set[str] = get_neutral_aas()
basic_aas: Set[str] = get_basic_aas()
acidic_aas: Set[str] = get_acidic_aas()

BOND_LENGTHS: NamedTuple = namedtuple("BondLengths", "min max")


class ResidueNode(object):
    """
    Residue (class)
    ---------------
    Provides accessory methods to allow for the
    addition of bonds between amino acid residues.

    params
    ------
    name(str): Name of the residue
    residue(Residue): Residue descriptor as provided by Biopython
    num_peptide_bonds (int): Number of peptide bonds formed
    num_disulfide_bonds (int): Number of disulfide bonds formed
    disulfide (bool): Determines if disulfide bonds can be formed
    acidic_side_chain (bool): Determines if the side chain is acidic
    basic_side_chain (bool): Determines if the side chain is basic

    """

    def __init__(self, name: str, residue: Residue):
        self.name: str = name
        self.residue: Residue = residue
        self.num_peptide_bonds: int = 0
        self.num_disulfide_bonds: int = 0
        self.MAX_PEPTIDE_BONDS = 2
        self.MAX_DISULFIDE_BONDS = 1

        self.disulfide: bool = self._has_thiol_side_chain()
        self.acidic_side_chain: bool = self._has_acidic_side_chain()
        self.basic_side_chain: bool = self._has_basic_side_chain()

    def get_coordinates(self, element_type="SG") -> List:
        if element_type == "carbon":
            return [
                np.array(atom.get_coord())
                for atom in self.residue
                if atom.get_name().startswith("C")
            ]
        elif element_type == "alpha-carbon":
            return [
                np.array(atom.get_coord())
                for atom in self.residue
                if atom.get_name() == "CA" or atom.get_name() == "C1"
            ]
        else:  # All carbon types. More useful for non-canonicals or amino acids which can form peptide bonds with their side chain.
            return [
                np.array(atom.get_coord())
                for atom in self.residue
                if atom.get_name() == "SG"
            ]

    def _get_side_chain_atoms(self):
        backbone_atoms: List = (
            "CA",
            "C1",
            "HA",
            "N",
            "HN",
            "H",
            "C",
            "O",
            "H1",
            "H2",
            "H3",
            "NH1",
            "NH2",
            "NH3",
            "OXT",
        )
        side_chain_atoms: List[str] = [
            atom.get_name()
            for atom in self.residue.get_atoms()
            if atom.get_name() not in (backbone_atoms)
        ]
        return side_chain_atoms

    def _detect_moeity(self, moeity: str):
        side_chain_atoms: List[str] = self._get_side_chain_atoms()
        return (
            True
            if [element for element in side_chain_atoms if element.startswith(moeity)]
            else False
        )

    def _has_thiol_side_chain(self) -> bool:
        return (
            True
            if [element for element in self.residue if element.get_name() == "SG"]
            else False
        )

    def is_acidic(self) -> bool:
        return self.acidic_side_chain

    def is_basic(self) -> bool:
        return self.basic_side_chain

    def can_form_disulfide_bond(self) -> bool:
        return (
            True
            if (self.num_disulfide_bonds < self.MAX_DISULFIDE_BONDS and self.disulfide)
            else False
        )

    def can_form_peptide_bond(self) -> bool:
        return True if self.num_peptide_bonds < self.MAX_PEPTIDE_BONDS else False

    def is_neutral(self):
        if self.residue.get_resname().upper() in neutral_side_chains:
            return True
        else:
            side_chain_atoms: Set[str] = self._get_side_chain_atoms()
            for atom in side_chain_atoms:
                if "C" in atom or "H" in atom:
                    pass
                else:
                    return False
            return True

    def _has_basic_side_chain(self) -> bool:
        if self.residue.get_resname().upper() in basic_aas:
            return True
        elif self.residue.get_resname().upper() in neutral_side_chains:
            return False
        else:
            if self._detect_moeity(moeity="N"):
                return True
        return False

    def _has_acidic_side_chain(self) -> bool:
        if self.residue.get_resname().upper() in acidic_aas:
            return True
        elif self.residue.get_resname().upper() in neutral_side_chains:
            return False
        else:
            if self._detect_moeity(moeity="O"):
                return True
        return False

    def _increment_num_peptide_bonds(self) -> None:
        self.num_peptide_bonds = self.num_peptide_bonds + 1

    def _increment_num_disulfide_bonds(self) -> None:
        self.num_disulfide_bonds = self.num_disulfide_bonds + 1

    def form_peptide_bond(self) -> None:
        if self.num_peptide_bonds < self.MAX_PEPTIDE_BONDS:
            self._increment_num_peptide_bonds()

    def form_disulfide_bond(self) -> None:
        if self.disulfide and self.num_disulfide_bonds < self.MAX_DISULFIDE_BONDS:
            self._increment_num_disulfide_bonds()


class PeptideGraph(object):  # Check if inheritance from nx.Graph would be better
    """ """

    def __init__(self, graph: nx.Graph, residues: List[Residue]):
        self.graph: nx.Graph = graph
        self.residues: List[Residue] = residues
        self.peptide_length: NamedTuple = BOND_LENGTHS(
            2.899, 4.18
        )  # bond lengths should be ammended
        self.disulfide_length: NamedTuple = BOND_LENGTHS(1.999, 2.1)
        self.ionic_radii_length: NamedTuple = BOND_LENGTHS(0.3, 30)

    def build_graph(self) -> None:
        # It would be better to obtain statistics about bond lenghts,
        # such that we do not include bonds that are too short to be realistic
        # or exist

        # I need to think about how to identify the places a cycle can occur, rather than
        # blindly adding them
        # 1. 1st and last amino acid
        # 2. 1st/last and some amino acid in the middle
        # 3. Two amino acid in the middle (Apart of CYS, two similar canonical amino acids can only link together through a peptide bond)

        try:
            for res_one, res_two in combinations(self.residues, 2):
                residue_one_node: ResidueNode = None
                residue_two_node: ResidueNode = None

                name_one, name_two = (
                    res_one.get_resname().upper(),
                    res_two.get_resname().upper(),
                )
                name_one, name_two = "".join(name_one.split()), "".join(
                    name_two.split()
                )

                first_residue_idx, last_residue_idx = (
                    self.residues[0].id[1],
                    self.residues[-1].id[1],
                )
                res_one_in_graph: bool = self._find_node(
                    attribute="tag", value=f"{name_one}_{res_one.id[1]}"
                )
                res_two_in_graph: bool = self._find_node(
                    attribute="tag", value=f"{name_two}_{res_two.id[1]}"
                )

                if not res_one_in_graph and not res_two_in_graph:
                    residue_one_node = ResidueNode(name=name_one, residue=res_one)
                    residue_two_node = ResidueNode(name=name_two, residue=res_two)
                elif not res_one_in_graph:
                    residue_one_node = self._get_node(
                        attribute="tag", value=f"{name_one}_{res_one.id[1]}"
                    )
                    residue_two_node = ResidueNode(name=name_one, residue=res_one)
                elif not res_two_in_graph:
                    residue_two_node = self._get_node(
                        attribute="tag", value=f"{name_two}_{res_two.id[1]}"
                    )
                    residue_one_node = ResidueNode(name=name_two, residue=res_two)
                else:
                    residue_one_node = self._get_node(
                        attribute="tag", value=f"{name_one}_{res_one.id[1]}"
                    )
                    residue_two_node = self._get_node(
                        attribute="tag", value=f"{name_two}_{res_two.id[1]}"
                    )

                # # MET does not form disulfide bonds
                if (
                    residue_one_node.can_form_disulfide_bond()
                    and residue_two_node.can_form_disulfide_bond()
                ):
                    self._add_disulfide_edge(residue_one_node, residue_two_node)

                if (
                    res_one.id[1] == first_residue_idx
                    or res_two.id[1] == first_residue_idx
                ) and (
                    res_one.id[1] == last_residue_idx
                    or res_two.id[1] == last_residue_idx
                ):  # First and last residues
                    self._add_peptide_edge(residue_one_node, residue_two_node)

                elif (
                    (res_one.id[1] == first_residue_idx)
                    and (res_two.id[1] != last_residue_idx)
                    and (res_one.id[1] not in (res_two.id[1] - 1, res_one.id[1] + 1))
                ):
                    self._add_peptide_edge(residue_one_node, residue_two_node)

                elif (
                    (res_two.id[1] == first_residue_idx)
                    and (res_one.id[1] != last_residue_idx)
                    and (res_one.id[1] not in (res_two.id[1] - 1, res_one.id[1] + 1))
                ):
                    self._add_peptide_edge(residue_one_node, residue_two_node)

                elif (
                    (res_one.id[1] == last_residue_idx)
                    and (res_two.id[1] != first_residue_idx)
                    and (res_one.id[1] not in (res_two.id[1] - 1, res_one.id[1] + 1))
                ):
                    self._add_peptide_edge(residue_one_node, residue_two_node)

                elif (
                    (res_two.id[1] == last_residue_idx)
                    and (res_one.id[1] != first_residue_idx)
                    and (res_one.id[1] not in (res_two.id[1] - 1, res_one.id[1] + 1))
                ):
                    self._add_peptide_edge(residue_one_node, residue_two_node)

                elif res_one.id[1] in (
                    res_two.id[1] - 1,
                    res_two.id[1] + 1,
                ) or res_two.id[1] in (res_one.id[1] - 1, res_one.id[1] + 1):
                    # sequential amino acids
                    self._add_peptide_edge(residue_one_node, residue_two_node)
                else:
                    if (name_one != name_two) and not (
                        residue_one_node.is_neutral() or residue_two_node.is_neutral()
                    ):
                        # make sure the amino acid across is different and is neutral
                        self._add_peptide_edge(residue_one_node, residue_two_node)
        except IndexError as e:
            return f"Error message: {e.__repr__()}"

    def get_graph(self) -> nx.Graph:
        return self.graph

    def _find_node(self, attribute: str, value: str) -> bool:
        try:
            return any(
                [
                    node
                    for node in self.graph.nodes(data=True)
                    if node[1][attribute] == value
                ]
            )
        except KeyError as e:
            print(f"Key not found: {e.__repr__()}")
        finally:
            return False

    def _get_node(self, attribute: str, value: str) -> ResidueNode:
        try:
            return [
                node
                for node in self.graph.nodes(data=True)
                if node[1][attribute] == value
            ][0]
        except IndexError as e:
            print(f"Key not found: {e.__repr__()}")
        # Find appropriate way to handle the finally clause here

    def _add_peptide_edge(
        self, residue_one: ResidueNode, residue_two: ResidueNode
    ) -> None:
        if residue_one.can_form_peptide_bond() and residue_two.can_form_peptide_bond():

            coordinates_one: List
            coordinates_two: List

            add_node: bool = False

            # if either of residue_one or residue_two cannot form peptide bonds with their side chain,
            # we only get the Ca carbon else, we check all other carbons
            if (
                residue_one.name in all_side_chains
                or residue_two.name in all_side_chains
            ):
                coordinates_one: List = residue_one.get_coordinates(
                    element_type="alpha-carbon"
                )
                coordinates_two: List = residue_two.get_coordinates(
                    element_type="alpha-carbon"
                )
            else:
                coordinates_one: List = residue_one.get_coordinates(
                    element_type="carbon"
                )
                coordinates_two: List = residue_two.get_coordinates(
                    element_type="carbon"
                )

            for carbon_one, carbon_two in product(coordinates_one, coordinates_two):
                distance = np.linalg.norm(carbon_one - carbon_two)
                bond_angle = self._theta(carbon_one, carbon_two)

                if (
                    distance >= self.peptide_length.min
                    and distance <= self.peptide_length.max
                ) and ((-0.1 <= bond_angle <= 0.1) or (179.9 <= bond_angle <= 180.9)):
                    add_node = True
                    break

            if add_node:
                residue_one_node: ResidueNode = None
                residue_two_node: ResidueNode = None

                node_one_present: bool = self._find_node(
                    attribute="tag",
                    value=f"{residue_one.name}_{residue_one.residue.id[1]}",
                )
                node_two_present: bool = self._find_node(
                    attribute="tag",
                    value=f"{residue_two.name}_{residue_two.residue.id[1]}",
                )

                if not node_one_present:
                    self._add_node(residue_one)

                if not node_two_present:
                    self._add_node(residue_two)

                residue_one_node = self._get_node(
                    attribute="tag",
                    value=f"{residue_one.name}_{residue_one.residue.id[1]}",
                )
                residue_two_node = self._get_node(
                    attribute="tag",
                    value=f"{residue_two.name}_{residue_two.residue.id[1]}",
                )

                if residue_one_node and residue_two_node:
                    residue_one.form_peptide_bond()
                    residue_two.form_peptide_bond()
                    self._add_edge(residue_one_node[0], residue_two_node[0])

    def _add_disulfide_edge(
        self, residue_one: ResidueNode, residue_two: ResidueNode
    ) -> None:
        coordinates_one: List = residue_one.get_coordinates()[0]
        coordinates_two: List = residue_two.get_coordinates()[0]

        atomic_distance: float = distance.euclidean(
            coordinates_one, coordinates_two
        )  # For some reason, numy.linalg.norm does not work! Interesting
        bond_angle: float = self._theta(coordinates_one, coordinates_two)

        if (
            self.disulfide_length.min <= atomic_distance <= self.disulfide_length.max
        ) and ((-0.1 <= bond_angle <= 0.1) or (179.9 <= bond_angle <= 180.9)):

            residue_one_node: ResidueNode = None
            residue_two_node: ResidueNode = None

            node_one_present: bool = self._find_node(
                attribute="tag", value=f"{residue_one.name}_{residue_one.residue.id[1]}"
            )
            node_two_present: bool = self._find_node(
                attribute="tag", value=f"{residue_two.name}_{residue_two.residue.id[1]}"
            )

            if not node_one_present:
                self._add_node(residue_one)

            if not node_two_present:
                self._add_node(residue_two)

            residue_one_node = self._get_node(
                attribute="tag", value=f"{residue_one.name}_{residue_one.residue.id[1]}"
            )
            residue_two_node = self._get_node(
                attribute="tag", value=f"{residue_two.name}_{residue_two.residue.id[1]}"
            )

            if residue_one_node and residue_two_node:
                residue_one.form_disulfide_bond()
                residue_two.form_disulfide_bond()
                self._add_edge(residue_one_node[0], residue_two_node[0])

    def _add_edge(self, residue_one, residue_two, bond_type: str = "peptide") -> None:
        if bond_type == "peptide":
            self.graph.add_edge(residue_one, residue_two, length=4.0)

        if bond_type == "disulfide":
            self.graph.add_edge(residue_one, residue_two, length=2.0)

    def _add_node(self, node):
        self.graph.add_node(
            f"{node.name}-{node.residue.id[1]}",
            name=str(node.residue.id[1]),
            type="node",
            tag=f"{node.name}_{node.residue.id[1]}",
        )

    def _theta(self, v: np.array, w: np.array) -> float:
        return np.arccos(v.dot(w) / (np.linalg.norm(v) * np.linalg.norm(w)))

    def _add_ionic_interaction(self, residue_one, residue_two) -> bool:
        # calculate ionic radius.
        # make sure it falls between range
        pass

    def is_cyclic(self, source: str = None, orientation: str = None) -> bool:
        """Determines if a peptide is cyclic or not"""

        try:
            cycle = nx.find_cycle(self.graph, source=source, orientation=orientation)
            return True if cycle else False
        except Exception as e:
            return False

    def get_cycle_members(self, source: str = None, orientation: str = None) -> List:
        """Provides all residues in the cycle"""
        try:
            cycle = nx.find_cycle(self.graph, source=source, orientation=orientation)
            return list(cycle)
        except nx.exception.NetworkXNoCycle as e:
            return f"{e.__repr__()}: No Cycle found"
