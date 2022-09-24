### peptide.py

"""
Provides methods to obtain peptide descriptors.

    1. cyclic or linear
    2. amino acid sequence
    3. frequency of canonicals vs non-canonicals in the sequence
"""

from itertools import combinations
from typing import Dict, List

import networkx as nx
import numpy as np
from Bio.PDB.Chain import Chain
from Bio.PDB.QCPSuperimposer import QCPSuperimposer

from .utils.aas import get_codon_table
from .utils.structures import PeptideGraph

aa_code: Dict[str, str] = get_codon_table()


class Peptide(object):
    """
    Peptide (class)
    ---------------
    Provides methods to model peptides and methods
    to obtain physico-chemical properties.

    params:
    chain (Chain): Chain than references peptide sequence
    residues (List): List of residues in the chain
    graph (nx.Graph): Undirected graph

    """

    def __init__(self, chain: Chain):
        self.chain: Chain = chain
        self.residues = [
            x
            for x in self.chain.get_residues()
            if x.get_resname() != "NH2" and x.get_resname() != "HOH"
        ]
        self.graph: nx.Graph = None

    def get_sequence(self) -> str:
        """returns the 3 letter code sequence"""
        return self._get_sequence()

    def get_one_letter_sequence(self) -> str:
        """
        returns the one-letter code sequence of a peptide.
        non_canonicals are replaced with X
        """
        tmp: List[str] = self._get_sequence().split("-")
        one_seq: str = ""

        for t in tmp:
            aa = aa_code[t] if t in aa_code else "X"
            one_seq = "".join([one_seq, aa])

        return one_seq

    def _get_sequence(self) -> str:
        """Extracts the amino acid sequence from the protein chain"""
        return "-".join(x.get_resname().upper() for x in self.residues)

    def generate_peptide_graph(self) -> None:
        peptide_graph: PeptideGraph = PeptideGraph(nx.Graph(), self.residues)
        peptide_graph.build_graph()
        self.graph = peptide_graph.get_graph()

    def get_graph(self):
        return self.graph

    def get_rsmd(self, reference_structure: np.ndarray, backbone: bool = True) -> float:

        try:
            if backbone:  # Align atoms based on Ca (backbone) atoms
                ca_coords = np.ndarray(
                    [
                        atom.get_coord()
                        for atom in self.residues
                        if atom.get_name() == "CA"
                    ]
                )

                assert ca_coords.shape == reference_structure.shape
                qcp_superimposer: QCPSuperimposer = QCPSuperimposer()
                qcp_superimposer.set(reference_structure, ca_coords)
                qcp_superimposer.run()

                return qcp_superimposer.rms
            else:
                cb_coords = np.ndarray(
                    [
                        atom.get_coord()
                        for atom in self.residues
                        if atom.get_name() == "CB"
                    ]
                )

                assert ca_coords.shape == reference_structure.shape
                qcp_superimposer: QCPSuperimposer = QCPSuperimposer()
                qcp_superimposer.set(reference_structure, cb_coords)
                qcp_superimposer.run()

                return qcp_superimposer.rms

        except Exception as e:
            raise e(f"{e.__repr__()}: Coordinates are of different dimensions")

    def cyclization_type(self) -> str:
        pass

    def _disulfide_cycle(self) -> bool:
        cycle_members: List = self.get_cycle_members()
        cysteines: List = [
            member for member in cycle_members if member.get_resname().upper() == "CYS"
        ]

        if len(cysteines) < 2:
            return False

        for res_one, res_two in combinations(cysteines, 2):
            sulphur_coords_one: List = [
                atom for atom in res_one.get_atoms() if atom.get_name() == "S"
            ][0]
            sulphur_coords_two: List = [
                atom for atom in res_two.get_atoms() if atom.get_name() == "S"
            ][0]

            distance: float = np.linalg.norm(sulphur_coords_one - sulphur_coords_two)

            if distance <= self.disulfide_length:
                return True

        return False

    def _n_to_c_cycle(self) -> bool:
        pass
