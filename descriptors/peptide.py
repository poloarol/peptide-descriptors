### peptide.py

"""
Provides methods to obtain peptide descriptors.

    1. cyclic or linear
    2. amino acid sequence
    3. frequency of canonicals vs non-canonicals in the sequence
"""

from typing import List, Dict
from itertools import combinations

import numpy as np
import networkx as nx

from Bio.PDB import PPBuilder
from Bio.PDB.Chain import Chain
from Bio.SeqUtils.ProtParam import ProteinAnalysis

aa_code: Dict[str, str] = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
                            'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
                            'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
                            'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

class Peptide(object):
    """
    class Peptide
        
    Args:
        
    """
    
    def __init__(self, chain: Chain, pep_bond: float = 3.0):
        self.chain = chain
        self.bond_length = pep_bond
        self.graph = nx.Graph()
        
    def __post_init__(self):
        self.residues: List = list(self.chain.get_residues())
        self._build_graph()
        
    def get_sequence(self) -> str:
        """ returns the 3 letter code sequence """
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
            one_seq = "".join(one_seq, aa)
            
        return one_seq
    
    def _get_sequence(self) -> str:
        """ Extracts the amino acid sequence from the protein chain """
        return "-".join(x.get_resname().upper() for x in self.residues)
    
    def _get_num_canonicals(self) -> int:
        sequence: str = self.get_one_letter_sequence()
        num_non_canonicals: int = self._get_num_non_canonicals()
        return len(sequence) - num_non_canonicals
    
    def _get_num_non_canonicals(self) -> int:
        sequence: str = self.get_one_letter_sequence()
        return sequence.count("X")
    
    def get_num_canonical_non_canonical(self) -> Dict[str, int]:
        sequence: str = self.get_one_letter_sequence()
        
        return {"canonical": self._get_num_canonicals(), "non_canonical": self._get_num_non_canonicals()}
    
    def _get_polypeptide(self) -> PPBuilder:
        return PPBuilder(self.chain)

    def get_dihedral_angles(self) -> List:
        
        dihedrals: List = []
        polypeptide: PPBuilder = self._get_polypeptide()
        
        for _, peptide in enumerate(polypeptide):
            dihedrals.append(peptide.get_phi_psi_list())
            
        return dihedrals
    
    def get_tau_angles(self) -> List:
        """ List of tau torsions angles for all 4 consecutive Calpha atoms """
        
        tau_angles: List = []
        polypeptide: PPBuilder = self._get_polypeptide()
        
        for _, peptide in enumerate(polypeptide):
            tau_angles.append(peptide.get_tau_angles())
        
        return tau_angles
    
    def get_theta_angles(self) -> List:
        """ List of theta angles for all 3 consecutive Calpha atoms """
        
        theta_angles: List = []
        polypeptide: PPBuilder = self._get_polypeptide()
        
        for _, peptide in enumerate(polypeptide):
            theta_angles.append(peptide.get_theta_list())
        
        return theta_angles
           
    def _build_graph(self) -> None:
        for res_one, res_two in combinations(self.residues, 2):
            ca_coords_one: np.ndarray = [atom for atom in res_one.get_atoms() if atom.get_name() == 'CA'][0]
            ca_coords_two: np.ndarray = [atom for atom in res_two.get_atoms() if atom.get_name() == 'CA'][0]
            
            distance: float = np.linalg.norm(ca_coords_one - ca_coords_two)
            
            if distance <= self.bond_length:
                self.graph.add_nodes([res_one, res_two])
                self.graph.add_edges_from([(res_one, res_two), (res_two, res_one)])
                
    def get_graph(self) -> nx.Graph:
        return self.graph
                   
    def is_cyclic(self, source: str = None, orientation: str = None) -> bool:
        """ Determines if a peptide is cyclic or not"""
        
        try:
            cycle = nx.find_cycle(self.graph, source=source, orientation=orientation)
            return True if cycle else False
        except Exception as e:
            return False
    
    def get_cycle_members(self, source: str = None, orientation: str = None) -> List:
        """ Provides all residues in the cycle"""
        try:
            cycle = nx.find_cycle(self.graph, source=source, orientation=orientation)
            return list(cycle)
        except nx.exception.NetworkXNoCycle as e:
            raise f"{e.__repr__()}: No Cycle found"

