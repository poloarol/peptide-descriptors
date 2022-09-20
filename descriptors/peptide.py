### peptide.py

"""
Provides methods to obtain peptide descriptors.

    1. cyclic or linear
    2. amino acid sequence
    3. frequency of canonicals vs non-canonicals in the sequence
"""

from typing import List, Dict, NamedTuple
from collections import namedtuple
from itertools import combinations, product

import numpy as np
import networkx as nx

from peptides import Peptide as PTDS
from Bio.PDB import PPBuilder
from Bio.PDB.Chain import Chain
from Bio.PDB.QCPSuperimposer import QCPSuperimposer

aa_code: Dict[str, str] = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
                            'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
                            'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
                            'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}


BOND_LENGTHS: NamedTuple = namedtuple('BondLengths', 'min max')

class Peptide(object):
    """
    class Peptide
        
    Args:
        
    """
    
    def __init__(self, chain: Chain):
        self.chain: Chain = chain
        self.peptide_length: NamedTuple = BOND_LENGTHS(2.9, 3.799) # bond lengths should be ammended
        self.disulfide_length: NamedTuple = BOND_LENGTHS(1.9, 2.05)
        self.graph = nx.Graph()
        self.residues = [x for x in self.chain.get_residues() \
                    if x.get_resname() != "NH2" and x.get_resname() != "HOH"]
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
            one_seq = "".join([one_seq, aa])
            
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
        return PPBuilder().build_peptides(self.chain)

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
            tau_angles.append(peptide.get_tau_list())
        
        return tau_angles
    
    def get_theta_angles(self) -> List:
        """ List of theta angles for all 3 consecutive Calpha atoms """
        
        theta_angles: List = []
        polypeptide: PPBuilder = self._get_polypeptide()
        
        for _, peptide in enumerate(polypeptide):
            theta_angles.append(peptide.get_theta_list())
        
        return theta_angles
           
    def _build_graph(self) -> None:
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
                
                sulphur_coords_one: List = [atom for atom in res_one.get_atoms() if atom.get_name() == 'SG']
                sulphur_coords_two: List = [atom for atom in res_two.get_atoms() if atom.get_name() == 'SG']
                
                if (res_one.id[1] == 1 or res_two.id[1] == 1) and (res_one.id[1] == len(self.residues) or res_two.id[1] == len(self.residues)): # first and last amino acid
                    self._add_edges(res_one, res_two)
                    
                elif res_one.id[1] == 1 and res_two.id[1] != len(self.residues): # first amino acid and any other amino acid except the last
                    # both of them should be different, except they are cysteines or sulphur containing amino acids
                    self._add_edges(res_one, res_two)
                
                elif res_two.id[1] == 1 and res_one.id[1] != len(self.residues): # first amino acid and any other amino acid except the last
                    self._add_edges(res_one, res_two)
                
                elif res_one.id[1] == len(self.residues) and res_two.id[1] != 1: # last amino acid and any other amino acid except the last
                    # both of them should be different, except they are cysteines or sulphur containing amino acids
                    self._add_edges(res_one, res_two)
                          
                elif res_two.id[1] == len(self.residues) and res_one.id[1] != 1: # last amino acid and any other amino acid except the last
                    # both of them should be different, except they are cysteines or sulphur containing amino acids
                    self._add_edges(res_one, res_two)
                    
                elif sulphur_coords_one and sulphur_coords_two and res_one.get_resname().upper() != "MET" or res_two.get_resname().upper() != "MET": # check all sulfur containing amino acids
                    self._add_edges(res_one, res_two)
                
                else: # neither the first nor last amino acid
                    
                    if res_one.id[1] in (res_two.id[1]-1, res_two.id[1]+1) or res_two.id[1] in (res_one.id[1]-1, res_one.id[1]+1):
                        # sequential amino acids
                        self._add_edges(res_one, res_two)
                    elif res_one.get_resname() != res_two.get_resname():
                        self._add_edges(res_one, res_two)

                    
        except IndexError as e:
            return f"Error message: {e.__repr__()}"
                
    def get_graph(self) -> nx.Graph:
        return self.graph
        
    def _add_edges(self, residue_one, residue_two) -> None:
        
        coordinates_one = [atom for atom in residue_one.get_atoms() if "C" in atom.get_name()]
        coordinates_two = [atom for atom in residue_two.get_atoms() if "C" in atom.get_name()]
        
        sulphur_coords_one: List = [atom for atom in residue_one.get_atoms() if atom.get_name() == 'SG']
        sulphur_coords_two: List = [atom for atom in residue_two.get_atoms() if atom.get_name() == 'SG']
        
        for carbon_one, carbon_two in product(coordinates_one, coordinates_two):
            distance: float = np.linalg.norm(carbon_one - carbon_two)

            if distance >= self.peptide_length.min and distance <= self.peptide_length.max:
                node_one_name, node_two_name = f"{residue_one.get_resname()}_{residue_one.id[1]}", f"{residue_two.get_resname()}_{residue_two.id[1]}"
                self.graph.add_nodes_from([node_one_name, node_two_name])
                self.graph.add_weighted_edges_from([(node_one_name, node_two_name, 3), (node_one_name, node_two_name, 3)])
        
            if sulphur_coords_one and sulphur_coords_two: # Checks the presence of disulfide bonds if SG is present in AA
                # Needs some refining because two methionines cannot form a disulfide bond
                distance: float = round(np.linalg.norm(sulphur_coords_one[0] - sulphur_coords_two[0]), 3)

                if distance >= self.disulfide_length.min and distance <= self.disulfide_length.max:
                    node_one_name, node_two_name = f"{residue_one.get_resname()}_{residue_one.id[1]}", f"{residue_two.get_resname()}_{residue_two.id[1]}"
                    self.graph.add_nodes_from([node_one_name, node_two_name])
                    self.graph.add_weighted_edges_from([(node_one_name, node_two_name, 2), (node_one_name, node_two_name, 2)])
                               
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
            return f"{e.__repr__()}: No Cycle found"

    def get_rsmd(self, reference_structure: np.ndarray, backbone: bool = True) -> float:
        
        try:
            if backbone: # Align atoms based on Ca (backbone) atoms
                ca_coords = np.ndarray([atom.get_coords() for atom in self.residues if atom.get_name() == "CA"])
                
                assert ca_coords.shape == reference_structure.shape
                qcp_superimposer: QCPSuperimposer = QCPSuperimposer()
                qcp_superimposer.set(reference_structure, ca_coords)
                qcp_superimposer.run()
                
                return qcp_superimposer.rms
            else:
                cb_coords = np.ndarray([atom.get_coords() for atom in self.residues if atom.get_name() == "CB"])
                
                assert ca_coords.shape == reference_structure.shape
                qcp_superimposer: QCPSuperimposer = QCPSuperimposer()
                qcp_superimposer.set(reference_structure, cb_coords)
                qcp_superimposer.run()
                
                return qcp_superimposer.rms

        except Exception as e:
            raise e(f"{e.__repr__()}: Coordinates are of different dimensions")
        
    def _get_peptide_descriptor(self):
        return PTDS(sequence=self.get_one_letter_sequence())
    
    def get_mass(self) -> float:
        pep = self._get_peptide_descriptor()
        return pep.molecular_weight()
    
    def isoelectric_point(self, pKscale: str = "EMBOSS") -> float:
        pep = self._get_peptide_descriptor()
        return pep.isoelectric_point(pKscale=pKscale)
    
    def instability_index(self) -> float:
        pep = self._get_peptide_descriptor()
        return pep.instability_index()
    
    def hydrophobicity(self, scale: str = "KyteDoolittle") -> float:
        pep = self._get_peptide_descriptor()
        return pep.hydrophobicity(scale=scale)

    def boman(self) -> float:
        pep = self._get_peptide_descriptor()
        return pep.boman()
    
    def charge(self, pH: float = 7, pKscale: str = "Lehninger") -> float:
        pep = self._get_peptide_descriptor()
        return pep.charge(pH=pH, pKscale=pKscale)
    
    def aliphatic_index(self) -> float:
        pep = self._get_peptide_descriptor()
        return pep.aliphatic_index()

    def cyclization_type(self) -> str:
        pass
    
    def _disulfide_cycle(self) -> bool:
        cycle_members: List = self.get_cycle_members()
        cysteines: List = [member for member in cycle_members if member.get_resname().upper() == "CYS"]
        
        if len(cysteines) < 2:
            return False
        
        for res_one, res_two in combinations(cysteines, 2):
            sulphur_coords_one: List = [atom for atom in res_one.get_atoms() if atom.get_name() == 'S'][0]
            sulphur_coords_two: List = [atom for atom in res_two.get_atoms() if atom.get_name() == 'S'][0]
            
            distance: float = np.linalg.norm(sulphur_coords_one - sulphur_coords_two)
            
            if distance <= self.disulfide_length:
                return True
        
        return False
            
    
    def _n_to_c_cycle(self) -> bool:
        pass