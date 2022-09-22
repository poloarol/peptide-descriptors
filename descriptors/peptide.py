### peptide.py

"""
Provides methods to obtain peptide descriptors.

    1. cyclic or linear
    2. amino acid sequence
    3. frequency of canonicals vs non-canonicals in the sequence
"""

from audioop import add
from platform import node
from typing import List, Dict, NamedTuple, Tuple
from collections import namedtuple
from itertools import combinations, product

import numpy as np
import networkx as nx

from peptides import Peptide as PTDS
from Bio.PDB import PPBuilder
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.QCPSuperimposer import QCPSuperimposer

aa_code: Dict[str, str] = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
                            'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
                            'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
                            'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

basic_aa: Tuple[str] = ("HIS", "ARG", "LYS")
acidic_aa: Tuple[str] = ("ASP", "GLU")
non_bonding_side_chains: Tuple[str] = ("GLY", "SER", "LEU", "ILE", "MET", "TYR", "ALA", "VAL", "PHE", "TRP", "PRO", "GLN", "ASN")
backbone_atoms: Tuple[str] = ('CA', 'C1', 'HA', 'N', 'HN', 'H', 'C', 'O', 'H1', 'H2', 'H3', 'NH1', 'NH2', 'NH3', 'OXT')

BOND_LENGTHS: NamedTuple = namedtuple('BondLengths', 'min max')

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
        self.disulfide: bool = False
        self.acidic_side_chain: bool = False
        self.basic_side_chain: bool = False
        
        self._has_acidic_side_chain()
        self._has_basic_side_chain()
    
    def get_coordinates(self, element_type = "SG") -> List:
        if element_type == "SG":
            return [np.array(atom.get_coord()) for atom in self.residue if atom.get_name() == element_type]
        elif element_type == "CA":
            return [np.array(atom.get_coord()) for atom in self.residue if atom.get_name() == element_type]
        else: # All carbon types. More useful for non-canonicals or amino acids which can form peptide bonds with their side chain.
            return [np.array(atom.get_coord()) for atom in self.residue if atom.get_name().startswith(element_type)]
        
    def _get_side_chain_atoms(self):
        backbone_atoms: List = ('CA', 'C1', 'HA', 'N', 'HN', 'H', 'C', 'O', 'H1', 'H2', 'H3', 'NH1', 'NH2', 'NH3', 'OXT')
        side_chain_atoms: List[str] = [atom.get_name() for atom in self.residue.get_atoms() if atom.get_name() not in (backbone_atoms)]
        return side_chain_atoms
    
    def _detect_moeity(self, moeity: str):
        side_chain_atoms: List[str] = self._get_side_chain_atoms()
        return True if [element for element in side_chain_atoms if element.startswith(moeity)] else False
        
    def _has_thiol_side_chain(self) -> bool:
        side_chain_atoms: List[str] = self._get_side_chain_atoms()
        return True if [element for element in side_chain_atoms if element == "SG"] else False
    
    def can_form_disulfide(self) -> bool:
        return self._has_thiol_side_chain()
    
    def is_acidic(self) -> bool:
        return self._has_acidic_side_chain()
    
    def is_basic(self) -> bool:
        return self._has_basic_side_chain()
    
    def can_form_disulfide_bond(self) -> bool:
        self.disulfide = self._has_thiol_side_chain()
        return self.disulfide
    
    def can_form_peptide_bond(self) -> bool:
        return True if self.num_peptide_bonds < 2 else False

    def _has_basic_side_chain(self) -> bool:
        if self.residue.get_resname().upper() in basic_aa:
            self.basic_side_chain = True
            return self.basic_side_chain
        elif self.residue.get_resname().upper() in acidic_aa or self.residue.get_resname().upper() in non_bonding_side_chains:
            return False
        else:
            if self._detect_moeity(moeity="N"):
                self.basic_side_chain = True
                return self.basic_side_chain
            return False
    
    def _has_acidic_side_chain(self) -> bool:
        if self.residue.get_resname().upper() in acidic_aa:
            self.acidic_side_chain = True
            return self.acidic_side_chain
        elif self.residue.get_resname().upper() in basic_aa or self.residue.get_resname().upper() in non_bonding_side_chains:
            return False
        else:
            if self._detect_moeity(moeity="O"):
                self.acidic_side_chain = True
            return self.acidic_side_chain
    
    def _increment_num_peptide_bonds(self) -> None:
        self.num_peptide_bonds = self.num_peptide_bonds + 1
    
    def _increment_num_disulfide_bonds(self) -> None:
        self.num_disulfide_bonds = self.num_disulfide_bonds + 1
        
    def form_peptide_bond(self) -> bool:
        return True if self.num_peptide_bonds < 2 else False
    
    def form_disulfide_bond(self) -> bool:
        if self.disulfide:
            return True if self.num_disulfide_bonds < 1 else False
        return False
                
                 
                 
class Peptide(object):
    """
    Peptide (class)
    ---------------
    Provides methods to model peptides and methods 
    to obtain physico-chemical properties.
        
    params:
    chain (Chain): Chain than references peptide sequence
    residues (List): List of residues in the chain
    peptide_length (namedtuple): Min-Max bond length
    disulfide_length (namedtuple): Min-Max bond length
    graph (nx.Graph): Undirected graph
        
    """
    
    def __init__(self, chain: Chain):
        self.chain: Chain = chain
        self.peptide_length: NamedTuple = BOND_LENGTHS(2.899, 4.25) # bond lengths should be ammended
        self.disulfide_length: NamedTuple = BOND_LENGTHS(1.999, 2.05)
        self.ionic_radii_length: NamedTuple = BOND_LENGTHS(0.3, 30)
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
        return {"canonical": self._get_num_canonicals(), "non_canonical": self._get_num_non_canonicals()}
    
    def _get_polypeptide(self) -> PPBuilder:
        return PPBuilder().build_peptides(self.chain)

    def get_dihedral_angles(self) -> List:
        
        dihedrals: List = []
        polypeptide: PPBuilder = self._get_polypeptide()
        
        for _ , peptide in enumerate(polypeptide):
            dihedrals.append(peptide.get_phi_psi_list())
            
        return dihedrals
    
    def get_tau_angles(self) -> List:
        """ List of tau torsions angles for all 4 consecutive Calpha atoms """
        
        tau_angles: List = []
        polypeptide: PPBuilder = self._get_polypeptide()
        
        for _ , peptide in enumerate(polypeptide):
            tau_angles.append(peptide.get_tau_list())
        
        return tau_angles
    
    def get_theta_angles(self) -> List:
        """ List of theta angles for all 3 consecutive Calpha atoms """
        
        theta_angles: List = []
        polypeptide: PPBuilder = self._get_polypeptide()
        
        for _ , peptide in enumerate(polypeptide):
            theta_angles.append(peptide.get_theta_list())
        
        return theta_angles
    
    def _get_side_chain_atoms(self, residue):
        side_chain_atoms: List[str] = [atom.get_name() for atom in residue.get_atoms() if atom.get_name() not in (backbone_atoms)]
        return side_chain_atoms
    
    def _reset_peptide_indices(self) -> None:
        for i, residue in enumerate(self.residues):
            residue.id[1] = i
           
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
                residue_one_node: ResidueNode = None
                residue_two_node: ResidueNode = None
                name_one, name_two = res_one.get_resname().upper(), res_two.get_resname().upper()
                
                # Before adding every node to the graph
                # 1. check if disulfide bond can be created between the two nodes
                # 2. Either create node or not
                
                first_residue_idx, last_residue_idx = self.residues[0].id[1], self.residues[-1].id[1]
                
                res_one_in_graph: bool = self._find_node(attribute="tag", value=f"{name_one}_{res_one.id[1]}")
                res_two_in_graph: bool = self._find_node(attribute="tag", value=f"{name_two}_{res_two.id[1]}")
                                
                if not res_one_in_graph and not res_two_in_graph:
                    residue_one_node = ResidueNode(name=name_one, residue=res_one)
                    residue_two_node = ResidueNode(name=name_two, residue=res_two)
                    
                elif not res_one_in_graph:
                    residue_one_node = self._get_node(attribute="tag", value=f"{name_one}_{res_one.id[1]}")
                    if residue_two_node.num_peptide_bonds < 2: # neutral amino acid
                        residue_one_node = ResidueNode(name=name_one, residue=res_one)
                    elif name_one =="CYS" and name_two == "CYS":
                            residue_one_node = ResidueNode(name=name_one, residue=res_one)
                elif not res_two_in_graph:
                    if residue_one_node.num_peptide_bonds < 2:
                        residue_two_node = ResidueNode(name=name_two, residue=res_two)
                    elif name_one =="CYS" and name_two == "CYS":
                        residue_two_node = ResidueNode(name=name_two, residue=res_two)
                else:
                    residue_one_node = self._get_node(attribute="tag", value=f"{name_one}_{res_one.id[1]}")
                    residue_two_node = self._get_node(attribute="tag", value=f"{name_two}_{res_two.id[1]}")
                    if name_one in acidic_aa and name_two in basic_aa:
                        self._add_ionic_interaction(residue_one_node, residue_two_node)
                    elif name_two in basic_aa and name_one in acidic_aa:
                        self._add_ionic_interaction(residue_one_node, residue_two_node)
                                                                
                # MET does not form disulfide bonds
                if (residue_one_node.can_form_disulfide() and residue_two_node.can_form_disulfide()) and (name_one != "MET" and name_two != "MET"):
                    self._add_edges(residue_one_node, residue_two_node, edge_type="disulfide")
                
                if(res_one.id[1] == first_residue_idx or res_two.id[1] == first_residue_idx) and (res_one.id[1] == last_residue_idx or res_two.id[1] == last_residue_idx): # First and last residues
                    self._add_edges(residue_one_node, residue_two_node, edge_type="peptide")
                
                elif (res_one.id[1] == first_residue_idx) and (res_two.id[1] != last_residue_idx) and (res_one.id[1] not in (res_two.id[1]-1, res_one.id[1]+1)): 
                    self._add_edges(residue_one_node, residue_two_node, edge_type="peptide")
                
                elif (res_two.id[1] == first_residue_idx) and (res_one.id[1] != last_residue_idx) and (res_one.id[1] not in (res_two.id[1]-1, res_one.id[1]+1)):
                    self._add_edges(residue_one_node, residue_two_node, edge_type="peptide")
                
                elif (res_one.id[1] == last_residue_idx) and (res_two.id[1] != first_residue_idx) and (res_one.id[1] not in (res_two.id[1]-1, res_one.id[1]+1)):
                    self._add_edges(residue_one_node, residue_two_node, edge_type="peptide")
                
                elif (res_two.id[1] == last_residue_idx) and (res_one.id[1] != first_residue_idx) and (res_one.id[1] not in (res_two.id[1]-1, res_one.id[1]+1)):
                    self._add_edges(residue_one_node, residue_two_node, edge_type="peptide")  
                    
                else:
                    if res_one.id[1] in (res_two.id[1]-1, res_two.id[1]+1) or res_two.id[1] in (res_one.id[1]-1, res_one.id[1]+1):
                        # sequential amino acids
                        self._add_edges(residue_one_node, residue_two_node, edge_type="peptide")
                    # elif res_one.get_resname() != res_two.get_resname():
                    self._add_edges(residue_one_node, residue_two_node, edge_type="peptide")
                           

        except IndexError as e:
            return f"Error message: {e.__repr__()}"
                
    def get_graph(self) -> nx.Graph:
        return self.graph
        
    def _check_disulfide_bond(self, coordinates_one: np.array, coordinates_two: np.array) -> bool:
        distance: float = np.linalg.norm(coordinates_one - coordinates_two)
        return distance >= self.disulfide_length.min and distance <= self.disulfide_length.max
    
    def _check_peptide_bond(self, coordinates_one: np.array, coordinates_two: np.array) -> bool:
        distance: float = np.linalg.norm(coordinates_one - coordinates_two)
        return distance >= self.peptide_length.min and distance <= self.peptide_length.max
    
    def _theta(self, v: np.array, w: np.array) -> float:
        return np.arccos(v.dot(w)/(np.linalg.norm(v)*np.linalg.norm(w)))
    
    def _add_edge(self, residue_one, residue_two, bond_type: str = "peptide") -> None:
        if bond_type == "peptide":
            self.graph.add_edge(residue_one, residue_two, length=4.0)
        elif bond_type == "disulfide":
            self.graph.add_edge(residue_one, residue_two, length=2.0)
    
    def _find_node(self, attribute: str, value: str) -> bool:
        try:
            return any([node for node in self.graph.nodes(data=True) if node[1][attribute] == value])
        except KeyError as e:
            print(f"Key not found: {e.__repr__()}")
        finally:
            return False
    
    def _get_node(self, attribute: str, value: str) -> ResidueNode:
        try:
            return [node for node in self.graph.nodes(data=True) if node[1][attribute] == value][0]
        except IndexError as e:
            print(f"Key not found: {e.__repr__()}")
        # Find appropriate way to handle the finally clause here
    
    def _add_node(self, node):
        self.graph.add_node(f"{node.name}-{node.residue.id[1]}", name=str(node.residue.id[1]), type='node', tag=f"{node.name}_{node.residue.id[1]}")
        
    def _add_disulfide_edge(self, residue_one: ResidueNode, residue_two: ResidueNode) -> bool:
        if residue_one.form_disulfide_bond() and residue_two.form_disulfide_bond():
            coordinates_one: List = residue_one.get_coordinates()
            coordinates_two: List = residue_two.get_coordinates()
            
            distance: float = np.linalg.norm(coordinates_one, coordinates_two)
            bond_angle: float = self._theta(coordinates_one, coordinates_two)
            
            if (distance >= self.disulfide_length.min and distance <= self.disulfide_length.max) and (179.9 <= bond_angle <= 180.9):
                self._add_node(residue_one)
                self._add_node(residue_two)
                self._add_edge(residue_one, residue_two)
                return True
            return False
        return False
    
    def _add_peptide_edge(self, residue_one: ResidueNode, residue_two: ResidueNode):
        if residue_one.can_form_peptide_bond() and residue_two.can_form_peptide_bond():
            
            coordinates_one: List
            coordinates_two: List
            
            add_node: bool = False
            
            # if either of residue_one or residue_two cannot form peptide bonds with their side chain, 
            # we only get the Ca carbon else, we check all other carbons
            if residue_one.name in non_bonding_side_chains or residue_two.name in non_bonding_side_chains:
                coordinates_one: List = residue_one.get_coordinates(element_type="CA")
                coordinates_two: List = residue_two.get_coordinates(element_type="CA")
            else:
                coordinates_one: List = residue_one.get_coordinates(element_type="C")
                coordinates_two: List = residue_two.get_coordinates(element_type="C")
            
            for carbon_one, carbon_two in product(coordinates_one, coordinates_two):
                distance = np.linalg.norm(carbon_one - carbon_two)
                bond_angle = self._theta(carbon_one, carbon_two)
                                            
                if (distance >= self.peptide_length.min and distance <= self.peptide_length.max) and ((-0.1 <= bond_angle <= 0.1) or (179.9 <= bond_angle <= 180.9)):
                    add_node = True
                    break

            if add_node:
                residue_one_node: ResidueNode = None
                residue_two_node: ResidueNode = None
                
                node_one_present: bool = self._find_node(attribute="tag", value=f"{residue_one.name}_{residue_one.residue.id[1]}")
                node_two_present: bool = self._find_node(attribute="tag", value=f"{residue_two.name}_{residue_two.residue.id[1]}")
                
                if not node_one_present:
                    self._add_node(residue_one)
                
                if not node_two_present:
                    self._add_node(residue_two)
                    
                residue_one_node = self._get_node(attribute="tag", value=f"{residue_one.name}_{residue_one.residue.id[1]}")
                residue_two_node = self._get_node(attribute="tag", value=f"{residue_two.name}_{residue_two.residue.id[1]}")

                if (residue_one.is_acidic() and residue_two.is_acidic()) and residue_one.is_basic() and residue_two.is_basic():
                    pass
                else:  
                    if residue_one_node and residue_two_node:
                        if  residue_one.can_form_peptide_bond() and residue_two.can_form_peptide_bond():
                            if residue_one.form_peptide_bond() and residue_two.form_peptide_bond():
                                self._add_edge(residue_one_node[0], residue_two_node[0])
                       

    def _add_ionic_interaction(self, residue_one, residue_two) -> bool:
        # calculate ionic radius.
        # make sure it falls between range
        pass
                
    def _add_edges(self, residue_one: ResidueNode, residue_two: ResidueNode, edge_type: str) -> None:
        if edge_type == "disulfide":
            self._add_disulfide_edge(residue_one, residue_two)
        elif edge_type == "peptide":
            self._add_peptide_edge(residue_one, residue_two)
        elif edge_type == "ionic":
            self._add_ionic_interaction(residue_one, residue_two)
         
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
                ca_coords = np.ndarray([atom.get_coord() for atom in self.residues if atom.get_name() == "CA"])
                
                assert ca_coords.shape == reference_structure.shape
                qcp_superimposer: QCPSuperimposer = QCPSuperimposer()
                qcp_superimposer.set(reference_structure, ca_coords)
                qcp_superimposer.run()
                
                return qcp_superimposer.rms
            else:
                cb_coords = np.ndarray([atom.get_coord() for atom in self.residues if atom.get_name() == "CB"])
                
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