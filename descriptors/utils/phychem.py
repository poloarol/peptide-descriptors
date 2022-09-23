### phychem.py
### Author Paul Arol Wambo

### Provides Physico-Chemical properties of peptide

### TO DO:
### Look into DFT to calculate physico-chemical properties, depending
### how chemically intensive it is.
### This is because non-cannonicals are not modelled here.
### I have no idea, but maybe DFT can be useful here.

from typing import Dict, List

from peptides import Peptide

from Bio.PDB import PPBuilder
from Bio.PDB.Chain import Chain

class PhysicoChemical(object):
    """
    """
    
    def __init__(self, chain: Chain, sequence: str):
        self.chain = chain
        self.peptide = Peptide(sequence=sequence)
    
    def get_mass(self) -> float:
        # Non-canonicals are not modelled accurately, so might need something more robust
        # Actually, I can implement it myself. 
        # 1. Get all atoms in a residue
        # 2. get their monoisotopic 
        # 3. Sum all of them an you get the molecular mass of the protein.
        # if they are non-noncanonicals infer to biopython, else, calculate the mass of 
        # non canonical as stated above
        return self.peptide.molecular_weight()
    
    def isoelectric_point(self, pKscale: str = "EMBOSS") -> float:
        return self.peptide.isoelectric_point(pKscale=pKscale)
    
    def instability_index(self) -> float:
        return self.peptide.instability_index()
    
    def hydrophobicity(self, scale: str = "KyteDoolittle") -> float:
        return self.peptide.hydrophobicity(scale=scale)

    def boman(self) -> float:
        return self.peptide.boman()
    
    def charge(self, pH: float = 7, pKscale: str = "Lehninger") -> float:
        return self.peptide.charge(pH=pH, pKscale=pKscale)
    
    def aliphatic_index(self) -> float:
        return self.peptide.aliphatic_index()
    
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