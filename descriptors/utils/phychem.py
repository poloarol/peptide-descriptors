### phychem.py
### Author Paul Arol Wambo

### Provides Physico-Chemical properties of peptide

### TO DO:
### Look into DFT to calculate physico-chemical properties, depending
### how chemically intensive it is.
### This is because non-cannonicals are not modelled here.
### I have no idea, but maybe DFT can be useful here.

from typing import Dict, List

from Bio.PDB import PPBuilder
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from molmass import Formula
from peptides import Peptide

from .aas import get_codon_table


class PhysicoChemical(object):
    """ """

    def __init__(self, chain: Chain, residues: List[Residue], sequence: str):
        self.chain = chain
        self.residues = residues
        self.peptide = Peptide(sequence=sequence)
        self.codon_table: Dict[str, str] = get_codon_table()

    def _get_all_atoms(self, compounds: List) -> Dict[str, int]:
        """
        Given a list of atoms, as descibed within the BioPython,
        counts the number of atoms present within the molecule (C, H, N, O, S),
        and returns their count as a dictionary.

        param
        -----
        compound: List of atoms

        return
        ------
        abundance of each atom in the compound

        """
        elmts: Dict[str, int] = {"C": 0, "H": 0, "O": 0, "S": 0, "N": 0}

        for compound in compounds:
            for ind_atom in compound.get_name():
                if ind_atom.isalpha():
                    if ind_atom == "H":
                        elmts["H"] += 1
                    elif ind_atom == "C":
                        elmts["C"] += 1
                    elif ind_atom == "O":
                        elmts["O"] += 1
                    elif ind_atom == "S":
                        elmts["S"] += 1
                    elif ind_atom == "N":
                        elmts["N"] += 1

        return elmts

    def get_mass(self) -> float:
        """
        Calculates the mass of the amino acid sequence, by extracting all atoms
        and multiplying their count by their monoisotopic mass

        return
        -----
        mass (float)
        """
        mass: float = 0.0
        for residue in self.residues:
            elements: Dict[str, int] = self._get_all_atoms(residue.get_atoms())
            tmp: str = "".join(
                [f"{x}{y}" for x, y in zip(elements.keys(), elements.values())]
            )
            formulae: str = "".join([i for i in tmp if not i.isdigit()])

            mass = mass + Formula(formulae).isotope.mass

        return mass

    def isoelectric_point(self, pKscale: str = "EMBOSS") -> float:
        """Provides the peptide isoelectric point on a specific scale"""
        return self.peptide.isoelectric_point(pKscale=pKscale)

    def instability_index(self) -> float:
        """Provides the peptides instability index"""
        return self.peptide.instability_index()

    def hydrophobicity(self, scale: str = "KyteDoolittle") -> float:
        """Provides the peptide hydrophobicity on a specific scale"""
        return self.peptide.hydrophobicity(scale=scale)

    def boman(self) -> float:
        return self.peptide.boman()

    def charge(self, pH: float = 7, pKscale: str = "Lehninger") -> float:
        """Provides the peptide charge on a specific scale"""
        return self.peptide.charge(pH=pH, pKscale=pKscale)

    def aliphatic_index(self) -> float:
        """Provides the peptide's aliphatic index"""
        return self.peptide.aliphatic_index()

    def _get_num_canonicals(self) -> int:
        """determines the number of canonical amino acids"""
        sequence: str = self.get_one_letter_sequence()
        num_non_canonicals: int = self._get_num_non_canonicals()
        return len(sequence) - num_non_canonicals

    def _get_num_non_canonicals(self) -> int:
        """determines the number of non-canonical amino acids"""
        sequence: str = self.get_one_letter_sequence()
        return sequence.count("X")

    def get_num_canonical_non_canonical(self) -> Dict[str, int]:
        """Provides a count of canonical vs non-canonical amino acids"""
        return {
            "canonical": self._get_num_canonicals(),
            "non_canonical": self._get_num_non_canonicals(),
        }

    def _get_polypeptide(self) -> PPBuilder:
        """
        Builds a polypeptide chain to be ready for BioPython analysis.
        """
        return PPBuilder().build_peptides(self.chain)

    def get_dihedral_angles(self) -> List:
        """
        Generates a list of dihedral angles

        return:
        List
        """

        dihedrals: List = []
        polypeptide: PPBuilder = self._get_polypeptide()

        for _, peptide in enumerate(polypeptide):
            dihedrals.append(peptide.get_phi_psi_list())

        return dihedrals

    def get_tau_angles(self) -> List:
        """
        List of tau torsions angles for all 4 consecutive Calpha atoms
        """

        tau_angles: List = []
        polypeptide: PPBuilder = self._get_polypeptide()

        for _, peptide in enumerate(polypeptide):
            tau_angles.append(peptide.get_tau_list())

        return tau_angles

    def get_theta_angles(self) -> List:
        """
        List of theta angles for all 3 consecutive Calpha atoms
        """

        theta_angles: List = []
        polypeptide: PPBuilder = self._get_polypeptide()

        for _, peptide in enumerate(polypeptide):
            theta_angles.append(peptide.get_theta_list())

        return theta_angles
    
    def calculate_coulombic_charge(self) -> float:
        pass
