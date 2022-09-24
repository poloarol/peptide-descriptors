### aas.py
### Author: Paul A. Wambo


from typing import Dict, Set, Tuple

aa_code: Dict[str, str] = {
    "CYS": "C",
    "ASP": "D",
    "SER": "S",
    "GLN": "Q",
    "LYS": "K",
    "ILE": "I",
    "PRO": "P",
    "THR": "T",
    "PHE": "F",
    "ASN": "N",
    "GLY": "G",
    "HIS": "H",
    "LEU": "L",
    "ARG": "R",
    "TRP": "W",
    "ALA": "A",
    "VAL": "V",
    "GLU": "E",
    "TYR": "Y",
    "MET": "M",
}

basic_aa: Tuple[str] = ("HIS", "ARG", "LYS")
acidic_aa: Tuple[str] = ("ASP", "GLU")
all_side_chains: Tuple[str] = list(aa_code.keys())
backbone_atoms: Tuple[str] = (
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


def get_codon_table() -> Dict[str, str]:
    return aa_code


def get_acidic_aas() -> Set[str]:
    return acidic_aa


def get_basic_aas() -> Set[str]:
    return basic_aa


def get_neutral_aas() -> Set[str]:
    neutral_aa: Set[str] = set(all_side_chains) - set(basic_aa + acidic_aa)
    return neutral_aa


def get_backbone_aas() -> Set[str]:
    return backbone_atoms


def get_all_aas() -> Set[str]:
    return set(aa_code.keys())


def get_canonical_mass(amino_acid: str) -> float:
    mass_table_aas: Dict[str, float] = {
        "ALA": 71.03711,
        "CYS": 103.00919,
        "ASP": 115.02694,
        "GLU": 129.04259,
        "PHE": 147.06841,
        "GLY": 57.02146,
        "HIS": 137.05891,
        "ILE": 113.08406,
        "LYS": 128.09496,
        "LEU": 113.08406,
        "MET": 131.04049,
        "GLN": 114.04293,
        "PRO": 97.05276,
        "ASN": 128.05858,
        "ARG": 156.10111,
        "SER": 87.03203,
        "THR": 101.04768,
        "VAL": 99.06841,
        "TRP": 186.07931,
        "TYR": 163.06333,
    }

    return mass_table_aas[amino_acid]
