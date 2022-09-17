# app.py

import os
from typing import List

from Bio.PDB import PDBParser

def load_structure(filename: str):
    complex = parser.get_structure("XXX", filename)
    peptide_chain: List = [chain for chain in complex.get_chains() if len(list(chain.get_residues())) < 51][0]
    
    return peptide_chain



parser: PDBParser = PDBParser(PERMISSIVE=1)

if __name__ == '__main__':
    path: str = os.path.join(os.getcwd(), "data/2c9t.pdb")
    pep_chain = load_structure(path)
    
    print(pep_chain, type(pep_chain))