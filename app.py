# app.py

import os
import tempfile
import urllib
from typing import List

import networkx as nx
import matplotlib.pyplot as plt

from Bio.PDB import PDBParser
from descriptors.peptide import Peptide
from descriptors.utils.phychem import PhysicoChemical

def load_structure(filename: str):
    try:
        complex = parser.get_structure("XXX", filename)
        peptide_chain: List = [chain for chain in complex.get_chains() if len(list(chain.get_residues())) < 51]
        
        return peptide_chain
    except Exception as e:
        pass



parser: PDBParser = PDBParser(PERMISSIVE=1)

if __name__ == '__main__':
    with open(os.path.join(os.getcwd(), "data/pdbs.txt"), 'r') as lines:
        for line in lines:
            pdb = line.strip("\n").strip()
            file_name: str
            
            with tempfile.NamedTemporaryFile(delete=False) as tmp:
                url: str = f"https://files.rcsb.org/download/{pdb}.pdb"
                urllib.request.urlretrieve(url, tmp.name)
                file_name = tmp.name
                
            chains: List = load_structure(filename=file_name)
            for i, chain in enumerate(chains):
                cycle_found: bool

                peptide: Peptide = Peptide(chain=chain)
                peptide.generate_peptide_graph()
                graph = peptide.get_graph()

                physico_chemical = PhysicoChemical(chain=chain, residues=chain, sequence=peptide.get_one_letter_sequence())

                try:
                    cycle_found = any(nx.find_cycle(graph, orientation="ignore"))
                except:
                    cycle_found = False
                
                print(pdb, cycle_found, peptide.get_sequence(), physico_chemical.get_mass())
            os.remove(file_name)
