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
                peptide: Peptide = Peptide(chain=chain)
                peptide.generate_peptide_graph()
                graph = peptide.get_graph()
                
                physico_chemical = PhysicoChemical(chain=chain, sequence=peptide.get_one_letter_sequence())
                
                print(pdb, graph, peptide.get_sequence(), physico_chemical.get_mass())
                # print(pdb, peptide.is_cyclic(), peptide.get_sequence())
                # graph = peptide.get_graph()
                
                # if peptide.is_cyclic():
                #     nx.draw(graph, with_labels=True)
                #     plt.savefig(os.path.join(os.getcwd(), f"graphs/cycle/{pdb}_{i}.png"))
                #     plt.clf()
                # else:
                #     nx.draw(graph, with_labels=True)
                #     plt.savefig(os.path.join(os.getcwd(), f"graphs/no_cycle/{pdb}_{i}.png"))
                #     plt.clf()
            os.remove(file_name)
