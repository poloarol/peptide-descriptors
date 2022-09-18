# app.py

import os
import tempfile
import urllib
from typing import List

import networkx as nx
import matplotlib.pyplot as plt

from Bio.PDB import PDBParser
from descriptors.peptide import Peptide

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
                
            pep_chains = load_structure(file_name)
            for i, chain in enumerate(pep_chains):
                if pdb == "4MQ9":
                    print(chain)
                if len(list(chain.get_residues())) > 3:
                    peptide = Peptide(chain=chain)
                    graph = peptide.get_graph()
                    
                    print(peptide.get_sequence())
                        # print(peptide.get_one_letter_sequence())
                        # print(peptide.get_num_canonical_non_canonical())
                        # print(peptide.get_dihedral_angles())
                        # print(peptide.get_tau_angles())
                        # print(peptide.get_theta_angles())
                        
                    if peptide.is_cyclic():
                        print(peptide.get_cycle_members())
                            
                            # print(peptide.get_mass())
                            # print(peptide.isoelectric_point())
                            # print(peptide.hydrophobicity())
                            # print(peptide.boman())
                            # print(peptide.charge())
                            # print(peptide.aliphatic_index())
                            
                        nx.draw(graph, with_labels=True)
                        plt.savefig(os.path.join(os.getcwd(), f"graphs/cycle/{pdb}_{i}.png"))
                        plt.clf()
                    else:
                        nx.draw(graph, with_labels=True)
                        plt.savefig(os.path.join(os.getcwd(), f"graphs/no_cycle/{pdb}_{i}.png"))
                        plt.clf()
            os.remove(file_name)