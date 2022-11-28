# peptide-descriptors
Provides methods to handle PDB files, to obtain both molecular and strcutural descriptors for peptides.

## Current Status

1. Determine is a peptide is cyclic or linear
2. Extract side chain angles: phi, psi and tau, theta
3. Count number of non-cannonicals vs cannonicals
4. For peptides containing only cannonical amino acids, it calculates:
    - Isoelectric Point at a given pH and pK scale
    - Boman Index at a given scale
    - Aliphatic Index at a given scale
    - Charge at a given pH and pK scale
    - Hydrophobicity at a given scale

## Currently working on

1. Calculate the coulombic charge of peptide
3. Implementing DFT to calculate properties of non-cannonical amino acids
