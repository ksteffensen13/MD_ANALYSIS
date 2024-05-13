from pymol import cmd
from pymol import math
import numpy as np
import math
import sys
import math


# Extract command-line arguments. Here, we can set the start/end distance values for the docking function every time we call this script
if len(sys.argv) != 14:
    print(
        "Usage: python tpm_translations.py <mutant> <Factin_or_Gactin> <nucleotide_state> <cation> <reference_structure> <path_for_reference_structure> <chains_to_delete> <tpm_chains> <path_for_translation_output> <path_for_pymol_axes_script> <path_to_cluster_chain> <length_of_filament_to_build> <reference_starting_actin_chain>")
    sys.exit(1)

mutant = sys.argv[1]
Factin_or_Gactin = sys.argv[2]
nucleotide_state = sys.argv[3]
cation = sys.argv[4]
reference_structure = sys.argv[5]
path_for_reference_structure = sys.argv[6]
chains_to_delete = sys.argv[7]
tpm_chains = sys.argv[8]
path_for_translation_output = sys.argv[9]
path_for_pymol_axes_script = sys.argv[10]
path_to_cluster_chain = sys.argv[11]
length_of_filament_to_build = int(sys.argv[12])
reference_starting_actin_chain = sys.argv[13]

actin_chains = ('A' 'B' 'C' 'D' 'E' 'F' 'G' 'H' 'I' 'J' 'K' 'L' 'M' 'N' 'O' 'P' 'Q' 'R' 'S' 'T' 'U' 'V' 'W' 'X' 'Y' 'Z'
                '1' '2' '3' '4' '5' '6' '7' '8' '9')

reference_starting_actin_chain_index = actin_chains.index(reference_starting_actin_chain)

if nucleotide_state == 'ADP-Pi' or nucleotide_state == 'ADP':
    nucleotide = 'ADP'
else:
    nucleotide = 'ATP'


# create mutant filament from cluster analysis
# want top 5 clusters, so set range 1-6 (not inclusive of 6)
for i in range(1, 6):
    # load and isolate chains
    cmd.load(f'{path_for_reference_structure}{reference_structure}.pdb')
    #remove everything but actin filament
    cmd.select(f'c. {chains_to_delete}+{tpm_chains}')
    cmd.remove('sele')

    translation_output_per_cluster_path = f'{path_for_translation_output}cluster{i}/'
    cmd.load(f'{path_to_cluster_chain}{mutant}_combined_cluster{i}.pdb')

    s = length_of_filament_to_build
    for index, actin_chain in enumerate(
            actin_chains[reference_starting_actin_chain_index:reference_starting_actin_chain_index + s]):
        # pymol and gromacs are kind of dumb together. Pymol won't let me set chains with multiple letters/numbers (ie AA)
        # but GROMACS requires protein/ligands to be in separate chains, which instantly caps the number of chains we can analyze
        # so the only way to get around that cap is to group all similar ligands within the same chain (ie all ADP in one chain)
        # BUT gromacs will then combine all of the ligands within that into just 1 ligand unless we give them unique residue numbers
        #so we have to create new chains, give nucleotide/cation unique numbering, then group all nucleotides into a chain
        nucleotide_numbering = 400 + index
        cation_numbering = 500 + index
        PO4_numbering = 600 + index

        # create copies of mutant chain
        cmd.create(f"{mutant}_combined_cluster{i}_chain{actin_chain}", f"{mutant}_combined_cluster{i}")
        # structural superposition
        cmd.super(f'{mutant}_combined_cluster{i}_chain{actin_chain}', f'{reference_structure} & chain {actin_chain}')
        # change the chain to match the reference structure chain
        cmd.alter(selection=f'{mutant}_combined_cluster{i}_chain{actin_chain}', expression=f"chain='{actin_chain}'")
        #give nucleotide and cation within the current chain a unique number, which increases each iteration
        cmd.alter(selection=f'chain {actin_chain} and resn {nucleotide}', expression=f"resi='{nucleotide_numbering}'")
        cmd.alter(selection=f'chain {actin_chain} and resn {cation}', expression=f"resi='{cation_numbering}'")
        if nucleotide_state == 'ADP-Pi':
            cmd.alter(selection=f'chain {actin_chain} and resn PO4', expression=f"resi='{PO4_numbering}'")

    # change the chain for ALL nucletoides and cations, so that all nucleotide are in same chain (with different residue numbers from above)
    #also change chain for all cations so that they're in the same chain
    #First, calculate what chain the nucleotide will be in (number of actin chains + 1)
    nucleotide_index = s + 1
    #determine what letter corresponds to that
    nucleotide_chain = actin_chains[nucleotide_index]
    #repeat for cation
    cation_index = s + 2
    cation_chain = actin_chains[cation_index]
    #actually alter the chains
    cmd.alter(selection=f"resn {nucleotide}", expression=f"chain='{nucleotide_chain}'")
    cmd.alter(selection=f"resn {cation}", expression=f"chain='{cation_chain}'")
    if nucleotide_state == 'ADP-Pi':
        PO4_index = s + 3
        PO4_chain = actin_chains[PO4_index]
        cmd.alter(selection=f"resn PO4", expression=f"chain='{PO4_chain}'")
    #remove the initial actin chain that was copied
    cmd.delete(f"{mutant}_combined_cluster{i}")
    cmd.delete(f'{reference_structure}')
    cmd.save(f'{translation_output_per_cluster_path}{mutant}_base_filament_cluster{i}.pdb')
    cmd.delete('all')

