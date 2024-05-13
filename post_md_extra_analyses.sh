#extra distance analyses after looking at structure/key areas

#!/bin/bash
#SBATCH --time=0-12:00:00           # time limit (D-HH:MM:SS)
#SBATCH --account=def-jfdawson
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=0
#SBATCH --output=R312_extra_analyses.out
#SBATCH --job-name=R312_extra_analyses
#SBATCH --mail-user=ksteffen@uoguelph.ca
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=ALL

module load StdEnv/2020 gcc/9.3.0 cuda/11.0 openmpi/4.0.3 gromacs/2021.2
export GMXLIB=/home/ksteffen/gromacs-2021.5/share/top

###################################ENTER MANUALLY##################################################################################

mutants=("WT" "R312C" "R312H")
Factin_or_Gactin="Factin"
directory_path="/home/ksteffen/scratch/MUTANT_FILAMENTS"
#path to files from network analysis tutorial that are used to generate our network
network_tutorial_path="/home/ksteffen/NETWORK_ANALYSIS/network-tutorial/tutorial_files/common"
#path to networkview plugin scripts for network analysis
path_to_networkview="/home/ksteffen/vmd-1.9.4a55/plugins/noarch/tcl/networkview1.41/"
# usually 375, might change due to deletions like F90delta or AL
number_of_residues=375
# chain you want to analyze
number_of_chains=5
target_chain='C'
# what nucleotide is in structure? ATP, ADP-Pi, or ADP
nucleotide_state='ADP-Pi'
# bound cation: MG or CAL
cation='MG'
# simulation length in nanoseconds
simulation_length=100
# check index files from the equilibration step to get the number of indices
# THIS DOES NOT INCLUDE THE INDICES YOU MADE LIKE "PROTEIN_MG_ADP", JUST THE DEFAULT ONES
#to figure this out, take your STARTING PDB, and use gmx make_ndx -f STARTING_PDB.pdb -o TEST.ndx and look at the last index number, then cancel the command
number_of_indices=19
#set time you want to consider the simulation at equilibrium (usually 70% of simulation time)
eq_start_time=70



#BELOW HERE IS AUTOMATED AND ONLY REQUIRES CHANGES IF YOURE USING A VERY DIFFERENT SYSTEM


##########################################FOR CALCULATIONS##############################################################

chains=('A' 'B' 'C' 'D' 'E' 'F' 'G' 'H' 'I' 'J' 'K')

# Find the index of the target_chain in the array
# this looks through the array 'chains', at every index of the array, and if that entry is the same as 'target chain',
# then it sets that index number to be the variable 'target_chain_number' which can be called later on
#note: 'A' in chains is index number 0 (numbering starts at 0 not 1), so to get target chain number for calculations have to add 1
for i in "${!chains[@]}"; do
    if [[ "${chains[i]}" == "$target_chain" ]]; then
        target_chain_number=$((i + 1))
        target_chain_index_number="$i"
        break
    fi
done

#based on the target chain, calculate the letters and numbers of the neighbouring chains which can be called for analyses later
if [[ "$target_chain_number" -gt 2 ]]; then
  interacting_chain_number_minus2=$((target_chain_index_number- 2))
  interacting_chain_minus2="${chains[$interacting_chain_number_minus2]}"
  interacting_chain_number_minus1=$((target_chain_index_number - 1))
  interacting_chain_minus1="${chains[$interacting_chain_number_minus1]}"
  interacting_chain_number_plus1=$((target_chain_index_number + 1))
  interacting_chain_plus1="${chains[$interacting_chain_number_plus1]}"
  interacting_chain_number_plus2=$((target_chain_index_number + 2))
  interacting_chain_plus2="${chains[$interacting_chain_number_plus2]}"
elif [[ "$target_chain_number" -eq 2 ]]; then
  interacting_chain_number_minus1=1
  interacting_chain_minus1='A'
  interacting_chain_number_plus1=3
  interacting_chain_plus1='C'
  interacting_chain_number_plus2=4
  interacting_chain_plus2='D'
elif [[ "$target_chain_number" -eq 1 ]]; then
  interacting_chain_number_plus1=2
  interacting_chain_plus1='B'
  interacting_chain_number_plus2=3
  interacting_chain_plus2='C'
fi


# to automatically generate index files for analyzing specific chains, we need to calculate the index numbers
# This depends on how you set up the starting structure.
# usually in Pymol (depending on the starting PDB), the order of chains goes Protein, ADP, MG, PO4
# In that case, protein indexes will appear first, then nucleotide (ADP), then MG, then (if ADP-Pi nucleotide state) PO4 last
# based on that order, the below calculates index numbers first of protein, THEN ADP assuming it comes right after protein
# THEN MG assuming it comes right after ADP, then (if it's in the structure) PO4 last
# You will need to change the order of these calculations if your sequence in your structure is different

target_chain_start=$((target_chain_number * number_of_residues - 374))
target_chain_end=$((target_chain_number * number_of_residues))
target_chain_residues="$target_chain_start-$target_chain_end"

# the nucleotide residue number in the index file comes after all of the protein residues
# So, we calculate the number of the last residue, then add a number corresponding to the target chain
# for example, if we want the ADP of chain C (chain #3), in a 5 protomer filament, then ADP index=(5*375) + 3 = 1878
# this is because all 5 protein chains come first, then each ADP for each chain
nucleotide_index_number=$((number_of_chains * number_of_residues + target_chain_number))
# similar process for cation. Comes after all ADP, so calculate the number of protein residues, then add the total number of ADP, then add the target chain number
cation_index_number=$((number_of_chains * number_of_residues + number_of_chains + target_chain_number))
# if ADP-Pi, the PO4 usually comes last
if [[ "$nucleotide_state" == "ADP-Pi" ]]; then
  po4_index_number=$((number_of_chains * number_of_residues + 2 * number_of_chains + target_chain_number))
  nucleotide='ADP'
elif [[ "$nucleotide_state" == "ADP" ]]; then
  nucleotide='ADP'
else
  nucleotide='ATP'
fi

start_index=$((number_of_indices + 1))

# calculate subdomain residue numbers for the desired chain
subdomain1_residue_start1=$((1 + (number_of_residues) * (target_chain_number - 1)))
subdomain1_residue_end1=$((32 + (number_of_residues) * (target_chain_number - 1)))
subdomain1_residue_start2=$((70 + (number_of_residues) * (target_chain_number - 1)))
subdomain1_residue_end2=$((137 + (number_of_residues) * (target_chain_number - 1)))
subdomain1_residue_start3=$((337 + (number_of_residues) * (target_chain_number - 1)))
subdomain1_residue_end3=$(((number_of_residues) + (number_of_residues) * (target_chain_number - 1)))

subdomain2_residue_start=$((33 + (number_of_residues) * (target_chain_number - 1)))
subdomain2_residue_end=$((69 + (number_of_residues) * (target_chain_number - 1)))

subdomain3_residue_start1=$((138 + (number_of_residues) * (target_chain_number - 1)))
subdomain3_residue_end1=$((178 + (number_of_residues) * (target_chain_number - 1)))
subdomain3_residue_start2=$((274 + (number_of_residues) * (target_chain_number - 1)))
subdomain3_residue_end2=$((336 + (number_of_residues) * (target_chain_number - 1)))

subdomain4_residue_start=$((179 + (number_of_residues) * (target_chain_number - 1)))
subdomain4_residue_end=$((273 + (number_of_residues) * (target_chain_number - 1)))

# generate text strings for subdomains that can be called for making index files
target_chain_subdomain1_residues1="$subdomain1_residue_start1-$subdomain1_residue_end1"
target_chain_subdomain1_residues2="$subdomain1_residue_start2-$subdomain1_residue_end2"
target_chain_subdomain1_residues3="$subdomain1_residue_start3-$subdomain1_residue_end3"

target_chain_subdomain2_residues="$subdomain2_residue_start-$subdomain2_residue_end"

target_chain_subdomain3_residues1="$subdomain3_residue_start1-$subdomain3_residue_end1"
target_chain_subdomain3_residues2="$subdomain3_residue_start2-$subdomain3_residue_end2"

target_chain_subdomain4_residues="$subdomain4_residue_start-$subdomain4_residue_end"

##########################################ACTUAL SCRIPT################################################################

for mutant in "${mutants[@]}"; do
  for replicate in Replicate1 Replicate2 Replicate3; do
    cd "$directory_path"/"$mutant"/"$replicate"/CONTINUATION

    if [ "$replicate" == "Replicate1" ]; then
        rep="rep1"
    elif [ "$replicate" == "Replicate2" ]; then
        rep="rep2"
    elif [ "$replicate" == "Replicate3" ]; then
        rep="rep3"
    fi

    mkdir DISTANCES2/
    mkdir DISTANCES2/INTERPROTOMER2/

    ChainC_D244=$((244 + (number_of_residues * (target_chain_number - 1))))
    ChainA_M325=$((325 + (number_of_residues * (target_chain_number - 3))))
    ChainA_R290=$((290 + (number_of_residues * (target_chain_number - 3))))
    ChainA_I287=$((287 + (number_of_residues * (target_chain_number - 3))))

    ChainC_E205=$((205 + (number_of_residues * (target_chain_number - 1))))
    #ChainA_I287=$((287 + (number_of_residues * (target_chain_number - 3))))

    ChainC_R62=$((62 + (number_of_residues * (target_chain_number - 1))))
    ChainA_D288=$((288 + (number_of_residues * (target_chain_number - 3))))

    ChainC_G63=$((63 + (number_of_residues * (target_chain_number - 1))))
    ChainA_D286=$((286 + (number_of_residues * (target_chain_number - 3))))
    #ChainA_D288=$((288 + (number_of_residues * (target_chain_number - 3))))

    ChainC_I289=$((289 + (number_of_residues * (target_chain_number - 1))))
    ChainC_D292=$((292 + (number_of_residues * (target_chain_number - 1))))
    ChainC_Y166=$((166 + (number_of_residues * (target_chain_number - 1))))


    ChainE_D244=$((244 + (number_of_residues * (target_chain_number + 1))))
    ChainC_M325=$((325 + (number_of_residues * (target_chain_number - 1))))
    ChainC_R290=$((290 + (number_of_residues * (target_chain_number - 1))))
    ChainC_I287=$((287 + (number_of_residues * (target_chain_number - 1))))

    ChainE_E205=$((205 + (number_of_residues * (target_chain_number + 1))))
    #ChainC_I287=$((287 + (number_of_residues * (target_chain_number - 3))))

    ChainE_R62=$((62 + (number_of_residues * (target_chain_number + 1))))
    ChainC_D288=$((288 + (number_of_residues * (target_chain_number - 1))))

    ChainE_G63=$((63 + (number_of_residues * (target_chain_number + 1))))
    ChainC_D286=$((286 + (number_of_residues * (target_chain_number - 1))))
    #ChainC_D288=$((288 + (number_of_residues * (target_chain_number - 3))))


    echo -e 'ri '$ChainC_D244'\n name '$start_index 'ChainC_D244\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX2.ndx -nobackup
    echo -e 'ri '$ChainA_M325'\n name '$((start_index + 1)) 'ChainA_M325\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX2.ndx -n "$mutant"_"$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX2.ndx -nobackup
    echo -e 'ri '$ChainA_R290'\n name '$((start_index + 2)) 'ChainA_R290\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX2.ndx -n "$mutant"_"$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX2.ndx -nobackup
    echo -e 'ri '$ChainA_I287'\n name '$((start_index + 3)) 'ChainA_I287\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX2.ndx -n "$mutant"_"$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX2.ndx -nobackup
    echo -e 'ri '$ChainC_E205'\n name '$((start_index + 4)) 'ChainC_E205\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX2.ndx -n "$mutant"_"$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX2.ndx -nobackup
    echo -e 'ri '$ChainC_R62'\n name '$((start_index + 5)) 'ChainC_R62\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX2.ndx -n "$mutant"_"$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX2.ndx -nobackup
    echo -e 'ri '$ChainA_D288'\n name '$((start_index + 6)) 'ChainA_D288\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX2.ndx -n "$mutant"_"$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX2.ndx -nobackup
    echo -e 'ri '$ChainC_G63'\n name '$((start_index + 7)) 'ChainC_G63\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX2.ndx -n "$mutant"_"$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX2.ndx -nobackup
    echo -e 'ri '$ChainA_D286'\n name '$((start_index + 8)) 'ChainA_D286\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX2.ndx -n "$mutant"_"$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX2.ndx -nobackup
    echo -e 'ri '$ChainC_I289'\n name '$((start_index + 9)) 'ChainC_I289\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX2.ndx -n "$mutant"_"$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX2.ndx -nobackup
    echo -e 'ri '$ChainC_D292'\n name '$((start_index + 10)) 'ChainC_D292\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX2.ndx -n "$mutant"_"$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX2.ndx -nobackup
    echo -e 'ri '$ChainC_Y166'\n name '$((start_index + 11)) 'ChainC_Y166\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX2.ndx -n "$mutant"_"$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX2.ndx -nobackup

    echo -e 'ri '$ChainE_D244'\n name '$((start_index + 12)) 'ChainE_D244\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX2.ndx -n "$mutant"_"$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX2.ndx -nobackup
    echo -e 'ri '$ChainC_M325'\n name '$((start_index + 13)) 'ChainC_M325\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX2.ndx -n "$mutant"_"$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX2.ndx -nobackup
    echo -e 'ri '$ChainC_R290'\n name '$((start_index + 14)) 'ChainC_R290\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX2.ndx -n "$mutant"_"$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX2.ndx -nobackup
    echo -e 'ri '$ChainC_I287'\n name '$((start_index + 15)) 'ChainC_I287\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX2.ndx -n "$mutant"_"$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX2.ndx -nobackup
    echo -e 'ri '$ChainE_E205'\n name '$((start_index + 16)) 'ChainE_E205\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX2.ndx -n "$mutant"_"$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX2.ndx -nobackup
    echo -e 'ri '$ChainE_R62'\n name '$((start_index + 17)) 'ChainE_R62\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX2.ndx -n "$mutant"_"$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX2.ndx -nobackup
    echo -e 'ri '$ChainC_D288'\n name '$((start_index + 18)) 'ChainC_D288\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX2.ndx -n "$mutant"_"$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX2.ndx -nobackup
    echo -e 'ri '$ChainE_G63'\n name '$((start_index + 19)) 'ChainE_G63\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX2.ndx -n "$mutant"_"$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX2.ndx -nobackup
    echo -e 'ri '$ChainC_D286'\n name '$((start_index + 20)) 'ChainC_D286\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX2.ndx -n "$mutant"_"$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX2.ndx -nobackup




    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX2.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/INTERPROTOMER2/"$mutant"_"$Factin_or_Gactin"_Chain_C_D244_Chain_A_M325_"$rep".xvg -ref ChainC_D244 -sel ChainA_M325
    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX2.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/INTERPROTOMER2/"$mutant"_"$Factin_or_Gactin"_Chain_C_D244_Chain_A_R290_"$rep".xvg -ref ChainC_D244 -sel ChainA_R290
    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX2.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/INTERPROTOMER2/"$mutant"_"$Factin_or_Gactin"_Chain_C_D244_Chain_A_I287_"$rep".xvg -ref ChainC_D244 -sel ChainA_I287
    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX2.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/INTERPROTOMER2/"$mutant"_"$Factin_or_Gactin"_Chain_C_E205_Chain_A_I287_"$rep".xvg -ref ChainC_E205 -sel ChainA_I287
    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX2.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/INTERPROTOMER2/"$mutant"_"$Factin_or_Gactin"_Chain_C_R62_Chain_A_D288_"$rep".xvg -ref ChainC_R62 -sel ChainA_D288
    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX2.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/INTERPROTOMER2/"$mutant"_"$Factin_or_Gactin"_Chain_C_G63_Chain_A_D286_"$rep".xvg -ref ChainC_G63 -sel ChainA_D286
    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX2.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/INTERPROTOMER2/"$mutant"_"$Factin_or_Gactin"_Chain_C_G63_Chain_A_D288_"$rep".xvg -ref ChainC_G63 -sel ChainA_D288

    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX2.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/INTERPROTOMER2/"$mutant"_"$Factin_or_Gactin"_Chain_C_D292_Chain_C_Y166_"$rep".xvg -ref ChainC_D292 -sel ChainC_Y166
    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX2.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/INTERPROTOMER2/"$mutant"_"$Factin_or_Gactin"_Chain_C_I289_Chain_C_Y166_"$rep".xvg -ref ChainC_I289 -sel ChainC_Y166

    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX2.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/INTERPROTOMER2/"$mutant"_"$Factin_or_Gactin"_Chain_E_D244_Chain_C_M325_"$rep".xvg -ref ChainE_D244 -sel ChainC_M325
    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX2.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/INTERPROTOMER2/"$mutant"_"$Factin_or_Gactin"_Chain_E_D244_Chain_C_R290_"$rep".xvg -ref ChainE_D244 -sel ChainC_R290
    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX2.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/INTERPROTOMER2/"$mutant"_"$Factin_or_Gactin"_Chain_E_D244_Chain_C_I287_"$rep".xvg -ref ChainE_D244 -sel ChainC_I287
    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX2.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/INTERPROTOMER2/"$mutant"_"$Factin_or_Gactin"_Chain_E_E205_Chain_C_I287_"$rep".xvg -ref ChainE_E205 -sel ChainC_I287
    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX2.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/INTERPROTOMER2/"$mutant"_"$Factin_or_Gactin"_Chain_E_R62_Chain_C_D288_"$rep".xvg -ref ChainE_R62 -sel ChainC_D288
    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX2.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/INTERPROTOMER2/"$mutant"_"$Factin_or_Gactin"_Chain_E_G63_Chain_C_D286_"$rep".xvg -ref ChainE_G63 -sel ChainC_D286
    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX2.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/INTERPROTOMER2/"$mutant"_"$Factin_or_Gactin"_Chain_E_G63_Chain_C_D288_"$rep".xvg -ref ChainE_G63 -sel ChainC_D288


    mkdir DISTANCES2/R312_SD4_TPM_BINDING/


    ChainC_T318=$((318 + (number_of_residues * (target_chain_number - 1))))
    ChainC_I327=$((327 + (number_of_residues * (target_chain_number - 1))))
    ChainC_I317=$((317 + (number_of_residues * (target_chain_number - 1))))
    ChainC_Q314=$((314 + (number_of_residues * (target_chain_number - 1))))
    ChainC_I329=$((329 + (number_of_residues * (target_chain_number - 1))))
    ChainC_M313=$((313 + (number_of_residues * (target_chain_number - 1))))
    ChainC_A310=$((310 + (number_of_residues * (target_chain_number - 1))))
    ChainC_A331=$((331 + (number_of_residues * (target_chain_number - 1))))
    ChainC_M305=$((305 + (number_of_residues * (target_chain_number - 1))))
    ChainC_R335=$((335 + (number_of_residues * (target_chain_number - 1))))
    ChainC_K336=$((336 + (number_of_residues * (target_chain_number - 1))))


    echo -e 'ri '$ChainC_T318'\n name '$start_index 'ChainC_T318\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_R312_SD4_TPM_BINDING_INDEX.ndx -nobackup
    echo -e 'ri '$ChainC_I327'\n name '$((start_index + 1)) '$ChainC_I327\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_R312_SD4_TPM_BINDING_INDEX.ndx -n "$mutant"_"$Factin_or_Gactin"_R312_SD4_TPM_BINDING_INDEX.ndx -nobackup
    echo -e 'ri '$ChainC_I317'\n name '$((start_index + 2)) '$ChainC_I317\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_R312_SD4_TPM_BINDING_INDEX.ndx -n "$mutant"_"$Factin_or_Gactin"_R312_SD4_TPM_BINDING_INDEX.ndx -nobackup
    echo -e 'ri '$ChainC_Q314'\n name '$((start_index + 3)) '$ChainC_Q314\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_R312_SD4_TPM_BINDING_INDEX.ndx -n "$mutant"_"$Factin_or_Gactin"_R312_SD4_TPM_BINDING_INDEX.ndx -nobackup
    echo -e 'ri '$ChainC_I329'\n name '$((start_index + 4)) '$ChainC_I329\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_R312_SD4_TPM_BINDING_INDEX.ndx -n "$mutant"_"$Factin_or_Gactin"_R312_SD4_TPM_BINDING_INDEX.ndx -nobackup
    echo -e 'ri '$ChainC_M313'\n name '$((start_index + 5)) '$ChainC_M313\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_R312_SD4_TPM_BINDING_INDEX.ndx -n "$mutant"_"$Factin_or_Gactin"_R312_SD4_TPM_BINDING_INDEX.ndx -nobackup
    echo -e 'ri '$ChainC_A310'\n name '$((start_index + 6)) '$ChainC_A310\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_R312_SD4_TPM_BINDING_INDEX.ndx -n "$mutant"_"$Factin_or_Gactin"_R312_SD4_TPM_BINDING_INDEX.ndx -nobackup
    echo -e 'ri '$ChainC_A331'\n name '$((start_index + 7)) '$ChainC_A331\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_R312_SD4_TPM_BINDING_INDEX.ndx -n "$mutant"_"$Factin_or_Gactin"_R312_SD4_TPM_BINDING_INDEX.ndx -nobackup
    echo -e 'ri '$ChainC_M305'\n name '$((start_index + 8)) '$ChainC_M305\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_R312_SD4_TPM_BINDING_INDEX.ndx -n "$mutant"_"$Factin_or_Gactin"_R312_SD4_TPM_BINDING_INDEX.ndx -nobackup
    echo -e 'ri '$ChainC_R335'\n name '$((start_index + 9)) '$ChainC_R335\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_R312_SD4_TPM_BINDING_INDEX.ndx -n "$mutant"_"$Factin_or_Gactin"_R312_SD4_TPM_BINDING_INDEX.ndx -nobackup
    echo -e 'ri '$ChainC_K336'\n name '$((start_index + 10)) '$ChainC_K336\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_R312_SD4_TPM_BINDING_INDEX.ndx -n "$mutant"_"$Factin_or_Gactin"_R312_SD4_TPM_BINDING_INDEX.ndx -nobackup

    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_R312_SD4_TPM_BINDING_INDEX.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/R312_SD4_TPM_BINDING/"$mutant"_"$Factin_or_Gactin"_Chain_C_T318_Chain_C_I327_"$rep".xvg -ref ChainC_T318 -sel ChainC_I327
    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_R312_SD4_TPM_BINDING_INDEX.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/R312_SD4_TPM_BINDING/"$mutant"_"$Factin_or_Gactin"_Chain_C_I317_Chain_C_I327_"$rep".xvg -ref ChainC_I317 -sel ChainC_I327
    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_R312_SD4_TPM_BINDING_INDEX.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/R312_SD4_TPM_BINDING/"$mutant"_"$Factin_or_Gactin"_Chain_C_Q314_Chain_C_I327_"$rep".xvg -ref ChainC_Q314 -sel ChainC_I327
    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_R312_SD4_TPM_BINDING_INDEX.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/R312_SD4_TPM_BINDING/"$mutant"_"$Factin_or_Gactin"_Chain_C_Q314_Chain_C_I329_"$rep".xvg -ref ChainC_Q314 -sel ChainC_I329
    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_R312_SD4_TPM_BINDING_INDEX.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/R312_SD4_TPM_BINDING/"$mutant"_"$Factin_or_Gactin"_Chain_C_M313_Chain_C_I329_"$rep".xvg -ref ChainC_M313 -sel ChainC_I329
    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_R312_SD4_TPM_BINDING_INDEX.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/R312_SD4_TPM_BINDING/"$mutant"_"$Factin_or_Gactin"_Chain_C_A310_Chain_C_I329_"$rep".xvg -ref ChainC_A310 -sel ChainC_I329
    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_R312_SD4_TPM_BINDING_INDEX.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/R312_SD4_TPM_BINDING/"$mutant"_"$Factin_or_Gactin"_Chain_C_A310_Chain_C_A331_"$rep".xvg -ref ChainC_A310 -sel ChainC_A331
    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_R312_SD4_TPM_BINDING_INDEX.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/R312_SD4_TPM_BINDING/"$mutant"_"$Factin_or_Gactin"_Chain_C_M305_Chain_C_R335_"$rep".xvg -ref ChainC_M305 -sel ChainC_R335
    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_R312_SD4_TPM_BINDING_INDEX.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/R312_SD4_TPM_BINDING/"$mutant"_"$Factin_or_Gactin"_Chain_C_M305_Chain_C_K336_"$rep".xvg -ref ChainC_M305 -sel ChainC_K336



    mkdir DISTANCES2/R312_NEARBY_INTERACTIONS/

    ChainC_R312=$((312 + (number_of_residues * (target_chain_number - 1))))
    ChainC_K315=$((315 + (number_of_residues * (target_chain_number - 1))))
    ChainC_E316=$((316 + (number_of_residues * (target_chain_number - 1))))
    ChainC_F262=$((262 + (number_of_residues * (target_chain_number - 1))))
    ChainC_E259=$((259 + (number_of_residues * (target_chain_number - 1))))
    ChainC_F223=$((223 + (number_of_residues * (target_chain_number - 1))))
    ChainC_L221=$((221 + (number_of_residues * (target_chain_number - 1))))
    ChainC_V219=$((219 + (number_of_residues * (target_chain_number - 1))))
    ChainC_G308=$((308 + (number_of_residues * (target_chain_number - 1))))
    ChainC_D222=$((222 + (number_of_residues * (target_chain_number - 1))))
    ChainC_A220=$((220 + (number_of_residues * (target_chain_number - 1))))

    echo -e 'ri '$ChainC_R312'\n name '$start_index 'ChainC_R312\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_R312_NEARBY_INTERACTIONS_INDEX.ndx -nobackup
    echo -e 'ri '$ChainC_K315'\n name '$((start_index + 1)) '$ChainC_K315\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_R312_NEARBY_INTERACTIONS_INDEX.ndx -n "$mutant"_"$Factin_or_Gactin"_R312_NEARBY_INTERACTIONS_INDEX.ndx -nobackup
    echo -e 'ri '$ChainC_E316'\n name '$((start_index + 2)) '$ChainC_E316\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_R312_NEARBY_INTERACTIONS_INDEX.ndx -n "$mutant"_"$Factin_or_Gactin"_R312_NEARBY_INTERACTIONS_INDEX.ndx -nobackup
    echo -e 'ri '$ChainC_F262'\n name '$((start_index + 3)) '$ChainC_F262\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_R312_NEARBY_INTERACTIONS_INDEX.ndx -n "$mutant"_"$Factin_or_Gactin"_R312_NEARBY_INTERACTIONS_INDEX.ndx -nobackup
    echo -e 'ri '$ChainC_E259'\n name '$((start_index + 4)) '$ChainC_E259\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_R312_NEARBY_INTERACTIONS_INDEX.ndx -n "$mutant"_"$Factin_or_Gactin"_R312_NEARBY_INTERACTIONS_INDEX.ndx -nobackup
    echo -e 'ri '$ChainC_F223'\n name '$((start_index + 5)) '$ChainC_F223\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_R312_NEARBY_INTERACTIONS_INDEX.ndx -n "$mutant"_"$Factin_or_Gactin"_R312_NEARBY_INTERACTIONS_INDEX.ndx -nobackup
    echo -e 'ri '$ChainC_L221'\n name '$((start_index + 6)) '$ChainC_L221\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_R312_NEARBY_INTERACTIONS_INDEX.ndx -n "$mutant"_"$Factin_or_Gactin"_R312_NEARBY_INTERACTIONS_INDEX.ndx -nobackup
    echo -e 'ri '$ChainC_V219'\n name '$((start_index + 7)) '$ChainC_V219\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_R312_NEARBY_INTERACTIONS_INDEX.ndx -n "$mutant"_"$Factin_or_Gactin"_R312_NEARBY_INTERACTIONS_INDEX.ndx -nobackup
    echo -e 'ri '$ChainC_G308'\n name '$((start_index + 8)) '$ChainC_G308\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_R312_NEARBY_INTERACTIONS_INDEX.ndx -n "$mutant"_"$Factin_or_Gactin"_R312_NEARBY_INTERACTIONS_INDEX.ndx -nobackup
    echo -e 'ri '$ChainC_D222'\n name '$((start_index + 9)) '$ChainC_D222\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_R312_NEARBY_INTERACTIONS_INDEX.ndx -n "$mutant"_"$Factin_or_Gactin"_R312_NEARBY_INTERACTIONS_INDEX.ndx -nobackup
    echo -e 'ri '$ChainC_A220'\n name '$((start_index + 10)) '$ChainC_A220\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_R312_NEARBY_INTERACTIONS_INDEX.ndx -n "$mutant"_"$Factin_or_Gactin"_R312_NEARBY_INTERACTIONS_INDEX.ndx -nobackup

    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_R312_NEARBY_INTERACTIONS_INDEX.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/R312_NEARBY_INTERACTIONS/"$mutant"_"$Factin_or_Gactin"_Chain_C_R312_Chain_C_K315_"$rep".xvg -ref ChainC_R312 -sel ChainC_K315
    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_R312_NEARBY_INTERACTIONS_INDEX.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/R312_NEARBY_INTERACTIONS/"$mutant"_"$Factin_or_Gactin"_Chain_C_R312_Chain_C_E316_"$rep".xvg -ref ChainC_R312 -sel ChainC_E316
    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_R312_NEARBY_INTERACTIONS_INDEX.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/R312_NEARBY_INTERACTIONS/"$mutant"_"$Factin_or_Gactin"_Chain_C_R312_Chain_C_F262_"$rep".xvg -ref ChainC_R312 -sel ChainC_F262
    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_R312_NEARBY_INTERACTIONS_INDEX.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/R312_NEARBY_INTERACTIONS/"$mutant"_"$Factin_or_Gactin"_Chain_C_R312_Chain_C_E259_"$rep".xvg -ref ChainC_R312 -sel ChainC_E259
    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_R312_NEARBY_INTERACTIONS_INDEX.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/R312_NEARBY_INTERACTIONS/"$mutant"_"$Factin_or_Gactin"_Chain_C_R312_Chain_C_F223_"$rep".xvg -ref ChainC_R312 -sel ChainC_F223
    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_R312_NEARBY_INTERACTIONS_INDEX.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/R312_NEARBY_INTERACTIONS/"$mutant"_"$Factin_or_Gactin"_Chain_C_R312_Chain_C_L221_"$rep".xvg -ref ChainC_R312 -sel ChainC_L221
    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_R312_NEARBY_INTERACTIONS_INDEX.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/R312_NEARBY_INTERACTIONS/"$mutant"_"$Factin_or_Gactin"_Chain_C_R312_Chain_C_V219_"$rep".xvg -ref ChainC_R312 -sel ChainC_V219
    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_R312_NEARBY_INTERACTIONS_INDEX.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/R312_NEARBY_INTERACTIONS/"$mutant"_"$Factin_or_Gactin"_Chain_C_R312_Chain_C_G308_"$rep".xvg -ref ChainC_R312 -sel ChainC_G308
    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_R312_NEARBY_INTERACTIONS_INDEX.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/R312_NEARBY_INTERACTIONS/"$mutant"_"$Factin_or_Gactin"_Chain_C_R312_Chain_C_D222_"$rep".xvg -ref ChainC_R312 -sel ChainC_D222
    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_R312_NEARBY_INTERACTIONS_INDEX.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/R312_NEARBY_INTERACTIONS/"$mutant"_"$Factin_or_Gactin"_Chain_C_R312_Chain_C_A220_"$rep".xvg -ref ChainC_R312 -sel ChainC_A220
    

    mkdir DISTANCES2/R312_A230_HELIX/

    ChainC_K315=$((315 + (number_of_residues * (target_chain_number - 1))))
    ChainC_L221=$((221 + (number_of_residues * (target_chain_number - 1))))
    ChainC_I309=$((309 + (number_of_residues * (target_chain_number - 1))))
    ChainC_V219=$((219 + (number_of_residues * (target_chain_number - 1))))
    ChainC_G308=$((308 + (number_of_residues * (target_chain_number - 1))))
    ChainC_P307=$((307 + (number_of_residues * (target_chain_number - 1))))
    ChainC_Y218=$((218 + (number_of_residues * (target_chain_number - 1))))
    ChainC_Y306=$((306 + (number_of_residues * (target_chain_number - 1))))
    ChainC_C217=$((217 + (number_of_residues * (target_chain_number - 1))))
    ChainC_E259=$((259 + (number_of_residues * (target_chain_number - 1))))
    ChainC_F223=$((223 + (number_of_residues * (target_chain_number - 1))))

    echo -e 'ri '$ChainC_K315'\n name '$start_index 'ChainC_K315\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_R312_A230_HELIX_INDEX.ndx -nobackup
    echo -e 'ri '$ChainC_L221'\n name '$((start_index + 1)) '$ChainC_L221\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_R312_A230_HELIX_INDEX.ndx -n "$mutant"_"$Factin_or_Gactin"_R312_A230_HELIX_INDEX.ndx -nobackup
    echo -e 'ri '$ChainC_I309'\n name '$((start_index + 2)) '$ChainC_I309\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_R312_A230_HELIX_INDEX.ndx -n "$mutant"_"$Factin_or_Gactin"_R312_A230_HELIX_INDEX.ndx -nobackup
    echo -e 'ri '$ChainC_V219'\n name '$((start_index + 3)) '$ChainC_V219\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_R312_A230_HELIX_INDEX.ndx -n "$mutant"_"$Factin_or_Gactin"_R312_A230_HELIX_INDEX.ndx -nobackup
    echo -e 'ri '$ChainC_G308'\n name '$((start_index + 4)) '$ChainC_G308\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_R312_A230_HELIX_INDEX.ndx -n "$mutant"_"$Factin_or_Gactin"_R312_A230_HELIX_INDEX.ndx -nobackup
    echo -e 'ri '$ChainC_P307'\n name '$((start_index + 5)) '$ChainC_P307\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_R312_A230_HELIX_INDEX.ndx -n "$mutant"_"$Factin_or_Gactin"_R312_A230_HELIX_INDEX.ndx -nobackup
    echo -e 'ri '$ChainC_Y218'\n name '$((start_index + 6)) '$ChainC_Y218\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_R312_A230_HELIX_INDEX.ndx -n "$mutant"_"$Factin_or_Gactin"_R312_A230_HELIX_INDEX.ndx -nobackup
    echo -e 'ri '$ChainC_Y306'\n name '$((start_index + 7)) '$ChainC_Y306\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_R312_A230_HELIX_INDEX.ndx -n "$mutant"_"$Factin_or_Gactin"_R312_A230_HELIX_INDEX.ndx -nobackup
    echo -e 'ri '$ChainC_C217'\n name '$((start_index + 8)) '$ChainC_C217\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_R312_A230_HELIX_INDEX.ndx -n "$mutant"_"$Factin_or_Gactin"_R312_A230_HELIX_INDEX.ndx -nobackup
    echo -e 'ri '$ChainC_E259'\n name '$((start_index + 9)) '$ChainC_E259\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_R312_A230_HELIX_INDEX.ndx -n "$mutant"_"$Factin_or_Gactin"_R312_A230_HELIX_INDEX.ndx -nobackup
    echo -e 'ri '$ChainC_F223'\n name '$((start_index + 10)) '$ChainC_F223\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_R312_A230_HELIX_INDEX.ndx -n "$mutant"_"$Factin_or_Gactin"_R312_A230_HELIX_INDEX.ndx -nobackup

    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_R312_A230_HELIX_INDEX.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/R312_A230_HELIX/"$mutant"_"$Factin_or_Gactin"_Chain_C_K315_Chain_C_L221_"$rep".xvg -ref ChainC_K315 -sel ChainC_L221
    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_R312_A230_HELIX_INDEX.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/R312_A230_HELIX/"$mutant"_"$Factin_or_Gactin"_Chain_C_I309_Chain_C_V219_"$rep".xvg -ref ChainC_I309 -sel ChainC_V219
    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_R312_A230_HELIX_INDEX.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/R312_A230_HELIX/"$mutant"_"$Factin_or_Gactin"_Chain_C_G308_Chain_C_V219_"$rep".xvg -ref ChainC_G308 -sel ChainC_V219
    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_R312_A230_HELIX_INDEX.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/R312_A230_HELIX/"$mutant"_"$Factin_or_Gactin"_Chain_C_P307_Chain_C_V219_"$rep".xvg -ref ChainC_P307 -sel ChainC_V219
    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_R312_A230_HELIX_INDEX.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/R312_A230_HELIX/"$mutant"_"$Factin_or_Gactin"_Chain_C_P307_Chain_C_Y218_"$rep".xvg -ref ChainC_P307 -sel ChainC_Y218
    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_R312_A230_HELIX_INDEX.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/R312_A230_HELIX/"$mutant"_"$Factin_or_Gactin"_Chain_C_Y306_Chain_C_V219_"$rep".xvg -ref ChainC_Y306 -sel ChainC_V219
    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_R312_A230_HELIX_INDEX.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/R312_A230_HELIX/"$mutant"_"$Factin_or_Gactin"_Chain_C_Y306_Chain_C_C217_"$rep".xvg -ref ChainC_Y306 -sel ChainC_C217
    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_R312_A230_HELIX_INDEX.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/R312_A230_HELIX/"$mutant"_"$Factin_or_Gactin"_Chain_C_E259_Chain_C_F223_"$rep".xvg -ref ChainC_E259 -sel ChainC_F223


    mkdir DISTANCES2/R312_E270_LOOP/

    ChainC_R312=$((312 + (number_of_residues * (target_chain_number - 1))))
    ChainC_F262=$((262 + (number_of_residues * (target_chain_number - 1))))
    ChainC_E316=$((316 + (number_of_residues * (target_chain_number - 1))))
    ChainC_H275=$((275 + (number_of_residues * (target_chain_number - 1))))
    ChainC_I267=$((267 + (number_of_residues * (target_chain_number - 1))))
    ChainC_D184=$((184 + (number_of_residues * (target_chain_number - 1))))
    ChainC_Y188=$((188 + (number_of_residues * (target_chain_number - 1))))
    ChainC_L180=$((180 + (number_of_residues * (target_chain_number - 1))))
    ChainC_T260=$((260 + (number_of_residues * (target_chain_number - 1))))
    ChainC_L178=$((178 + (number_of_residues * (target_chain_number - 1))))
    ChainC_S271=$((271 + (number_of_residues * (target_chain_number - 1))))
    ChainC_P264=$((264 + (number_of_residues * (target_chain_number - 1))))
    ChainC_L261=$((261 + (number_of_residues * (target_chain_number - 1))))
    ChainC_M269=$((269 + (number_of_residues * (target_chain_number - 1))))
    ChainC_F266=$((266 + (number_of_residues * (target_chain_number - 1))))

    echo -e 'ri '$ChainC_R312'\n name '$start_index 'ChainC_R312\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_R312_E270_LOOP_INDEX.ndx -nobackup
    echo -e 'ri '$ChainC_F262'\n name '$((start_index + 1)) '$ChainC_F262\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_R312_E270_LOOP_INDEX.ndx -n "$mutant"_"$Factin_or_Gactin"_R312_E270_LOOP_INDEX.ndx -nobackup
    echo -e 'ri '$ChainC_E316'\n name '$((start_index + 2)) '$ChainC_E316\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_R312_E270_LOOP_INDEX.ndx -n "$mutant"_"$Factin_or_Gactin"_R312_E270_LOOP_INDEX.ndx -nobackup
    echo -e 'ri '$ChainC_H275'\n name '$((start_index + 3)) '$ChainC_H275\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_R312_E270_LOOP_INDEX.ndx -n "$mutant"_"$Factin_or_Gactin"_R312_E270_LOOP_INDEX.ndx -nobackup
    echo -e 'ri '$ChainC_I267'\n name '$((start_index + 4)) '$ChainC_I267\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_R312_E270_LOOP_INDEX.ndx -n "$mutant"_"$Factin_or_Gactin"_R312_E270_LOOP_INDEX.ndx -nobackup
    echo -e 'ri '$ChainC_D184'\n name '$((start_index + 5)) '$ChainC_D184\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_R312_E270_LOOP_INDEX.ndx -n "$mutant"_"$Factin_or_Gactin"_R312_E270_LOOP_INDEX.ndx -nobackup
    echo -e 'ri '$ChainC_Y188'\n name '$((start_index + 6)) '$ChainC_Y188\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_R312_E270_LOOP_INDEX.ndx -n "$mutant"_"$Factin_or_Gactin"_R312_E270_LOOP_INDEX.ndx -nobackup
    echo -e 'ri '$ChainC_L180'\n name '$((start_index + 7)) '$ChainC_L180\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_R312_E270_LOOP_INDEX.ndx -n "$mutant"_"$Factin_or_Gactin"_R312_E270_LOOP_INDEX.ndx -nobackup
    echo -e 'ri '$ChainC_T260'\n name '$((start_index + 8)) '$ChainC_T260\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_R312_E270_LOOP_INDEX.ndx -n "$mutant"_"$Factin_or_Gactin"_R312_E270_LOOP_INDEX.ndx -nobackup
    echo -e 'ri '$ChainC_L178'\n name '$((start_index + 9)) '$ChainC_L178\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_R312_E270_LOOP_INDEX.ndx -n "$mutant"_"$Factin_or_Gactin"_R312_E270_LOOP_INDEX.ndx -nobackup
    echo -e 'ri '$ChainC_S271'\n name '$((start_index + 10)) '$ChainC_S271\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_R312_E270_LOOP_INDEX.ndx -n "$mutant"_"$Factin_or_Gactin"_R312_E270_LOOP_INDEX.ndx -nobackup
    echo -e 'ri '$ChainC_P264'\n name '$((start_index + 11)) '$ChainC_P264\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_R312_E270_LOOP_INDEX.ndx -n "$mutant"_"$Factin_or_Gactin"_R312_E270_LOOP_INDEX.ndx -nobackup
    echo -e 'ri '$ChainC_L261'\n name '$((start_index + 12)) '$ChainC_L261\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_R312_E270_LOOP_INDEX.ndx -n "$mutant"_"$Factin_or_Gactin"_R312_E270_LOOP_INDEX.ndx -nobackup
    echo -e 'ri '$ChainC_M269'\n name '$((start_index + 13)) '$ChainC_M269\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_R312_E270_LOOP_INDEX.ndx -n "$mutant"_"$Factin_or_Gactin"_R312_E270_LOOP_INDEX.ndx -nobackup
    echo -e 'ri '$ChainC_F266'\n name '$((start_index + 14)) '$ChainC_F266\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_R312_E270_LOOP_INDEX.ndx -n "$mutant"_"$Factin_or_Gactin"_R312_E270_LOOP_INDEX.ndx -nobackup

    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_R312_E270_LOOP_INDEX.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/R312_E270_LOOP/"$mutant"_"$Factin_or_Gactin"_Chain_C_R312_Chain_C_F262_"$rep".xvg -ref ChainC_R312 -sel ChainC_F262
    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_R312_E270_LOOP_INDEX.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/R312_E270_LOOP/"$mutant"_"$Factin_or_Gactin"_Chain_C_E316_Chain_C_H275_"$rep".xvg -ref ChainC_E316 -sel ChainC_H275
    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_R312_E270_LOOP_INDEX.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/R312_E270_LOOP/"$mutant"_"$Factin_or_Gactin"_Chain_C_F262_Chain_C_H275_"$rep".xvg -ref ChainC_F262 -sel ChainC_H275
    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_R312_E270_LOOP_INDEX.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/R312_E270_LOOP/"$mutant"_"$Factin_or_Gactin"_Chain_C_I267_Chain_C_D184_"$rep".xvg -ref ChainC_I267 -sel ChainC_D184
    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_R312_E270_LOOP_INDEX.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/R312_E270_LOOP/"$mutant"_"$Factin_or_Gactin"_Chain_C_I267_Chain_C_Y188_"$rep".xvg -ref ChainC_I267 -sel ChainC_Y188
    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_R312_E270_LOOP_INDEX.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/R312_E270_LOOP/"$mutant"_"$Factin_or_Gactin"_Chain_C_I267_Chain_C_L180_"$rep".xvg -ref ChainC_I267 -sel ChainC_L180
    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_R312_E270_LOOP_INDEX.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/R312_E270_LOOP/"$mutant"_"$Factin_or_Gactin"_Chain_C_I267_Chain_C_T260_"$rep".xvg -ref ChainC_I267 -sel ChainC_T260
    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_R312_E270_LOOP_INDEX.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/R312_E270_LOOP/"$mutant"_"$Factin_or_Gactin"_Chain_C_L178_Chain_C_S271_"$rep".xvg -ref ChainC_L178 -sel ChainC_S271
    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_R312_E270_LOOP_INDEX.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/R312_E270_LOOP/"$mutant"_"$Factin_or_Gactin"_Chain_C_L178_Chain_C_P264_"$rep".xvg -ref ChainC_L178 -sel ChainC_P264
    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_R312_E270_LOOP_INDEX.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/R312_E270_LOOP/"$mutant"_"$Factin_or_Gactin"_Chain_C_L180_Chain_C_L261_"$rep".xvg -ref ChainC_L180 -sel ChainC_L261
    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_R312_E270_LOOP_INDEX.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/R312_E270_LOOP/"$mutant"_"$Factin_or_Gactin"_Chain_C_L180_Chain_C_T260_"$rep".xvg -ref ChainC_L180 -sel ChainC_T260
    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_R312_E270_LOOP_INDEX.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/R312_E270_LOOP/"$mutant"_"$Factin_or_Gactin"_Chain_C_L180_Chain_C_I267_"$rep".xvg -ref ChainC_L180 -sel ChainC_I267
    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_R312_E270_LOOP_INDEX.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/R312_E270_LOOP/"$mutant"_"$Factin_or_Gactin"_Chain_C_L180_Chain_C_M269_"$rep".xvg -ref ChainC_L180 -sel ChainC_M269
    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_R312_E270_LOOP_INDEX.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/R312_E270_LOOP/"$mutant"_"$Factin_or_Gactin"_Chain_C_F266_Chain_C_T260_"$rep".xvg -ref ChainC_F266 -sel ChainC_T260

    mkdir DISTANCES2/R312_CROSS_DOMAIN/

    ChainC_L180=$((180 + (number_of_residues * (target_chain_number - 1))))
    ChainC_T160=$((160 + (number_of_residues * (target_chain_number - 1))))
    ChainC_A181=$((181 + (number_of_residues * (target_chain_number - 1))))
    ChainC_G158=$((158 + (number_of_residues * (target_chain_number - 1))))
    ChainC_D157=$((157 + (number_of_residues * (target_chain_number - 1))))
    ChainC_G182=$((182 + (number_of_residues * (target_chain_number - 1))))
    ChainC_R183=$((183 + (number_of_residues * (target_chain_number - 1))))
    ChainC_G156=$((156 + (number_of_residues * (target_chain_number - 1))))
    ChainC_T203=$((203 + (number_of_residues * (target_chain_number - 1))))
    ChainC_L67=$((67 + (number_of_residues * (target_chain_number - 1))))
    ChainC_T66=$((66 + (number_of_residues * (target_chain_number - 1))))
    
    echo -e 'ri '$ChainC_L180'\n name '$start_index 'ChainC_L180\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_R312_CROSS_DOMAIN_INDEX.ndx -nobackup
    echo -e 'ri '$ChainC_T160'\n name '$((start_index + 1)) '$ChainC_T160\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_R312_CROSS_DOMAIN_INDEX.ndx -n "$mutant"_"$Factin_or_Gactin"_R312_CROSS_DOMAIN_INDEX.ndx -nobackup
    echo -e 'ri '$ChainC_A181'\n name '$((start_index + 2)) '$ChainC_A181\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_R312_CROSS_DOMAIN_INDEX.ndx -n "$mutant"_"$Factin_or_Gactin"_R312_CROSS_DOMAIN_INDEX.ndx -nobackup
    echo -e 'ri '$ChainC_G158'\n name '$((start_index + 3)) '$ChainC_G158\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_R312_CROSS_DOMAIN_INDEX.ndx -n "$mutant"_"$Factin_or_Gactin"_R312_CROSS_DOMAIN_INDEX.ndx -nobackup
    echo -e 'ri '$ChainC_D157'\n name '$((start_index + 4)) '$ChainC_D157\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_R312_CROSS_DOMAIN_INDEX.ndx -n "$mutant"_"$Factin_or_Gactin"_R312_CROSS_DOMAIN_INDEX.ndx -nobackup
    echo -e 'ri '$ChainC_G182'\n name '$((start_index + 5)) '$ChainC_G182\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_R312_CROSS_DOMAIN_INDEX.ndx -n "$mutant"_"$Factin_or_Gactin"_R312_CROSS_DOMAIN_INDEX.ndx -nobackup
    echo -e 'ri '$ChainC_R183'\n name '$((start_index + 6)) '$ChainC_R183\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_R312_CROSS_DOMAIN_INDEX.ndx -n "$mutant"_"$Factin_or_Gactin"_R312_CROSS_DOMAIN_INDEX.ndx -nobackup
    echo -e 'ri '$ChainC_G156'\n name '$((start_index + 7)) '$ChainC_G156\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_R312_CROSS_DOMAIN_INDEX.ndx -n "$mutant"_"$Factin_or_Gactin"_R312_CROSS_DOMAIN_INDEX.ndx -nobackup
    echo -e 'ri '$ChainC_T203'\n name '$((start_index + 8)) '$ChainC_T203\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_R312_CROSS_DOMAIN_INDEX.ndx -n "$mutant"_"$Factin_or_Gactin"_R312_CROSS_DOMAIN_INDEX.ndx -nobackup
    echo -e 'ri '$ChainC_L67'\n name '$((start_index + 9)) '$ChainC_L67\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_R312_CROSS_DOMAIN_INDEX.ndx -n "$mutant"_"$Factin_or_Gactin"_R312_CROSS_DOMAIN_INDEX.ndx -nobackup
    echo -e 'ri '$ChainC_T66'\n name '$((start_index + 10)) '$ChainC_T66\n q' | gmx make_ndx -f PDBs/"$mutant"_"$Factin_or_Gactin"_FIRST_FRAME.pdb -o "$mutant"_"$Factin_or_Gactin"_R312_CROSS_DOMAIN_INDEX.ndx -n "$mutant"_"$Factin_or_Gactin"_R312_CROSS_DOMAIN_INDEX.ndx -nobackup

    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_R312_CROSS_DOMAIN_INDEX.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/R312_CROSS_DOMAIN/"$mutant"_"$Factin_or_Gactin"_Chain_C_L180_Chain_C_T160_"$rep".xvg -ref ChainC_L180 -sel ChainC_T160
    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_R312_CROSS_DOMAIN_INDEX.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/R312_CROSS_DOMAIN/"$mutant"_"$Factin_or_Gactin"_Chain_C_A181_Chain_C_G158_"$rep".xvg -ref ChainC_A181 -sel ChainC_G158
    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_R312_CROSS_DOMAIN_INDEX.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/R312_CROSS_DOMAIN/"$mutant"_"$Factin_or_Gactin"_Chain_C_D157_Chain_C_A181_"$rep".xvg -ref ChainC_D157 -sel ChainC_A181
    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_R312_CROSS_DOMAIN_INDEX.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/R312_CROSS_DOMAIN/"$mutant"_"$Factin_or_Gactin"_Chain_C_D157_Chain_C_G182_"$rep".xvg -ref ChainC_D157 -sel ChainC_G182
    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_R312_CROSS_DOMAIN_INDEX.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/R312_CROSS_DOMAIN/"$mutant"_"$Factin_or_Gactin"_Chain_C_D157_Chain_C_R183_"$rep".xvg -ref ChainC_D157 -sel ChainC_R183
    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_R312_CROSS_DOMAIN_INDEX.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/R312_CROSS_DOMAIN/"$mutant"_"$Factin_or_Gactin"_Chain_C_D157_Chain_C_G156_"$rep".xvg -ref ChainC_D157 -sel ChainC_G156
    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_R312_CROSS_DOMAIN_INDEX.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/R312_CROSS_DOMAIN/"$mutant"_"$Factin_or_Gactin"_Chain_C_T203_Chain_C_L67_"$rep".xvg -ref ChainC_T203 -sel ChainC_L67
    gmx pairdist -n "$mutant"_"$Factin_or_Gactin"_R312_CROSS_DOMAIN_INDEX.ndx -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES2/R312_CROSS_DOMAIN/"$mutant"_"$Factin_or_Gactin"_Chain_C_T203_Chain_C_T66_"$rep".xvg -ref ChainC_T203 -sel ChainC_T66


    echo -e 'Backbone\n Backbone_&_Chain_'$target_chain'_protein' | gmx rms -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -s "$mutant"_"$Factin_or_Gactin"_md_continue.tpr -n ../"$mutant"_"$Factin_or_Gactin"_analysis_index.ndx -o "$mutant"_Chain_"$target_chain"_rmsd_backbone_FIT_"$rep".xvg -fit rot+trans
    echo -e 'Backbone\n Backbone_&_Chain_'$target_chain'_protein' | gmx rms -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -s "$mutant"_"$Factin_or_Gactin"_md_continue.tpr -n ../"$mutant"_"$Factin_or_Gactin"_analysis_index.ndx -o "$mutant"_Chain_"$target_chain"_rmsd_backbone_FIT_PREV_"$rep".xvg -fit rot+trans -prev 1
    echo -e 'Backbone\n Backbone_&_Chain_'$target_chain'_protein' | gmx rms -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -s "$mutant"_"$Factin_or_Gactin"_md_continue.tpr -n ../"$mutant"_"$Factin_or_Gactin"_analysis_index.ndx -o "$mutant"_Chain_"$target_chain"_rmsd_backbone_PREV_"$rep".xvg -prev 1 
    
    echo -e 'Backbone\n Backbone_&_Chain_'$target_chain'_SD1' | gmx rms -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -s "$mutant"_"$Factin_or_Gactin"_md_continue.tpr -n ../"$mutant"_"$Factin_or_Gactin"_analysis_index.ndx -o "$mutant"_Chain_"$target_chain"_rmsd_backbone_FIT_SD1_"$rep".xvg -fit rot+trans
    echo -e 'Backbone\n Backbone_&_Chain_'$target_chain'_SD1' | gmx rms -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -s "$mutant"_"$Factin_or_Gactin"_md_continue.tpr -n ../"$mutant"_"$Factin_or_Gactin"_analysis_index.ndx -o "$mutant"_Chain_"$target_chain"_rmsd_backbone_FIT_PREV_SD1_"$rep".xvg -fit rot+trans -prev 1
    echo -e 'Backbone\n Backbone_&_Chain_'$target_chain'_SD1' | gmx rms -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -s "$mutant"_"$Factin_or_Gactin"_md_continue.tpr -n ../"$mutant"_"$Factin_or_Gactin"_analysis_index.ndx -o "$mutant"_Chain_"$target_chain"_rmsd_backbone_PREV_SD1_"$rep".xvg -prev 1 
    
    echo -e 'Backbone\n Backbone_&_Chain_'$target_chain'_SD2' | gmx rms -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -s "$mutant"_"$Factin_or_Gactin"_md_continue.tpr -n ../"$mutant"_"$Factin_or_Gactin"_analysis_index.ndx -o "$mutant"_Chain_"$target_chain"_rmsd_backbone_FIT_SD2_"$rep".xvg -fit rot+trans
    echo -e 'Backbone\n Backbone_&_Chain_'$target_chain'_SD2' | gmx rms -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -s "$mutant"_"$Factin_or_Gactin"_md_continue.tpr -n ../"$mutant"_"$Factin_or_Gactin"_analysis_index.ndx -o "$mutant"_Chain_"$target_chain"_rmsd_backbone_FIT_PREV_SD2_"$rep".xvg -fit rot+trans -prev 1
    echo -e 'Backbone\n Backbone_&_Chain_'$target_chain'_SD2' | gmx rms -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -s "$mutant"_"$Factin_or_Gactin"_md_continue.tpr -n ../"$mutant"_"$Factin_or_Gactin"_analysis_index.ndx -o "$mutant"_Chain_"$target_chain"_rmsd_backbone_PREV_SD2_"$rep".xvg -prev 1
    
    echo -e 'Backbone\n Backbone_&_Chain_'$target_chain'_SD3' | gmx rms -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -s "$mutant"_"$Factin_or_Gactin"_md_continue.tpr -n ../"$mutant"_"$Factin_or_Gactin"_analysis_index.ndx -o "$mutant"_Chain_"$target_chain"_rmsd_backbone_FIT_SD3_"$rep".xvg -fit rot+trans
    echo -e 'Backbone\n Backbone_&_Chain_'$target_chain'_SD3' | gmx rms -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -s "$mutant"_"$Factin_or_Gactin"_md_continue.tpr -n ../"$mutant"_"$Factin_or_Gactin"_analysis_index.ndx -o "$mutant"_Chain_"$target_chain"_rmsd_backbone_FIT_PREV_SD3_"$rep".xvg -fit rot+trans -prev 1
    echo -e 'Backbone\n Backbone_&_Chain_'$target_chain'_SD3' | gmx rms -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -s "$mutant"_"$Factin_or_Gactin"_md_continue.tpr -n ../"$mutant"_"$Factin_or_Gactin"_analysis_index.ndx -o "$mutant"_Chain_"$target_chain"_rmsd_backbone_PREV_SD3_"$rep".xvg -prev 1
    
    echo -e 'Backbone\n Backbone_&_Chain_'$target_chain'_SD4' | gmx rms -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -s "$mutant"_"$Factin_or_Gactin"_md_continue.tpr -n ../"$mutant"_"$Factin_or_Gactin"_analysis_index.ndx -o "$mutant"_Chain_"$target_chain"_rmsd_backbone_FIT_SD4_"$rep".xvg -fit rot+trans
    echo -e 'Backbone\n Backbone_&_Chain_'$target_chain'_SD4' | gmx rms -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -s "$mutant"_"$Factin_or_Gactin"_md_continue.tpr -n ../"$mutant"_"$Factin_or_Gactin"_analysis_index.ndx -o "$mutant"_Chain_"$target_chain"_rmsd_backbone_FIT_PREV_SD4_"$rep".xvg -fit rot+trans -prev 1
    echo -e 'Backbone\n Backbone_&_Chain_'$target_chain'_SD4' | gmx rms -f TRAJECTORIES/"$mutant"_"$Factin_or_Gactin"_md_trajec_fit.xtc -s "$mutant"_"$Factin_or_Gactin"_md_continue.tpr -n ../"$mutant"_"$Factin_or_Gactin"_analysis_index.ndx -o "$mutant"_Chain_"$target_chain"_rmsd_backbone_PREV_SD4_"$rep".xvg -prev 1 
  done
done
