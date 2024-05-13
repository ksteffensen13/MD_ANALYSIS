#!/bin/bash
#THIS SCRIPT FIRST GENERATES VARIABLES BASED ON THE COMMAND LINE INPUT, THEN USES THOSE VARIABLES TO GENERATE/RUN A SCRIPT WITH A JOB NAME CORRESPONDING TO THOSE VARIABLES
##FOR EXAMPLE, IF RUNNING ANALYSIS ON R312C F-ACTIN, IT WILL CREATE A SCRIPT THAT RUNS A SLURM JOB NAMED 'R312C_FACTIN_ANALYSIS'
##then, this will run the script for R312C F-actin, and basically every filename will be in the format "R312C_Factin_..."
#######NOTE: THERE IS A SECTION BELOW THAT REQUIRES MANUAL INPUTS
######TO RUN THIS, MAKE IT AN EXECUTABLE WITH CHMOD
######EXAMPLE: chmod +x submit_job.sh
#####THEN: ./submit_job.sh R312C Factin
#to run this file, need gromacs/vmd installed in your home directory. Also need the networkanalysis tutorial files

##FIRST, need to define the mutant so that we can input it in command line and call it in filenames
# Check if the "mutant" argument is provided
if [ -z "$1" ] || [ -z "$2" ]; then
    echo "Usage: $0 <mutant_value> <Factin_or_Gactin>"
    exit 1
fi
# Assign the value of the "mutant" argument to a variable
mutant="$1"
Factin_or_Gactin="$2"

#Here, we store the content of the script in a variable, then we'll call that variable to generate a script that we then run
#we do this because there is no way to include a variable in the SLURM job name, SO we use this to
script_content=$(cat <<EOF
#!/bin/bash
#SBATCH --time=2-00:00:00           # time limit (D-HH:MM:SS)
#SBATCH --account=rrg-jfdawson
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=4
#SBATCH --mem=0
#SBATCH --output=${mutant}_${Factin_or_Gactin}_post_MD.out
#SBATCH --job-name=${mutant}_${Factin_or_Gactin}_post_MD
#SBATCH --mail-user=ksteffen@uoguelph.ca
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=ALL

module load StdEnv/2020 gcc/9.3.0 cuda/11.0 openmpi/4.0.3 gromacs/2021.2
export GMXLIB=/home/ksteffen/gromacs-2021.5/share/top
export OMP_NUM_THREADS="\${SLURM_CPUS_PER_TASK:-1}"

###################################ENTER MANUALLY##################################################################################

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
simulation_length=200
# check index files from the equilibration step to get the number of indices
# THIS DOES NOT INCLUDE THE INDICES YOU MADE LIKE "PROTEIN_MG_ADP", JUST THE DEFAULT ONES
#to figure this out, take your STARTING PDB, and use gmx make_ndx -f STARTING_PDB.pdb -o TEST.ndx and look at the last index number, then cancel the command
number_of_indices=19
#set time you want to consider the simulation at equilibrium (usually 70% of simulation time)
eq_start_time=70



#BELOW HERE IS AUTOMATED AND ONLY REQUIRES CHANGES IF YOURE USING A VERY DIFFERENT SYSTEM

##FIRST, need to define the mutant so that we can input it in command line and call it in filenames
# Check if the "mutant" argument is provided
if [ -z "$1" ] || [ -z "$2" ]; then
    echo "Usage: $0 <mutant_value> <Factin_or_Gactin>"
    exit 1
fi
# Assign the value of the "mutant" argument to a variable
mutant="$1"
Factin_or_Gactin="$2"

##########################################FOR CALCULATIONS##############################################################

chains=('A' 'B' 'C' 'D' 'E' 'F' 'G' 'H' 'I' 'J' 'K')

# Find the index of the target_chain in the array
# this looks through the array 'chains', at every index of the array, and if that entry is the same as 'target chain',
# then it sets that index number to be the variable 'target_chain_number' which can be called later on
#note: 'A' in chains is index number 0 (numbering starts at 0 not 1), so to get target chain number for calculations have to add 1
for i in "\${!chains[@]}"; do
    if [[ "\${chains[i]}" == "\$target_chain" ]]; then
        target_chain_number=\$((i + 1))
        target_chain_index_number="\$i"
        break
    fi
done

#based on the target chain, calculate the letters and numbers of the neighbouring chains which can be called for analyses later
if [[ "\$target_chain_number" -gt 2 ]]; then
  interacting_chain_number_minus2=\$((target_chain_index_number- 2))
  interacting_chain_minus2="\${chains[\$interacting_chain_number_minus2]}"
  interacting_chain_number_minus1=\$((target_chain_index_number - 1))
  interacting_chain_minus1="\${chains[\$interacting_chain_number_minus1]}"
  interacting_chain_number_plus1=\$((target_chain_index_number + 1))
  interacting_chain_plus1="\${chains[\$interacting_chain_number_plus1]}"
  interacting_chain_number_plus2=\$((target_chain_index_number + 2))
  interacting_chain_plus2="\${chains[\$interacting_chain_number_plus2]}"
elif [[ "\$target_chain_number" -eq 2 ]]; then
  interacting_chain_number_minus1=1
  interacting_chain_minus1='A'
  interacting_chain_number_plus1=3
  interacting_chain_plus1='C'
  interacting_chain_number_plus2=4
  interacting_chain_plus2='D'
elif [[ "\$target_chain_number" -eq 1 ]]; then
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

target_chain_start=\$((target_chain_number * number_of_residues - 374))
target_chain_end=\$((target_chain_number * number_of_residues))
target_chain_residues="\$target_chain_start-\$target_chain_end"

# the nucleotide residue number in the index file comes after all of the protein residues
# So, we calculate the number of the last residue, then add a number corresponding to the target chain
# for example, if we want the ADP of chain C (chain #3), in a 5 protomer filament, then ADP index=(5*375) + 3 = 1878
# this is because all 5 protein chains come first, then each ADP for each chain
nucleotide_index_number=\$((number_of_chains * number_of_residues + target_chain_number))
# similar process for cation. Comes after all ADP, so calculate the number of protein residues, then add the total number of ADP, then add the target chain number
cation_index_number=\$((number_of_chains * number_of_residues + number_of_chains + target_chain_number))
# if ADP-Pi, the PO4 usually comes last
if [[ "\$nucleotide_state" == "ADP-Pi" ]]; then
  po4_index_number=\$((number_of_chains * number_of_residues + 2 * number_of_chains + target_chain_number))
  nucleotide='ADP'
elif [[ "\$nucleotide_state" == "ADP" ]]; then
  nucleotide='ADP'
else
  nucleotide='ATP'
fi

start_index=\$((number_of_indices + 1))

# calculate subdomain residue numbers for the desired chain
subdomain1_residue_start1=\$((1 + (number_of_residues) * (target_chain_number - 1)))
subdomain1_residue_end1=\$((32 + (number_of_residues) * (target_chain_number - 1)))
subdomain1_residue_start2=\$((70 + (number_of_residues) * (target_chain_number - 1)))
subdomain1_residue_end2=\$((137 + (number_of_residues) * (target_chain_number - 1)))
subdomain1_residue_start3=\$((337 + (number_of_residues) * (target_chain_number - 1)))
subdomain1_residue_end3=\$(((number_of_residues) + (number_of_residues) * (target_chain_number - 1)))

subdomain2_residue_start=\$((33 + (number_of_residues) * (target_chain_number - 1)))
subdomain2_residue_end=\$((69 + (number_of_residues) * (target_chain_number - 1)))

subdomain3_residue_start1=\$((138 + (number_of_residues) * (target_chain_number - 1)))
subdomain3_residue_end1=\$((178 + (number_of_residues) * (target_chain_number - 1)))
subdomain3_residue_start2=\$((274 + (number_of_residues) * (target_chain_number - 1)))
subdomain3_residue_end2=\$((336 + (number_of_residues) * (target_chain_number - 1)))

subdomain4_residue_start=\$((179 + (number_of_residues) * (target_chain_number - 1)))
subdomain4_residue_end=\$((273 + (number_of_residues) * (target_chain_number - 1)))

# generate text strings for subdomains that can be called for making index files
target_chain_subdomain1_residues1="\$subdomain1_residue_start1-\$subdomain1_residue_end1"
target_chain_subdomain1_residues2="\$subdomain1_residue_start2-\$subdomain1_residue_end2"
target_chain_subdomain1_residues3="\$subdomain1_residue_start3-\$subdomain1_residue_end3"

target_chain_subdomain2_residues="\$subdomain2_residue_start-\$subdomain2_residue_end"

target_chain_subdomain3_residues1="\$subdomain3_residue_start1-\$subdomain3_residue_end1"
target_chain_subdomain3_residues2="\$subdomain3_residue_start2-\$subdomain3_residue_end2"

target_chain_subdomain4_residues="\$subdomain4_residue_start-\$subdomain4_residue_end"

##########################################ACTUAL SCRIPT################################################################





######PROCESS STRUCTURE/TRAJECTORY
for replicate in Replicate1 Replicate2 Replicate3; do
  cd "\$directory_path"/"\$mutant"/"\$replicate"

  if [ "\$replicate" == "Replicate1" ]; then
      rep="rep1"
  elif [ "\$replicate" == "Replicate2" ]; then
      rep="rep2"
  elif [ "\$replicate" == "Replicate3" ]; then
      rep="rep3"
  fi

  ##MAKE INDEX FILES SO WE CAN EXTRACT TARGET CHAIN FOR SPECIFIC ANALYSES

  if [ "\$nucleotide_state" == "ADP-Pi" ]; then

    #first, make whole
    #NOTE: the order of the indexes 'Protein_MG_ADP_PO4' depends on what you have in your .mdp files
    echo -e 'Protein_'\$cation'_ADP_PO4' | gmx trjconv -s "\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_md.xtc -o "\$mutant"_"\$Factin_or_Gactin"_md_trajec_whole.xtc -n "\$mutant"_"\$Factin_or_Gactin"_index.ndx -pbc whole

    #second, dump initial structure and remove PBC jumps

    echo -e 'Protein_'\$cation'_ADP_PO4' | gmx trjconv -s "\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_whole.xtc -n "\$mutant"_"\$Factin_or_Gactin"_index.ndx -o "\$mutant"_"\$Factin_or_Gactin"_md_trajec_whole.gro -dump 0
    echo -e 'Protein_'\$cation'_ADP_PO4\n Protein_'\$cation'_ADP_PO4' | gmx trjconv -s "\$mutant"_"\$Factin_or_Gactin"_md_trajec_whole.gro -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_whole.xtc -o "\$mutant"_"\$Factin_or_Gactin"_md_trajec_nojump.xtc -n "\$mutant"_"\$Factin_or_Gactin"_index.ndx -center -pbc nojump

    #center

    echo -e 'Protein_'\$cation'_ADP_PO4\n Protein_'\$cation'_ADP_PO4' | gmx trjconv -s "\$mutant"_"\$Factin_or_Gactin"_md_trajec_whole.gro -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_nojump.xtc -o "\$mutant"_"\$Factin_or_Gactin"_md_trajec_center.xtc -n "\$mutant"_"\$Factin_or_Gactin"_index.ndx -center -pbc mol

    #fit

    echo -e 'Protein_'\$cation'_ADP_PO4\n Protein_'\$cation'_ADP_PO4' | gmx trjconv -s "\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_center.xtc -o "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -n "\$mutant"_"\$Factin_or_Gactin"_index.ndx -fit rot+trans
    echo -e 'Protein_'\$cation'_ADP_PO4\n Protein_'\$cation'_ADP_PO4' | gmx trjconv -s "\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_center.xtc -o "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -n "\$mutant"_"\$Factin_or_Gactin"_index.ndx -fit rot+trans -dump 0
    echo -e 'Protein_'\$cation'_ADP_PO4\n Protein_'\$cation'_ADP_PO4' | gmx trjconv -s "\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_center.xtc -o "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.gro -n "\$mutant"_"\$Factin_or_Gactin"_index.ndx -fit rot+trans -dump 0
    echo -e 'Protein_'\$cation'_ADP_PO4\n Protein_'\$cation'_ADP_PO4' | gmx trjconv -s "\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_center.xtc -o "\$mutant"_"\$Factin_or_Gactin"_FINAL_FRAME.pdb -n "\$mutant"_"\$Factin_or_Gactin"_index.ndx -fit rot+trans -tu ns -dump "\$simulation_length"

    echo -e 'ri '\$target_chain_residues'|ri '\$nucleotide_index_number'|ri '\$cation_index_number'|ri '\$po4_index_number'\n name '\$start_index 'Chain_'\$target_chain'\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -nobackup
    echo -e 'ri '\$target_chain_residues'\n name '\$((start_index + 1)) 'Chain_'\$target_chain'_protein\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -n "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -nobackup
    echo -e 'ri '\$target_chain_subdomain1_residues1'|ri '\$target_chain_subdomain1_residues2'|ri '\$target_chain_subdomain1_residues3'\n name '\$((start_index + 2)) 'Chain_'\$target_chain'_SD1\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -n "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -nobackup
    echo -e 'ri '\$target_chain_subdomain2_residues'\n name '\$((start_index + 3)) 'Chain_'\$target_chain'_SD2\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -n "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -nobackup
    echo -e 'ri '\$target_chain_subdomain3_residues1'|ri '\$target_chain_subdomain3_residues2'\n name '\$((start_index + 4)) 'Chain_'\$target_chain'_SD3\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -n "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -nobackup
    echo -e 'ri '\$target_chain_subdomain4_residues'\n name '\$((start_index + 5)) 'Chain_'\$target_chain'_SD4\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -n "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -nobackup
    echo -e '"Backbone" & "Chain_'\$target_chain'_protein"\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -n "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -nobackup
    echo -e '"Backbone" & "Chain_'\$target_chain'_SD1"\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -n "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -nobackup
    echo -e '"Backbone" & "Chain_'\$target_chain'_SD2"\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -n "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -nobackup
    echo -e '"Backbone" & "Chain_'\$target_chain'_SD3"\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -n "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -nobackup
    echo -e '"Backbone" & "Chain_'\$target_chain'_SD4"\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -n "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -nobackup
    echo -e '"C-alpha" & "Chain_'\$target_chain'_protein"\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -n "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -nobackup
    echo -e '"C-alpha" & "Chain_'\$target_chain'_SD1"\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -n "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -nobackup
    echo -e '"C-alpha" & "Chain_'\$target_chain'_SD2"\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -n "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -nobackup
    echo -e '"C-alpha" & "Chain_'\$target_chain'_SD3"\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -n "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -nobackup
    echo -e '"C-alpha" & "Chain_'\$target_chain'_SD4"\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -n "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -nobackup


  else

    #first, make whole
    #NOTE: the order of the indexes 'Protein_'\$cation'_ADP_PO4' depends on what you have in your .mdp files
    echo -e 'Protein_'\$cation'_'\$nucleotide | gmx trjconv -s "\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_md.xtc -o "\$mutant"_"\$Factin_or_Gactin"_md_trajec_whole.xtc -n "\$mutant"_"\$Factin_or_Gactin"_index.ndx -pbc whole

    #second, dump initial structure and remove PBC jumps

    echo -e 'Protein_'\$cation'_'\$nucleotide | gmx trjconv -s "\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_whole.xtc -n "\$mutant"_"\$Factin_or_Gactin"_index.ndx -o "\$mutant"_"\$Factin_or_Gactin"_md_trajec_whole.gro -dump 0
    echo -e 'Protein_'\$cation'_'\$nucleotide'\n Protein_'\$cation'_'\$nucleotide | gmx trjconv -s "\$mutant"_"\$Factin_or_Gactin"_md_trajec_whole.gro -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_whole.xtc -o "\$mutant"_"\$Factin_or_Gactin"_md_trajec_nojump.xtc -n "\$mutant"_"\$Factin_or_Gactin"_index.ndx -center -pbc nojump

    #center

    echo -e 'Protein_'\$cation'_'\$nucleotide'\n Protein_'\$cation'_'\$nucleotide | gmx trjconv -s "\$mutant"_"\$Factin_or_Gactin"_md_trajec_whole.gro -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_nojump.xtc -o "\$mutant"_"\$Factin_or_Gactin"_md_trajec_center.xtc -n "\$mutant"_"\$Factin_or_Gactin"_index.ndx -center -pbc mol

    #fit

    echo -e 'Protein_'\$cation'_'\$nucleotide'\n Protein_'\$cation'_'\$nucleotide | gmx trjconv -s "\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_center.xtc -o "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -n "\$mutant"_"\$Factin_or_Gactin"_index.ndx -fit rot+trans
    echo -e 'Protein_'\$cation'_'\$nucleotide'\n Protein_'\$cation'_'\$nucleotide | gmx trjconv -s "\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_center.xtc -o "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -n "\$mutant"_"\$Factin_or_Gactin"_index.ndx -fit rot+trans -dump 0
    echo -e 'Protein_'\$cation'_'\$nucleotide'\n Protein_'\$cation'_'\$nucleotide | gmx trjconv -s "\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_center.xtc -o "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.gro -n "\$mutant"_"\$Factin_or_Gactin"_index.ndx -fit rot+trans -dump 0
    echo -e 'Protein_'\$cation'_'\$nucleotide'\n Protein_'\$cation'_'\$nucleotide | gmx trjconv -s "\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_center.xtc -o "\$mutant"_"\$Factin_or_Gactin"_FINAL_FRAME.pdb -n "\$mutant"_"\$Factin_or_Gactin"_index.ndx -fit rot+trans -tu ns -dump "\$simulation_length"

    echo -e 'ri '\$target_chain_residues'|ri '\$nucleotide_index_number'|ri '\$cation_index_number'\n name '\$start_index 'Chain_'\$target_chain'\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -nobackup
    echo -e 'ri '\$target_chain_residues'\n name '\$((start_index + 1)) 'Chain_'\$target_chain'_protein\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -n "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -nobackup
    echo -e 'ri '\$target_chain_subdomain1_residues1'|ri '\$target_chain_subdomain1_residues2'|ri '\$target_chain_subdomain1_residues3'\n name '\$((start_index + 2)) 'Chain_'\$target_chain'_SD1\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -n "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -nobackup
    echo -e 'ri '\$target_chain_subdomain2_residues'\n name '\$((start_index + 3)) 'Chain_'\$target_chain'_SD2\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -n "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -nobackup
    echo -e 'ri '\$target_chain_subdomain3_residues1'|ri '\$target_chain_subdomain3_residues2'\n name '\$((start_index + 4)) 'Chain_'\$target_chain'_SD3\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -n "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -nobackup
    echo -e 'ri '\$target_chain_subdomain4_residues'\n name '\$((start_index + 5)) 'Chain_'\$target_chain'_SD4\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -n "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -nobackup
    echo -e '"Backbone" & "Chain_'\$target_chain'_protein"\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -n "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -nobackup
    echo -e '"Backbone" & "Chain_'\$target_chain'_SD1"\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -n "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -nobackup
    echo -e '"Backbone" & "Chain_'\$target_chain'_SD2"\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -n "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -nobackup
    echo -e '"Backbone" & "Chain_'\$target_chain'_SD3"\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -n "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -nobackup
    echo -e '"Backbone" & "Chain_'\$target_chain'_SD4"\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -n "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -nobackup
    echo -e '"C-alpha" & "Chain_'\$target_chain'_protein"\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -n "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -nobackup
    echo -e '"C-alpha" & "Chain_'\$target_chain'_SD1"\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -n "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -nobackup
    echo -e '"C-alpha" & "Chain_'\$target_chain'_SD2"\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -n "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -nobackup
    echo -e '"C-alpha" & "Chain_'\$target_chain'_SD3"\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -n "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -nobackup
    echo -e '"C-alpha" & "Chain_'\$target_chain'_SD4"\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -n "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -nobackup


  fi


  #export protein only

  echo -e 'Protein\n Protein' | gmx trjconv -s "\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o "\$mutant"_"\$Factin_or_Gactin"_md_trajec_protein.xtc -n "\$mutant"_"\$Factin_or_Gactin"_index.ndx -fit rot+trans
  echo -e 'Protein\n Protein' | gmx trjconv -s "\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o "\$mutant"_"\$Factin_or_Gactin"_md_trajec_protein.pdb -n "\$mutant"_"\$Factin_or_Gactin"_index.ndx -fit rot+trans -dump 0

  #fit the trajectory for the target chain only

  echo -e 'Chain_'\$target_chain'\n Chain_'\$target_chain | gmx trjconv -s "\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o "\$mutant"_Chain"\$target_chain"_trajec_fit.xtc -n "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -fit rot+trans
  echo -e 'Chain_'\$target_chain'\n Chain_'\$target_chain | gmx trjconv -s "\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o "\$mutant"_Chain"\$target_chain"_FIRST_FRAME.pdb -n "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -fit rot+trans -dump 0
  echo -e 'Chain_'\$target_chain'\n Chain_'\$target_chain | gmx trjconv -s "\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o "\$mutant"_Chain"\$target_chain"_FINAL_FRAME.pdb -n "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -fit rot+trans -tu ns -dump "\$simulation_length"

  #export protein only for target chain

  echo -e 'Chain_'\$target_chain'_protein\n Chain_'\$target_chain'_protein' | gmx trjconv -s "\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o "\$mutant"_Chain"\$target_chain"_trajec_protein.xtc -n "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -fit rot+trans
  echo -e 'Chain_'\$target_chain'_protein\n Chain_'\$target_chain'_protein' | gmx trjconv -s "\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o "\$mutant"_Chain"\$target_chain"_protein_FIRST_FRAME.pdb -n "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -fit rot+trans -dump 0
  echo -e 'Chain_'\$target_chain'_protein\n Chain_'\$target_chain'_protein' | gmx trjconv -s "\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o "\$mutant"_Chain"\$target_chain"_protein_FINAL_FRAME.pdb -n "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -fit rot+trans -tu ns -dump "\$simulation_length"



############################################ANALYSES#############################################################

  #RMSD calculations
  mkdir "\$directory_path"/"\$mutant"/"\$replicate"/RMSD
  echo -e 'Backbone\n Backbone_&_Chain_'\$target_chain'_protein' | gmx rms -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -s "\$mutant"_"\$Factin_or_Gactin"_md.tpr -n "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -o "\$mutant"_Chain_"\$target_chain"_protein_rmsd_backbone_"\$rep".xvg
  echo -e 'Backbone\n Backbone_&_Chain_'\$target_chain'_SD1' | gmx rms -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -s "\$mutant"_"\$Factin_or_Gactin"_md.tpr -n "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -o "\$mutant"_Chain_"\$target_chain"_SD1_rmsd_backbone_"\$rep".xvg
  echo -e 'Backbone\n Backbone_&_Chain_'\$target_chain'_SD2' | gmx rms -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -s "\$mutant"_"\$Factin_or_Gactin"_md.tpr -n "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -o "\$mutant"_Chain_"\$target_chain"_SD2_rmsd_backbone_"\$rep".xvg
  echo -e 'Backbone\n Backbone_&_Chain_'\$target_chain'_SD3' | gmx rms -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -s "\$mutant"_"\$Factin_or_Gactin"_md.tpr -n "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -o "\$mutant"_Chain_"\$target_chain"_SD3_rmsd_backbone_"\$rep".xvg
  echo -e 'Backbone\n Backbone_&_Chain_'\$target_chain'_SD4' | gmx rms -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -s "\$mutant"_"\$Factin_or_Gactin"_md.tpr -n "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -o "\$mutant"_Chain_"\$target_chain"_SD4_rmsd_backbone_"\$rep".xvg
  mv "\$directory_path"/"\$mutant"/"\$replicate"/*_rmsd_backbone_"\$rep".xvg "\$directory_path"/"\$mutant"/"\$replicate"/RMSD

###RMSF calculations for the whole filament. Calculate for whole simulation as well as the last 30% (representing equilibrium)
  mkdir "\$directory_path"/"\$mutant"/"\$replicate"/RMSF
  echo -e 'Protein' | gmx rmsf -s "\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -n "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -o "\$mutant"_Chain_"\$target_chain"_rmsf_"\$rep"_WHOLE_FILAMENT_Protein.xvg -oq "\$mutant"_Chain_"\$target_chain"_BFACTORS_"\$rep"_WHOLE_FILAMENT_Protein.pdb -od "\$mutant"_Chain_"\$target_chain"_rmsdev_"\$rep"_WHOLE_FILAMENT_Protein.xvg -res
  echo -e 'Protein' | gmx rmsf -s "\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -n "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -o "\$mutant"_Chain_"\$target_chain"_rmsf_"\$rep"_WHOLE_FILAMENT_Protein_EQ.xvg -oq "\$mutant"_Chain_"\$target_chain"_BFACTORS_"\$rep"_WHOLE_FILAMENT_Protein_EQ.pdb -od "\$mutant"_Chain_"\$target_chain"_rmsdev_"\$rep"_WHOLE_FILAMENT_Protein_EQ.xvg -res -b \$((eq_start_time*1000)) -e \$((simulation_length*1000))
  echo -e 'C-alpha' | gmx rmsf -s "\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -n "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -o "\$mutant"_Chain_"\$target_chain"_rmsf_"\$rep"_WHOLE_FILAMENT_C-alpha.xvg -oq "\$mutant"_Chain_"\$target_chain"_BFACTORS_"\$rep"_WHOLE_FILAMENT_C-alpha.pdb -od "\$mutant"_Chain_"\$target_chain"_rmsdev_"\$rep"_WHOLE_FILAMENT_C-alpha.xvg -res
  echo -e 'C-alpha' | gmx rmsf -s "\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -n "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -o "\$mutant"_Chain_"\$target_chain"_rmsf_"\$rep"_WHOLE_FILAMENT_C-alpha_EQ.xvg -oq "\$mutant"_Chain_"\$target_chain"_BFACTORS_"\$rep"_WHOLE_FILAMENT_C-alpha_EQ.pdb -od "\$mutant"_Chain_"\$target_chain"_rmsdev_"\$rep"_WHOLE_FILAMENT_C-alpha_EQ.xvg -res -b \$((eq_start_time*1000)) -e \$((simulation_length*1000))
  mv "\$directory_path"/"\$mutant"/"\$replicate"/*"\$rep"_WHOLE_FILAMENT_* "\$directory_path"/"\$mutant"/"\$replicate"/RMSF

  #hbonds
  mkdir "\$directory_path"/"\$mutant"/"\$replicate"/HBONDS
  #Hbonds within target chain
  echo -e 'Chain_'\$target_chain'_protein\n Chain_'\$target_chain'_protein' | gmx hbond -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -s "\$mutant"_"\$Factin_or_Gactin"_md.tpr -n "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -num "\$mutant"_"\$Factin_or_Gactin"_hbonds_Chain_"\$target_chain"_Chain_"\$target_chain"_"\$rep".xvg -dist "\$mutant"_"\$Factin_or_Gactin"_hbonds_distribution_Chain_"\$target_chain"_Chain_"\$target_chain"_"\$rep".xvg -ang "\$mutant"_"\$Factin_or_Gactin"_hbonds_angle_distribution_Chain_"\$target_chain"_Chain_"\$target_chain"_"\$rep".xvg -hx "\$mutant"_"\$Factin_or_Gactin"_hbonds_helices_Chain_"\$target_chain"_Chain_"\$target_chain"_"\$rep".xvg -dan "\$mutant"_"\$Factin_or_Gactin"_hbonds_donor_acceptor_Chain_"\$target_chain"_Chain_"\$target_chain"_"\$rep".xvg
  #Hbonds within target chain SD1
  echo -e 'Chain_'\$target_chain'_SD1\n Chain_'\$target_chain'_SD1' | gmx hbond -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -s "\$mutant"_"\$Factin_or_Gactin"_md.tpr -n "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -num "\$mutant"_"\$Factin_or_Gactin"_hbonds_Chain_"\$target_chain"_SD1_Chain_"\$target_chain"_SD1_"\$rep".xvg -dist "\$mutant"_"\$Factin_or_Gactin"_hbonds_distribution_Chain_"\$target_chain"_SD1_Chain_"\$target_chain"_SD1_"\$rep".xvg -ang "\$mutant"_"\$Factin_or_Gactin"_hbonds_angle_distribution_Chain_"\$target_chain"_SD1_Chain_"\$target_chain"_SD1_"\$rep".xvg -hx "\$mutant"_"\$Factin_or_Gactin"_hbonds_helices_Chain_"\$target_chain"_SD1_Chain_"\$target_chain"_SD1_"\$rep".xvg -dan "\$mutant"_"\$Factin_or_Gactin"_hbonds_donor_acceptor_Chain_"\$target_chain"_SD1_Chain_"\$target_chain"_SD1_"\$rep".xvg
  #Hbonds within target chain SD2
  echo -e 'Chain_'\$target_chain'_SD2\n Chain_'\$target_chain'_SD2' | gmx hbond -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -s "\$mutant"_"\$Factin_or_Gactin"_md.tpr -n "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -num "\$mutant"_"\$Factin_or_Gactin"_hbonds_Chain_"\$target_chain"_SD2_Chain_"\$target_chain"_SD2_"\$rep".xvg -dist "\$mutant"_"\$Factin_or_Gactin"_hbonds_distribution_Chain_"\$target_chain"_SD2_Chain_"\$target_chain"_SD2_"\$rep".xvg -ang "\$mutant"_"\$Factin_or_Gactin"_hbonds_angle_distribution_Chain_"\$target_chain"_SD2_Chain_"\$target_chain"_SD2_"\$rep".xvg -hx "\$mutant"_"\$Factin_or_Gactin"_hbonds_helices_Chain_"\$target_chain"_SD2_Chain_"\$target_chain"_SD2_"\$rep".xvg -dan "\$mutant"_"\$Factin_or_Gactin"_hbonds_donor_acceptor_Chain_"\$target_chain"_SD2_Chain_"\$target_chain"_SD2_"\$rep".xvg
  #Hbonds within target chain SD3
  echo -e 'Chain_'\$target_chain'_SD3\n Chain_'\$target_chain'_SD3' | gmx hbond -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -s "\$mutant"_"\$Factin_or_Gactin"_md.tpr -n "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -num "\$mutant"_"\$Factin_or_Gactin"_hbonds_Chain_"\$target_chain"_SD3_Chain_"\$target_chain"_SD3_"\$rep".xvg -dist "\$mutant"_"\$Factin_or_Gactin"_hbonds_distribution_Chain_"\$target_chain"_SD3_Chain_"\$target_chain"_SD3_"\$rep".xvg -ang "\$mutant"_"\$Factin_or_Gactin"_hbonds_angle_distribution_Chain_"\$target_chain"_SD3_Chain_"\$target_chain"_SD3_"\$rep".xvg -hx "\$mutant"_"\$Factin_or_Gactin"_hbonds_helices_Chain_"\$target_chain"_SD3_Chain_"\$target_chain"_SD3_"\$rep".xvg -dan "\$mutant"_"\$Factin_or_Gactin"_hbonds_donor_acceptor_Chain_"\$target_chain"_SD3_Chain_"\$target_chain"_SD3_"\$rep".xvg
  #Hbonds within target chain SD4
  echo -e 'Chain_'\$target_chain'_SD4\n Chain_'\$target_chain'_SD4' | gmx hbond -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -s "\$mutant"_"\$Factin_or_Gactin"_md.tpr -n "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -num "\$mutant"_"\$Factin_or_Gactin"_hbonds_Chain_"\$target_chain"_SD4_Chain_"\$target_chain"_SD4_"\$rep".xvg -dist "\$mutant"_"\$Factin_or_Gactin"_hbonds_distribution_Chain_"\$target_chain"_SD4_Chain_"\$target_chain"_SD4_"\$rep".xvg -ang "\$mutant"_"\$Factin_or_Gactin"_hbonds_angle_distribution_Chain_"\$target_chain"_SD4_Chain_"\$target_chain"_SD4_"\$rep".xvg -hx "\$mutant"_"\$Factin_or_Gactin"_hbonds_helices_Chain_"\$target_chain"_SD4_Chain_"\$target_chain"_SD4_"\$rep".xvg -dan "\$mutant"_"\$Factin_or_Gactin"_hbonds_donor_acceptor_Chain_"\$target_chain"_SD4_Chain_"\$target_chain"_SD4_"\$rep".xvg
  #Hbonds between nucleotide and target chain SD1
  echo -e '13\n Chain_'\$target_chain'_SD1' | gmx hbond -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -s "\$mutant"_"\$Factin_or_Gactin"_md.tpr -n "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -num "\$mutant"_"\$Factin_or_Gactin"_hbonds_"\$nucleotide"_Chain_"\$target_chain"_SD1_"\$rep".xvg -dist "\$mutant"_"\$Factin_or_Gactin"_hbonds_distribution_"\$nucleotide"_Chain_"\$target_chain"_SD1_"\$rep".xvg -ang "\$mutant"_"\$Factin_or_Gactin"_hbonds_angle_distribution_"\$nucleotide"_Chain_"\$target_chain"_SD1_"\$rep".xvg -hx "\$mutant"_"\$Factin_or_Gactin"_hbonds_helices_"\$nucleotide"_Chain_"\$target_chain"_SD1_"\$rep".xvg -dan "\$mutant"_"\$Factin_or_Gactin"_hbonds_donor_acceptor_"\$nucleotide"_Chain_"\$target_chain"_SD1_"\$rep".xvg
  #Hbonds between nucleotide and target chain SD2
  echo -e '13\n Chain_'\$target_chain'_SD2' | gmx hbond -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -s "\$mutant"_"\$Factin_or_Gactin"_md.tpr -n "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -num "\$mutant"_"\$Factin_or_Gactin"_hbonds_"\$nucleotide"_Chain_"\$target_chain"_SD2_"\$rep".xvg -dist "\$mutant"_"\$Factin_or_Gactin"_hbonds_distribution_"\$nucleotide"_Chain_"\$target_chain"_SD2_"\$rep".xvg -ang "\$mutant"_"\$Factin_or_Gactin"_hbonds_angle_distribution_"\$nucleotide"_Chain_"\$target_chain"_SD2_"\$rep".xvg -hx "\$mutant"_"\$Factin_or_Gactin"_hbonds_helices_"\$nucleotide"_Chain_"\$target_chain"_SD2_"\$rep".xvg -dan "\$mutant"_"\$Factin_or_Gactin"_hbonds_donor_acceptor_"\$nucleotide"_Chain_"\$target_chain"_SD2_"\$rep".xvg
  #Hbonds between nucleotide and target chain SD3
  echo -e '13\n Chain_'\$target_chain'_SD3' | gmx hbond -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -s "\$mutant"_"\$Factin_or_Gactin"_md.tpr -n "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -num "\$mutant"_"\$Factin_or_Gactin"_hbonds_"\$nucleotide"_Chain_"\$target_chain"_SD3_"\$rep".xvg -dist "\$mutant"_"\$Factin_or_Gactin"_hbonds_distribution_"\$nucleotide"_Chain_"\$target_chain"_SD3_"\$rep".xvg -ang "\$mutant"_"\$Factin_or_Gactin"_hbonds_angle_distribution_"\$nucleotide"_Chain_"\$target_chain"_SD3_"\$rep".xvg -hx "\$mutant"_"\$Factin_or_Gactin"_hbonds_helices_"\$nucleotide"_Chain_"\$target_chain"_SD3_"\$rep".xvg -dan "\$mutant"_"\$Factin_or_Gactin"_hbonds_donor_acceptor_"\$nucleotide"_Chain_"\$target_chain"_SD3_"\$rep".xvg
  #Hbonds between nucleotide and target chain SD4
  echo -e '13\n Chain_'\$target_chain'_SD4' | gmx hbond -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -s "\$mutant"_"\$Factin_or_Gactin"_md.tpr -n "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -num "\$mutant"_"\$Factin_or_Gactin"_hbonds_"\$nucleotide"_Chain_"\$target_chain"_SD4_"\$rep".xvg -dist "\$mutant"_"\$Factin_or_Gactin"_hbonds_distribution_"\$nucleotide"_Chain_"\$target_chain"_SD4_"\$rep".xvg -ang "\$mutant"_"\$Factin_or_Gactin"_hbonds_angle_distribution_"\$nucleotide"_Chain_"\$target_chain"_SD4_"\$rep".xvg -hx "\$mutant"_"\$Factin_or_Gactin"_hbonds_helices_"\$nucleotide"_Chain_"\$target_chain"_SD4_"\$rep".xvg -dan "\$mutant"_"\$Factin_or_Gactin"_hbonds_donor_acceptor_"\$nucleotide"_Chain_"\$target_chain"_SD4_"\$rep".xvg

    #Hbonds between PO4 and target chain SD1
  if [ "\$nucleotide_state" == 'ADP-Pi' ]; then
      #Hbonds between PO4 (index number 15) and target chain SD1
    echo -e '15\n Chain_'\$target_chain'_SD1' | gmx hbond -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -s "\$mutant"_"\$Factin_or_Gactin"_md.tpr -n "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -num "\$mutant"_"\$Factin_or_Gactin"_hbonds_PO4_Chain_"\$target_chain"_SD1_"\$rep".xvg -dist "\$mutant"_"\$Factin_or_Gactin"_hbonds_distribution_PO4_Chain_"\$target_chain"_SD1_"\$rep".xvg -ang "\$mutant"_"\$Factin_or_Gactin"_hbonds_angle_distribution_PO4_Chain_"\$target_chain"_SD1_"\$rep".xvg -hx "\$mutant"_"\$Factin_or_Gactin"_hbonds_helices_PO4_Chain_"\$target_chain"_SD1_"\$rep".xvg -dan "\$mutant"_"\$Factin_or_Gactin"_hbonds_donor_acceptor_PO4_Chain_"\$target_chain"_SD1_"\$rep".xvg
    #Hbonds between nucleotide and target chain SD2
    echo -e '15\n Chain_'\$target_chain'_SD2' | gmx hbond -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -s "\$mutant"_"\$Factin_or_Gactin"_md.tpr -n "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -num "\$mutant"_"\$Factin_or_Gactin"_hbonds_PO4_Chain_"\$target_chain"_SD2_"\$rep".xvg -dist "\$mutant"_"\$Factin_or_Gactin"_hbonds_distribution_PO4_Chain_"\$target_chain"_SD2_"\$rep".xvg -ang "\$mutant"_"\$Factin_or_Gactin"_hbonds_angle_distribution_PO4_Chain_"\$target_chain"_SD2_"\$rep".xvg -hx "\$mutant"_"\$Factin_or_Gactin"_hbonds_helices_PO4_Chain_"\$target_chain"_SD2_"\$rep".xvg -dan "\$mutant"_"\$Factin_or_Gactin"_hbonds_donor_acceptor_PO4_Chain_"\$target_chain"_SD2_"\$rep".xvg
    #Hbonds between nucleotide and target chain SD3
    echo -e '15\n Chain_'\$target_chain'_SD3' | gmx hbond -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -s "\$mutant"_"\$Factin_or_Gactin"_md.tpr -n "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -num "\$mutant"_"\$Factin_or_Gactin"_hbonds_PO4_Chain_"\$target_chain"_SD3_"\$rep".xvg -dist "\$mutant"_"\$Factin_or_Gactin"_hbonds_distribution_PO4_Chain_"\$target_chain"_SD3_"\$rep".xvg -ang "\$mutant"_"\$Factin_or_Gactin"_hbonds_angle_distribution_PO4_Chain_"\$target_chain"_SD3_"\$rep".xvg -hx "\$mutant"_"\$Factin_or_Gactin"_hbonds_helices_PO4_Chain_"\$target_chain"_SD3_"\$rep".xvg -dan "\$mutant"_"\$Factin_or_Gactin"_hbonds_donor_acceptor_PO4_Chain_"\$target_chain"_SD3_"\$rep".xvg
    #Hbonds between nucleotide and target chain SD4
    echo -e '15\n Chain_'\$target_chain'_SD4' | gmx hbond -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -s "\$mutant"_"\$Factin_or_Gactin"_md.tpr -n "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -num "\$mutant"_"\$Factin_or_Gactin"_hbonds_PO4_Chain_"\$target_chain"_SD4_"\$rep".xvg -dist "\$mutant"_"\$Factin_or_Gactin"_hbonds_distribution_PO4_Chain_"\$target_chain"_SD4_"\$rep".xvg -ang "\$mutant"_"\$Factin_or_Gactin"_hbonds_angle_distribution_PO4_Chain_"\$target_chain"_SD4_"\$rep".xvg -hx "\$mutant"_"\$Factin_or_Gactin"_hbonds_helices_PO4_Chain_"\$target_chain"_SD4_"\$rep".xvg -dan "\$mutant"_"\$Factin_or_Gactin"_hbonds_donor_acceptor_PO4_Chain_"\$target_chain"_SD4_"\$rep".xvg
  fi
  mv "\$directory_path"/"\$mutant"/"\$replicate"/*"\$mutant"_"\$Factin_or_Gactin"_hbonds_* "\$directory_path"/"\$mutant"/"\$replicate"/HBONDS




  #radius of gyration
  mkdir "\$directory_path"/"\$mutant"/"\$replicate"/RADIUS
  echo -e 'Protein' | gmx gyrate -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -s "\$mutant"_"\$Factin_or_Gactin"_md.tpr -p -o "\$mutant"_"\$Factin_or_Gactin"_radius.xvg -n "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx
  echo -e 'Chain_'\$target_chain'_protein' | gmx gyrate -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -s "\$mutant"_"\$Factin_or_Gactin"_md.tpr -p -o "\$mutant"_Chain_"\$target_chain"_radius.xvg -n "\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx
  mv "\$directory_path"/"\$mutant"/"\$replicate"/*radius.xvg* "\$directory_path"/"\$mutant"/"\$replicate"/RADIUS




  ###############distance analysis

  mkdir "\$directory_path"/"\$mutant"/"\$replicate"/DISTANCES
  mkdir "\$directory_path"/"\$mutant"/"\$replicate"/DISTANCES/INTERPROTOMER
  mkdir "\$directory_path"/"\$mutant"/"\$replicate"/DISTANCES/NUCLEOTIDE
  mkdir "\$directory_path"/"\$mutant"/"\$replicate"/DISTANCES/INTERSTRAND
  mkdir "\$directory_path"/"\$mutant"/"\$replicate"/DISTANCES/ALLOSTERY
  #calculate distances between interacting residues at interprotomer interface
  #set if loop so that target D-loop only if target chain number is greater than 2
  #if target chain is 1 or 2, either G-actin or pointed end
  #need to calculate the residue numbers
  #For interprotomer interface, want target chain's C-terminus and D-loop, then C-terminus/D-loop of neighbouring protomers
  if [[ "\$target_chain_number" -gt 2 && "\$number_of_chains" -gt 3 ]]; then
    target_F375=\$((number_of_residues * target_chain_number))
    target_T351=\$((351 + (number_of_residues * (target_chain_number - 1))))
    target_F352=\$((352 + (number_of_residues * (target_chain_number - 1))))
    target_L171=\$((171 + (number_of_residues * (target_chain_number - 1))))
    target_A170=\$((170 + (number_of_residues * (target_chain_number - 1))))
    target_Y169=\$((169 + (number_of_residues * (target_chain_number - 1))))
    target_E167=\$((167 + (number_of_residues * (target_chain_number - 1))))
    interacting_V43=\$((43 + (number_of_residues * (target_chain_number + 1))))
    interacting_Q41=\$((41 + (number_of_residues * (target_chain_number + 1))))
    interacting_V45=\$((45 + (number_of_residues * (target_chain_number + 1))))
    interacting_H40=\$((40 + (number_of_residues * (target_chain_number + 1))))
    interacting_Q49=\$((49 + (number_of_residues * (target_chain_number + 1))))
    interacting_L50=\$((50 + (number_of_residues * (target_chain_number + 1))))
    interacting_P38=\$((38 + (number_of_residues * (target_chain_number + 1))))

    interacting_F375=\$((375 + (number_of_residues * (target_chain_number - 3))))
    interacting_T351=\$((351 + (number_of_residues * (target_chain_number - 3))))
    interacting_F352=\$((352 + (number_of_residues * (target_chain_number - 3))))
    interacting_L171=\$((171 + (number_of_residues * (target_chain_number - 3))))
    interacting_A170=\$((170 + (number_of_residues * (target_chain_number - 3))))
    interacting_Y169=\$((169 + (number_of_residues * (target_chain_number - 3))))
    interacting_E167=\$((167 + (number_of_residues * (target_chain_number - 3))))
    target_V43=\$((43 + (number_of_residues * (target_chain_number - 1))))
    target_Q41=\$((41 + (number_of_residues * (target_chain_number - 1))))
    target_V45=\$((45 + (number_of_residues * (target_chain_number - 1))))
    target_H40=\$((40 + (number_of_residues * (target_chain_number - 1))))
    target_Q49=\$((49 + (number_of_residues * (target_chain_number - 1))))
    target_L50=\$((50 + (number_of_residues * (target_chain_number - 1))))
    target_P38=\$((38 + (number_of_residues * (target_chain_number - 1))))

    echo -e 'ri '\$target_F375'\n name '\$start_index 'Chain_'\$target_chain'_F375\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$target_T351'\n name '\$((start_index + 1)) 'Chain_'\$target_chain'_T351\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$target_F352'\n name '\$((start_index + 2)) 'Chain_'\$target_chain'_F352\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$target_L171'\n name '\$((start_index + 3)) 'Chain_'\$target_chain'_L171\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$target_A170'\n name '\$((start_index + 4)) 'Chain_'\$target_chain'_A170\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$target_Y169'\n name '\$((start_index + 5)) 'Chain_'\$target_chain'_Y169\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$target_E167'\n name '\$((start_index + 6)) 'Chain_'\$target_chain'_E167\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$interacting_F375'\n name '\$((start_index + 7)) 'Chain_'\$interacting_chain_minus2'_F375\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$interacting_T351'\n name '\$((start_index + 8)) 'Chain_'\$interacting_chain_minus2'_T351\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$interacting_F352'\n name '\$((start_index + 9)) 'Chain_'\$interacting_chain_minus2'_F352\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$interacting_L171'\n name '\$((start_index + 10)) 'Chain_'\$interacting_chain_minus2'_L171\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$interacting_A170'\n name '\$((start_index + 11)) 'Chain_'\$interacting_chain_minus2'_A170\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$interacting_Y169'\n name '\$((start_index + 12)) 'Chain_'\$interacting_chain_minus2'_Y169\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$interacting_E167'\n name '\$((start_index + 13)) 'Chain_'\$interacting_chain_minus2'_E167\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$target_V43'\n name '\$((start_index + 14)) 'Chain_'\$target_chain'_V43\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$target_Q41'\n name '\$((start_index + 15)) 'Chain_'\$target_chain'_Q41\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$target_V45'\n name '\$((start_index + 16)) 'Chain_'\$target_chain'_V45\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$target_H40'\n name '\$((start_index + 17)) 'Chain_'\$target_chain'_H40\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$target_Q49'\n name '\$((start_index + 18)) 'Chain_'\$target_chain'_Q49\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$target_L50'\n name '\$((start_index + 19)) 'Chain_'\$target_chain'_L50\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$target_P38'\n name '\$((start_index + 20)) 'Chain_'\$target_chain'_P38\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$interacting_V43'\n name '\$((start_index + 21)) 'Chain_'\$interacting_chain_plus2'_V43\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$interacting_Q41'\n name '\$((start_index + 22)) 'Chain_'\$interacting_chain_plus2'_Q41\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$interacting_V45'\n name '\$((start_index + 23)) 'Chain_'\$interacting_chain_plus2'_V45\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$interacting_H40'\n name '\$((start_index + 24)) 'Chain_'\$interacting_chain_plus2'_H40\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$interacting_Q49'\n name '\$((start_index + 25)) 'Chain_'\$interacting_chain_plus2'_Q49\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$interacting_L50'\n name '\$((start_index + 26)) 'Chain_'\$interacting_chain_plus2'_L50\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$interacting_P38'\n name '\$((start_index + 27)) 'Chain_'\$interacting_chain_plus2'_P38\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -nobackup

    gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/INTERPROTOMER/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_F375_Chain_"\$interacting_chain_plus2"_V43_"\$rep".xvg -ref Chain_"\$target_chain"_F375 -sel Chain_"\$interacting_chain_plus2"_V43
    gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/INTERPROTOMER/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_F375_Chain_"\$interacting_chain_plus2"_Q41_"\$rep".xvg -ref Chain_"\$target_chain"_F375 -sel Chain_"\$interacting_chain_plus2"_Q41
    gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/INTERPROTOMER/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_T351_Chain_"\$interacting_chain_plus2"_V45_"\$rep".xvg -ref Chain_"\$target_chain"_T351 -sel Chain_"\$interacting_chain_plus2"_V45
    gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/INTERPROTOMER/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_F352_Chain_"\$interacting_chain_plus2"_V45_"\$rep".xvg -ref Chain_"\$target_chain"_F352 -sel Chain_"\$interacting_chain_plus2"_V45
    gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/INTERPROTOMER/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_L171_Chain_"\$interacting_chain_plus2"_H40_"\$rep".xvg -ref Chain_"\$target_chain"_L171 -sel Chain_"\$interacting_chain_plus2"_H40
    gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/INTERPROTOMER/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_A170_Chain_"\$interacting_chain_plus2"_Q41_"\$rep".xvg -ref Chain_"\$target_chain"_A170 -sel Chain_"\$interacting_chain_plus2"_Q41
    gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/INTERPROTOMER/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_Y169_Chain_"\$interacting_chain_plus2"_Q49_"\$rep".xvg -ref Chain_"\$target_chain"_Y169 -sel Chain_"\$interacting_chain_plus2"_Q49
    gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/INTERPROTOMER/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_Y169_Chain_"\$interacting_chain_plus2"_L50_"\$rep".xvg -ref Chain_"\$target_chain"_Y169 -sel Chain_"\$interacting_chain_plus2"_L50
    gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/INTERPROTOMER/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_Y169_Chain_"\$interacting_chain_plus2"_P38_"\$rep".xvg -ref Chain_"\$target_chain"_Y169 -sel Chain_"\$interacting_chain_plus2"_P38
    gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/INTERPROTOMER/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_E167_Chain_"\$interacting_chain_plus2"_Q49_"\$rep".xvg -ref Chain_"\$target_chain"_E167 -sel Chain_"\$interacting_chain_plus2"_Q49

    gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/INTERPROTOMER/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$interacting_chain_minus2"_F375_Chain_"\$target_chain"_V43_"\$rep".xvg -ref Chain_"\$interacting_chain_minus2"_F375 -sel Chain_"\$target_chain"_V43
    gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/INTERPROTOMER/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$interacting_chain_minus2"_F375_Chain_"\$target_chain"_Q41_"\$rep".xvg -ref Chain_"\$interacting_chain_minus2"_F375 -sel Chain_"\$target_chain"_Q41
    gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/INTERPROTOMER/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$interacting_chain_minus2"_T351_Chain_"\$target_chain"_V45_"\$rep".xvg -ref Chain_"\$interacting_chain_minus2"_T351 -sel Chain_"\$target_chain"_V45
    gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/INTERPROTOMER/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$interacting_chain_minus2"_F352_Chain_"\$target_chain"_V45_"\$rep".xvg -ref Chain_"\$interacting_chain_minus2"_F352 -sel Chain_"\$target_chain"_V45
    gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/INTERPROTOMER/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$interacting_chain_minus2"_L171_Chain_"\$target_chain"_H40_"\$rep".xvg -ref Chain_"\$interacting_chain_minus2"_L171 -sel Chain_"\$target_chain"_H40
    gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/INTERPROTOMER/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$interacting_chain_minus2"_A170_Chain_"\$target_chain"_Q41_"\$rep".xvg -ref Chain_"\$interacting_chain_minus2"_A170 -sel Chain_"\$target_chain"_Q41
    gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/INTERPROTOMER/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$interacting_chain_minus2"_Y169_Chain_"\$target_chain"_Q49_"\$rep".xvg -ref Chain_"\$interacting_chain_minus2"_Y169 -sel Chain_"\$target_chain"_Q49
    gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/INTERPROTOMER/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$interacting_chain_minus2"_Y169_Chain_"\$target_chain"_L50_"\$rep".xvg -ref Chain_"\$interacting_chain_minus2"_Y169 -sel Chain_"\$target_chain"_L50
    gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/INTERPROTOMER/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$interacting_chain_minus2"_Y169_Chain_"\$target_chain"_P38_"\$rep".xvg -ref Chain_"\$interacting_chain_minus2"_Y169 -sel Chain_"\$target_chain"_P38
    gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/INTERPROTOMER/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$interacting_chain_minus2"_E167_Chain_"\$target_chain"_Q49_"\$rep".xvg -ref Chain_"\$interacting_chain_minus2"_E167 -sel Chain_"\$target_chain"_Q49

#otherwise (if the first loop doesnt go), for any chain number AND if number of chains is greater than 3 (meaning there will be at least 4 chains and therefore at least 2 interprotomer interfaces)
  elif [[ "\$target_chain_number" -gt 0 && "\$number_of_chains" -gt 3 ]]; then
    #this loop goes target chain is pointed end (chain A/B). Only looks at target chain's C-terminus and the D-loop below it
    target_F375=\$((number_of_residues * target_chain_number))
    target_T351=\$((351 + (number_of_residues * (target_chain_number - 1))))
    target_F352=\$((352 + (number_of_residues * (target_chain_number - 1))))
    target_L171=\$((171 + (number_of_residues * (target_chain_number - 1))))
    target_A170=\$((170 + (number_of_residues * (target_chain_number - 1))))
    target_Y169=\$((169 + (number_of_residues * (target_chain_number - 1))))
    target_E167=\$((167 + (number_of_residues * (target_chain_number - 1))))
    interacting_V43=\$((43 + (number_of_residues * (target_chain_number + 1))))
    interacting_Q41=\$((41 + (number_of_residues * (target_chain_number + 1))))
    interacting_V45=\$((45 + (number_of_residues * (target_chain_number + 1))))
    interacting_H40=\$((40 + (number_of_residues * (target_chain_number + 1))))
    interacting_Q49=\$((49 + (number_of_residues * (target_chain_number + 1))))
    interacting_L50=\$((50 + (number_of_residues * (target_chain_number + 1))))
    interacting_P38=\$((38 + (number_of_residues * (target_chain_number + 1))))

    echo -e 'ri '\$target_F375'\n name '\$start_index 'Chain_'\$target_chain'_F375\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$target_T351'\n name '\$((start_index + 1)) 'Chain_'\$target_chain'_T351\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$target_F352'\n name '\$((start_index + 2)) 'Chain_'\$target_chain'_F352\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$target_L171'\n name '\$((start_index + 3)) 'Chain_'\$target_chain'_L171\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$target_A170'\n name '\$((start_index + 4)) 'Chain_'\$target_chain'_A170\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$target_Y169'\n name '\$((start_index + 5)) 'Chain_'\$target_chain'_Y169\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$target_E167'\n name '\$((start_index + 6)) 'Chain_'\$target_chain'_E167\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$interacting_V43'\n name '\$((start_index + 21)) 'Chain_'\$interacting_chain_plus2'_V43\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$interacting_Q41'\n name '\$((start_index + 22)) 'Chain_'\$interacting_chain_plus2'_Q41\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$interacting_V45'\n name '\$((start_index + 23)) 'Chain_'\$interacting_chain_plus2'_V45\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$interacting_H40'\n name '\$((start_index + 24)) 'Chain_'\$interacting_chain_plus2'_H40\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$interacting_Q49'\n name '\$((start_index + 25)) 'Chain_'\$interacting_chain_plus2'_Q49\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$interacting_L50'\n name '\$((start_index + 26)) 'Chain_'\$interacting_chain_plus2'_L50\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$interacting_P38'\n name '\$((start_index + 27)) 'Chain_'\$interacting_chain_plus2'_P38\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -nobackup

    gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/INTERPROTOMER/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_F375_Chain_"\$interacting_chain_plus2"_V43_"\$rep".xvg -ref Chain_"\$target_chain"_F375 -sel Chain_"\$interacting_chain_plus2"_V43
    gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/INTERPROTOMER/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_F375_Chain_"\$interacting_chain_plus2"_Q41_"\$rep".xvg -ref Chain_"\$target_chain"_F375 -sel Chain_"\$interacting_chain_plus2"_Q41
    gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/INTERPROTOMER/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_T351_Chain_"\$interacting_chain_plus2"_V45_"\$rep".xvg -ref Chain_"\$target_chain"_T351 -sel Chain_"\$interacting_chain_plus2"_V45
    gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/INTERPROTOMER/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_F352_Chain_"\$interacting_chain_plus2"_V45_"\$rep".xvg -ref Chain_"\$target_chain"_F352 -sel Chain_"\$interacting_chain_plus2"_V45
    gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/INTERPROTOMER/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_L171_Chain_"\$interacting_chain_plus2"_H40_"\$rep".xvg -ref Chain_"\$target_chain"_L171 -sel Chain_"\$interacting_chain_plus2"_H40
    gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/INTERPROTOMER/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_A170_Chain_"\$interacting_chain_plus2"_Q41_"\$rep".xvg -ref Chain_"\$target_chain"_A170 -sel Chain_"\$interacting_chain_plus2"_Q41
    gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/INTERPROTOMER/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_Y169_Chain_"\$interacting_chain_plus2"_Q49_"\$rep".xvg -ref Chain_"\$target_chain"_Y169 -sel Chain_"\$interacting_chain_plus2"_Q49
    gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/INTERPROTOMER/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_Y169_Chain_"\$interacting_chain_plus2"_L50_"\$rep".xvg -ref Chain_"\$target_chain"_Y169 -sel Chain_"\$interacting_chain_plus2"_L50
    gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/INTERPROTOMER/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_Y169_Chain_"\$interacting_chain_plus2"_P38_"\$rep".xvg -ref Chain_"\$target_chain"_Y169 -sel Chain_"\$interacting_chain_plus2"_P38
    gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_INTERPROTOMER_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/INTERPROTOMER/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_E167_Chain_"\$interacting_chain_plus2"_Q49_"\$rep".xvg -ref Chain_"\$target_chain"_E167 -sel Chain_"\$interacting_chain_plus2"_Q49
  fi

  ######NUCLEOTIDE DISTANCES BETWEEN NUCLEOTIDE/PO4 AND NEARBY RESIDUES
  #nucleotide index
  target_S14=\$((14 + number_of_residues * (target_chain_number - 1)))
  target_G74=\$((74 + number_of_residues * (target_chain_number - 1)))
  target_Q137=\$((137 + number_of_residues * (target_chain_number - 1)))
  target_V159=\$((159 + number_of_residues * (target_chain_number - 1)))
  target_D157=\$((157 + number_of_residues * (target_chain_number - 1)))
  target_K18=\$((18 + number_of_residues * (target_chain_number - 1)))
  target_K336=\$((336 + number_of_residues * (target_chain_number - 1)))
  target_M305=\$((305 + number_of_residues * (target_chain_number - 1)))
  target_Y306=\$((306 + number_of_residues * (target_chain_number - 1)))
  target_T303=\$((303 + number_of_residues * (target_chain_number - 1)))
  target_E214=\$((214 + number_of_residues * (target_chain_number - 1)))
  target_K213=\$((213 + number_of_residues * (target_chain_number - 1)))
  target_R210=\$((210 + number_of_residues * (target_chain_number - 1)))
  target_L16=\$((16 + number_of_residues * (target_chain_number - 1)))

  echo -e 'ri '\$nucleotide_index_number'\n name '\$start_index 'Chain_'\$target_chain'_'\$nucleotide'\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_NUCLEOTIDE_DISTANCE_INDEX.ndx -nobackup
  echo -e 'ri '\$target_S14'\n name '\$((start_index + 1)) 'Chain_'\$target_chain'_S14\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_NUCLEOTIDE_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_NUCLEOTIDE_DISTANCE_INDEX.ndx -nobackup
  echo -e 'ri '\$target_G74'\n name '\$((start_index + 2)) 'Chain_'\$target_chain'_G74\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_NUCLEOTIDE_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_NUCLEOTIDE_DISTANCE_INDEX.ndx -nobackup
  echo -e 'ri '\$target_Q137'\n name '\$((start_index + 3)) 'Chain_'\$target_chain'_Q137\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_NUCLEOTIDE_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_NUCLEOTIDE_DISTANCE_INDEX.ndx -nobackup
  echo -e 'ri '\$target_V159'\n name '\$((start_index + 4)) 'Chain_'\$target_chain'_V159\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_NUCLEOTIDE_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_NUCLEOTIDE_DISTANCE_INDEX.ndx -nobackup
  echo -e 'ri '\$target_D157'\n name '\$((start_index + 5)) 'Chain_'\$target_chain'_D157\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_NUCLEOTIDE_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_NUCLEOTIDE_DISTANCE_INDEX.ndx -nobackup
  echo -e 'ri '\$target_K18'\n name '\$((start_index + 6)) 'Chain_'\$target_chain'_K18\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_NUCLEOTIDE_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_NUCLEOTIDE_DISTANCE_INDEX.ndx -nobackup
  echo -e 'ri '\$target_K336'\n name '\$((start_index + 7)) 'Chain_'\$target_chain'_K336\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_NUCLEOTIDE_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_NUCLEOTIDE_DISTANCE_INDEX.ndx -nobackup
  echo -e 'ri '\$target_M305'\n name '\$((start_index + 8)) 'Chain_'\$target_chain'_M305\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_NUCLEOTIDE_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_NUCLEOTIDE_DISTANCE_INDEX.ndx -nobackup
  echo -e 'ri '\$target_Y306'\n name '\$((start_index + 9)) 'Chain_'\$target_chain'_Y306\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_NUCLEOTIDE_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_NUCLEOTIDE_DISTANCE_INDEX.ndx -nobackup
  echo -e 'ri '\$target_T303'\n name '\$((start_index + 10)) 'Chain_'\$target_chain'_T303\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_NUCLEOTIDE_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_NUCLEOTIDE_DISTANCE_INDEX.ndx -nobackup
  echo -e 'ri '\$target_E214'\n name '\$((start_index + 11)) 'Chain_'\$target_chain'_E214\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_NUCLEOTIDE_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_NUCLEOTIDE_DISTANCE_INDEX.ndx -nobackup
  echo -e 'ri '\$target_K213'\n name '\$((start_index + 12)) 'Chain_'\$target_chain'_K213\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_NUCLEOTIDE_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_NUCLEOTIDE_DISTANCE_INDEX.ndx -nobackup
  echo -e 'ri '\$target_R210'\n name '\$((start_index + 13)) 'Chain_'\$target_chain'_R210\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_NUCLEOTIDE_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_NUCLEOTIDE_DISTANCE_INDEX.ndx -nobackup
  echo -e 'ri '\$target_L16'\n name '\$((start_index + 14)) 'Chain_'\$target_chain'_L16\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_NUCLEOTIDE_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_NUCLEOTIDE_DISTANCE_INDEX.ndx -nobackup

  gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_NUCLEOTIDE_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/NUCLEOTIDE/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_S14_Chain_"\$target_chain"_G74_"\$rep".xvg -ref Chain_"\$target_chain"_S14 -sel Chain_"\$target_chain"_G74
  gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_NUCLEOTIDE_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/NUCLEOTIDE/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_"\$nucleotide"_Chain_"\$target_chain"_K18_"\$rep".xvg -ref Chain_"\$target_chain"_"\$nucleotide" -sel Chain_"\$target_chain"_K18
  gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_NUCLEOTIDE_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/NUCLEOTIDE/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_"\$nucleotide"_Chain_"\$target_chain"_K336_"\$rep".xvg -ref Chain_"\$target_chain"_"\$nucleotide" -sel Chain_"\$target_chain"_K336
  gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_NUCLEOTIDE_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/NUCLEOTIDE/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_"\$nucleotide"_Chain_"\$target_chain"_M305_"\$rep".xvg -ref Chain_"\$target_chain"_"\$nucleotide" -sel Chain_"\$target_chain"_M305
  gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_NUCLEOTIDE_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/NUCLEOTIDE/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_"\$nucleotide"_Chain_"\$target_chain"_Y306_"\$rep".xvg -ref Chain_"\$target_chain"_"\$nucleotide" -sel Chain_"\$target_chain"_Y306
  gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_NUCLEOTIDE_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/NUCLEOTIDE/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_"\$nucleotide"_Chain_"\$target_chain"_T303_"\$rep".xvg -ref Chain_"\$target_chain"_"\$nucleotide" -sel Chain_"\$target_chain"_T303
  gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_NUCLEOTIDE_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/NUCLEOTIDE/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_"\$nucleotide"_Chain_"\$target_chain"_E214_"\$rep".xvg -ref Chain_"\$target_chain"_"\$nucleotide" -sel Chain_"\$target_chain"_E214
  gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_NUCLEOTIDE_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/NUCLEOTIDE/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_"\$nucleotide"_Chain_"\$target_chain"_K213_"\$rep".xvg -ref Chain_"\$target_chain"_"\$nucleotide" -sel Chain_"\$target_chain"_K213
  gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_NUCLEOTIDE_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/NUCLEOTIDE/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_"\$nucleotide"_Chain_"\$target_chain"_R210_"\$rep".xvg -ref Chain_"\$target_chain"_"\$nucleotide" -sel Chain_"\$target_chain"_R210
  gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_NUCLEOTIDE_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/NUCLEOTIDE/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_"\$nucleotide"_Chain_"\$target_chain"_D157_"\$rep".xvg -ref Chain_"\$target_chain"_"\$nucleotide" -sel Chain_"\$target_chain"_D157
  gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_NUCLEOTIDE_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/NUCLEOTIDE/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_"\$nucleotide"_Chain_"\$target_chain"_L16_"\$rep".xvg -ref Chain_"\$target_chain"_"\$nucleotide" -sel Chain_"\$target_chain"_L16

  if [ "\$nucleotide_state" == 'ADP-Pi' ]; then
    echo -e 'ri '\$po4_index_number'\n name '\$((start_index + 15)) 'Chain_'\$target_chain'_PO4\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_NUCLEOTIDE_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_NUCLEOTIDE_DISTANCE_INDEX.ndx -nobackup

    gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_NUCLEOTIDE_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/NUCLEOTIDE/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_PO4_Chain_"\$target_chain"_Q137_"\$rep".xvg -ref Chain_"\$target_chain"_PO4 -sel Chain_"\$target_chain"_Q137
    gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_NUCLEOTIDE_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/NUCLEOTIDE/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_PO4_Chain_"\$target_chain"_V159_"\$rep".xvg -ref Chain_"\$target_chain"_PO4 -sel Chain_"\$target_chain"_V159
    gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_NUCLEOTIDE_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/NUCLEOTIDE/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_PO4_Chain_"\$target_chain"_S14_"\$rep".xvg -ref Chain_"\$target_chain"_PO4 -sel Chain_"\$target_chain"_S14
    gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_NUCLEOTIDE_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/NUCLEOTIDE/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_PO4_Chain_"\$target_chain"_D157_"\$rep".xvg -ref Chain_"\$target_chain"_PO4 -sel Chain_"\$target_chain"_D157
  fi

  #INTERSTRAND DISTANCES. Looking at interstrand distances near/around interprotomer interface
  if [ "\$target_chain_number" -gt 2 ]; then
    interacting_minus2_K113=\$((113 + (number_of_residues * (target_chain_number - 3))))
    interacting_minus2_H173=\$((173 + (number_of_residues * (target_chain_number - 3))))
    interacting_minus2_P112=\$((112 + (number_of_residues * (target_chain_number - 3))))
    interacting_minus1_E195=\$((195 + (number_of_residues * (target_chain_number - 2))))
    interacting_minus1_E270=\$((270 + (number_of_residues * (target_chain_number - 2))))
    interacting_minus1_G268=\$((268 + (number_of_residues * (target_chain_number - 2))))
    interacting_minus1_I267=\$((267 + (number_of_residues * (target_chain_number - 2))))
    interacting_minus1_T194=\$((194 + (number_of_residues * (target_chain_number - 2))))
    target_T202=\$((202 + (number_of_residues * (target_chain_number - 1))))
    target_T203=\$((203 + (number_of_residues * (target_chain_number - 1))))
    target_T66=\$((66 + (number_of_residues * (target_chain_number - 1))))
    target_R39=\$((39 + (number_of_residues * (target_chain_number - 1))))
    target_K113=\$((113 + (number_of_residues * (target_chain_number - 1))))
    target_H173=\$((173 + (number_of_residues * (target_chain_number - 1))))
    target_P112=\$((112 + (number_of_residues * (target_chain_number - 1))))
    interacting_plus1_E195=\$((195 + (number_of_residues * (target_chain_number))))
    interacting_plus1_E270=\$((270 + (number_of_residues * (target_chain_number))))
    interacting_plus1_G268=\$((268 + (number_of_residues * (target_chain_number))))
    interacting_plus1_I267=\$((267 + (number_of_residues * (target_chain_number))))
    interacting_plus1_T194=\$((194 + (number_of_residues * (target_chain_number))))
    interacting_plus2_T202=\$((202 + (number_of_residues * (target_chain_number + 1))))
    interacting_plus2_T203=\$((203 + (number_of_residues * (target_chain_number + 1))))
    interacting_plus2_T66=\$((66 + (number_of_residues * (target_chain_number + 1))))
    interacting_plus2_R39=\$((39 + (number_of_residues * (target_chain_number + 1))))

    echo -e 'ri '\$interacting_minus2_K113'\n name '\$start_index 'Chain_'\$interacting_chain_minus2'_K113\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$interacting_minus2_H173'\n name '\$((start_index + 1)) 'Chain_'\$interacting_chain_minus2'_H173\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$interacting_minus2_P112'\n name '\$((start_index + 2)) 'Chain_'\$interacting_chain_minus2'_P112\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$interacting_minus1_E195'\n name '\$((start_index + 3)) 'Chain_'\$interacting_chain_minus1'_E195\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$interacting_minus1_E270'\n name '\$((start_index + 4)) 'Chain_'\$interacting_chain_minus1'_E270\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$interacting_minus1_G268'\n name '\$((start_index + 5)) 'Chain_'\$interacting_chain_minus1'_G268\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$interacting_minus1_I267'\n name '\$((start_index + 6)) 'Chain_'\$interacting_chain_minus1'_I267\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$interacting_minus1_T194'\n name '\$((start_index + 7)) 'Chain_'\$interacting_chain_minus1'_T194\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$target_T202'\n name '\$((start_index + 8)) 'Chain_'\$target_chain'_T202\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$target_T203'\n name '\$((start_index + 9)) 'Chain_'\$target_chain'_T203\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$target_T66'\n name '\$((start_index + 10)) 'Chain_'\$target_chain'_T66\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$target_R39'\n name '\$((start_index + 11)) 'Chain_'\$target_chain'_R39\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$target_K113'\n name '\$((start_index + 12)) 'Chain_'\$target_chain'_K113\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$target_H173'\n name '\$((start_index + 13)) 'Chain_'\$target_chain'_H173\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$target_P112'\n name '\$((start_index + 14)) 'Chain_'\$target_chain'_P112\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$interacting_plus1_E195'\n name '\$((start_index + 15)) 'Chain_'\$interacting_chain_plus1'_E195\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$interacting_plus1_E270'\n name '\$((start_index + 16)) 'Chain_'\$interacting_chain_plus1'_E270\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$interacting_plus1_G268'\n name '\$((start_index + 17)) 'Chain_'\$interacting_chain_plus1'_G268\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$interacting_plus1_I267'\n name '\$((start_index + 18)) 'Chain_'\$interacting_chain_plus1'_I267\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$interacting_plus1_T194'\n name '\$((start_index + 19)) 'Chain_'\$interacting_chain_plus1'_T194\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$interacting_plus2_T202'\n name '\$((start_index + 20)) 'Chain_'\$interacting_chain_plus2'_T202\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$interacting_plus2_T203'\n name '\$((start_index + 21)) 'Chain_'\$interacting_chain_plus2'_T203\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$interacting_plus2_T66'\n name '\$((start_index + 22)) 'Chain_'\$interacting_chain_plus2'_T66\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -nobackup
    echo -e 'ri '\$interacting_plus2_R39'\n name '\$((start_index + 23)) 'Chain_'\$interacting_chain_plus2'_R39\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -nobackup

    gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/INTERSTRAND/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$interacting_chain_minus2"_K113_Chain_"\$interacting_chain_minus1"_E195_"\$rep".xvg -ref Chain_"\$interacting_chain_minus2"_K113 -sel Chain_"\$interacting_chain_minus1"_E195
    gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/INTERSTRAND/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$interacting_chain_minus1"_E270_Chain_"\$target_chain"_T202_"\$rep".xvg -ref Chain_"\$interacting_chain_minus1"_E270 -sel Chain_"\$target_chain"_T202
    gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/INTERSTRAND/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$interacting_chain_minus1"_E270_Chain_"\$target_chain"_T203_"\$rep".xvg -ref Chain_"\$interacting_chain_minus1"_E270 -sel Chain_"\$target_chain"_T203
    gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/INTERSTRAND/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$interacting_chain_minus1"_E270_Chain_"\$target_chain"_T66_"\$rep".xvg -ref Chain_"\$interacting_chain_minus1"_E270 -sel Chain_"\$target_chain"_T66
    gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/INTERSTRAND/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$interacting_chain_minus1"_E270_Chain_"\$target_chain"_R39_"\$rep".xvg -ref Chain_"\$interacting_chain_minus1"_E270 -sel Chain_"\$target_chain"_R39
    gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/INTERSTRAND/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$interacting_chain_minus1"_G268_Chain_"\$target_chain"_R39_"\$rep".xvg -ref Chain_"\$interacting_chain_minus1"_G268 -sel Chain_"\$target_chain"_R39
    gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/INTERSTRAND/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$interacting_chain_minus2"_H173_Chain_"\$interacting_chain_minus1"_G268_"\$rep".xvg -ref Chain_"\$interacting_chain_minus2"_H173 -sel Chain_"\$interacting_chain_minus1"_G268
    gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/INTERSTRAND/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$interacting_chain_minus2"_H173_Chain_"\$interacting_chain_minus1"_I267_"\$rep".xvg -ref Chain_"\$interacting_chain_minus2"_H173 -sel Chain_"\$interacting_chain_minus1"_I267
    gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/INTERSTRAND/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$interacting_chain_minus2"_P112_Chain_"\$interacting_chain_minus1"_T194_"\$rep".xvg -ref Chain_"\$interacting_chain_minus2"_P112 -sel Chain_"\$interacting_chain_minus1"_T194

    gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/INTERSTRAND/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_K113_Chain_"\$interacting_chain_plus1"_E195_"\$rep".xvg -ref Chain_"\$target_chain"_K113 -sel Chain_"\$interacting_chain_plus1"_E195
    gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/INTERSTRAND/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$interacting_chain_plus1"_E270_Chain_"\$interacting_chain_plus2"_T202_"\$rep".xvg -ref Chain_"\$interacting_chain_plus1"_E270 -sel Chain_"\$interacting_chain_plus2"_T202
    gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/INTERSTRAND/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$interacting_chain_plus1"_E270_Chain_"\$interacting_chain_plus2"_T203_"\$rep".xvg -ref Chain_"\$interacting_chain_plus1"_E270 -sel Chain_"\$interacting_chain_plus2"_T203
    gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/INTERSTRAND/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$interacting_chain_plus1"_E270_Chain_"\$interacting_chain_plus2"_T66_"\$rep".xvg -ref Chain_"\$interacting_chain_plus1"_E270 -sel Chain_"\$interacting_chain_plus2"_T66
    gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/INTERSTRAND/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$interacting_chain_plus1"_E270_Chain_"\$interacting_chain_plus2"_R39_"\$rep".xvg -ref Chain_"\$interacting_chain_plus1"_E270 -sel Chain_"\$interacting_chain_plus2"_R39
    gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/INTERSTRAND/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$interacting_chain_plus1"_G268_Chain_"\$interacting_chain_plus2"_R39_"\$rep".xvg -ref Chain_"\$interacting_chain_plus1"_G268 -sel Chain_"\$interacting_chain_plus2"_R39
    gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/INTERSTRAND/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_H173_Chain_"\$interacting_chain_plus1"_G268_"\$rep".xvg -ref Chain_"\$target_chain"_H173 -sel Chain_"\$interacting_chain_plus1"_G268
    gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/INTERSTRAND/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_H173_Chain_"\$interacting_chain_plus1"_I267_"\$rep".xvg -ref Chain_"\$target_chain"_H173 -sel Chain_"\$interacting_chain_plus1"_I267
    gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_INTERSTRAND_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/INTERSTRAND/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_P112_Chain_"\$interacting_chain_plus1"_T194_"\$rep".xvg -ref Chain_"\$target_chain"_P112 -sel Chain_"\$interacting_chain_plus1"_T194
  fi



  #allosteric network index (network outlined in my literature review)
  target_F375=\$((number_of_residues * target_chain_number))
  target_R116=\$((116 + number_of_residues * (target_chain_number - 1)))
  target_H371=\$((371 + number_of_residues * (target_chain_number - 1)))
  target_E117=\$((117 + number_of_residues * (target_chain_number - 1)))
  target_E107=\$((107 + number_of_residues * (target_chain_number - 1)))
  target_Q137=\$((137 + number_of_residues * (target_chain_number - 1)))
  target_N115=\$((115 + number_of_residues * (target_chain_number - 1)))
  target_I76=\$((76 + number_of_residues * (target_chain_number - 1)))
  target_N111=\$((111 + number_of_residues * (target_chain_number - 1)))
  target_G74=\$((74 + number_of_residues * (target_chain_number - 1)))
  target_M123=\$((123 + number_of_residues * (target_chain_number - 1)))
  target_W86=\$((86 + number_of_residues * (target_chain_number - 1)))
  target_I122=\$((122 + number_of_residues * (target_chain_number - 1)))
  target_M119=\$((119 + number_of_residues * (target_chain_number - 1)))
  target_W79=\$((79 + number_of_residues * (target_chain_number - 1)))
  target_K118=\$((118 + number_of_residues * (target_chain_number - 1)))
  target_H88=\$((88 + number_of_residues * (target_chain_number - 1)))
  target_V54=\$((54 + number_of_residues * (target_chain_number - 1)))
  target_D56=\$((56 + number_of_residues * (target_chain_number - 1)))

  echo -e 'ri '\$target_F375'\n name '\$start_index 'Chain_'\$target_chain'_F375\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_ALLOSTERY_DISTANCE_INDEX.ndx -nobackup
  echo -e 'ri '\$target_R116'\n name '\$((start_index + 1)) 'Chain_'\$target_chain'_R116\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_ALLOSTERY_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_ALLOSTERY_DISTANCE_INDEX.ndx -nobackup
  echo -e 'ri '\$target_H371'\n name '\$((start_index + 2)) 'Chain_'\$target_chain'_H371\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_ALLOSTERY_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_ALLOSTERY_DISTANCE_INDEX.ndx -nobackup
  echo -e 'ri '\$target_E117'\n name '\$((start_index + 3)) 'Chain_'\$target_chain'_E117\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_ALLOSTERY_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_ALLOSTERY_DISTANCE_INDEX.ndx -nobackup
  echo -e 'ri '\$target_E107'\n name '\$((start_index + 4)) 'Chain_'\$target_chain'_E107\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_ALLOSTERY_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_ALLOSTERY_DISTANCE_INDEX.ndx -nobackup
  echo -e 'ri '\$target_Q137'\n name '\$((start_index + 5)) 'Chain_'\$target_chain'_Q137\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_ALLOSTERY_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_ALLOSTERY_DISTANCE_INDEX.ndx -nobackup
  echo -e 'ri '\$target_N115'\n name '\$((start_index + 6)) 'Chain_'\$target_chain'_N115\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_ALLOSTERY_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_ALLOSTERY_DISTANCE_INDEX.ndx -nobackup
  echo -e 'ri '\$target_I76'\n name '\$((start_index + 7)) 'Chain_'\$target_chain'_I76\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_ALLOSTERY_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_ALLOSTERY_DISTANCE_INDEX.ndx -nobackup
  echo -e 'ri '\$target_N111'\n name '\$((start_index + 8)) 'Chain_'\$target_chain'_N111\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_ALLOSTERY_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_ALLOSTERY_DISTANCE_INDEX.ndx -nobackup
  echo -e 'ri '\$target_G74'\n name '\$((start_index + 9)) 'Chain_'\$target_chain'_G74\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_ALLOSTERY_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_ALLOSTERY_DISTANCE_INDEX.ndx -nobackup
  echo -e 'ri '\$target_M123'\n name '\$((start_index + 10)) 'Chain_'\$target_chain'_M123\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_ALLOSTERY_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_ALLOSTERY_DISTANCE_INDEX.ndx -nobackup
  echo -e 'ri '\$target_W86'\n name '\$((start_index + 11)) 'Chain_'\$target_chain'_W86\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_ALLOSTERY_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_ALLOSTERY_DISTANCE_INDEX.ndx -nobackup
  echo -e 'ri '\$target_I122'\n name '\$((start_index + 12)) 'Chain_'\$target_chain'_I122\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_ALLOSTERY_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_ALLOSTERY_DISTANCE_INDEX.ndx -nobackup
  echo -e 'ri '\$target_M119'\n name '\$((start_index + 13)) 'Chain_'\$target_chain'_M119\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_ALLOSTERY_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_ALLOSTERY_DISTANCE_INDEX.ndx -nobackup
  echo -e 'ri '\$target_W79'\n name '\$((start_index + 14)) 'Chain_'\$target_chain'_W79\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_ALLOSTERY_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_ALLOSTERY_DISTANCE_INDEX.ndx -nobackup
  echo -e 'ri '\$target_K118'\n name '\$((start_index + 15)) 'Chain_'\$target_chain'_K118\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_ALLOSTERY_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_ALLOSTERY_DISTANCE_INDEX.ndx -nobackup
  echo -e 'ri '\$target_H88'\n name '\$((start_index + 16)) 'Chain_'\$target_chain'_H88\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_ALLOSTERY_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_ALLOSTERY_DISTANCE_INDEX.ndx -nobackup
  echo -e 'ri '\$target_V54'\n name '\$((start_index + 17)) 'Chain_'\$target_chain'_V54\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_ALLOSTERY_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_ALLOSTERY_DISTANCE_INDEX.ndx -nobackup
  echo -e 'ri '\$target_D56'\n name '\$((start_index + 18)) 'Chain_'\$target_chain'_D56\n q' | gmx make_ndx -f "\$mutant"_"\$Factin_or_Gactin"_FIRST_FRAME.pdb -o "\$mutant"_"\$Factin_or_Gactin"_ALLOSTERY_DISTANCE_INDEX.ndx -n "\$mutant"_"\$Factin_or_Gactin"_ALLOSTERY_DISTANCE_INDEX.ndx -nobackup

  gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_ALLOSTERY_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/ALLOSTERY/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_F375_Chain_"\$target_chain"_R116_"\$rep".xvg -ref Chain_"\$target_chain"_F375 -sel Chain_"\$target_chain"_R116
  gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_ALLOSTERY_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/ALLOSTERY/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_H371_Chain_"\$target_chain"_E117_"\$rep".xvg -ref Chain_"\$target_chain"_H371 -sel Chain_"\$target_chain"_E117
  gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_ALLOSTERY_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/ALLOSTERY/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_R116_Chain_"\$target_chain"_E107_"\$rep".xvg -ref Chain_"\$target_chain"_R116 -sel Chain_"\$target_chain"_E107
  gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_ALLOSTERY_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/ALLOSTERY/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_E107_Chain_"\$target_chain"_Q137_"\$rep".xvg -ref Chain_"\$target_chain"_E107 -sel Chain_"\$target_chain"_Q137
  gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_ALLOSTERY_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/ALLOSTERY/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_N115_Chain_"\$target_chain"_I76_"\$rep".xvg -ref Chain_"\$target_chain"_N115 -sel Chain_"\$target_chain"_I76
  gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_ALLOSTERY_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/ALLOSTERY/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_N111_Chain_"\$target_chain"_G74_"\$rep".xvg -ref Chain_"\$target_chain"_N111 -sel Chain_"\$target_chain"_G74
  gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_ALLOSTERY_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/ALLOSTERY/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_M123_Chain_"\$target_chain"_W86_"\$rep".xvg -ref Chain_"\$target_chain"_M123 -sel Chain_"\$target_chain"_W86
  gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_ALLOSTERY_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/ALLOSTERY/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_I122_Chain_"\$target_chain"_W86_"\$rep".xvg -ref Chain_"\$target_chain"_I122 -sel Chain_"\$target_chain"_W86
  gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_ALLOSTERY_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/ALLOSTERY/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_M119_Chain_"\$target_chain"_W79_"\$rep".xvg -ref Chain_"\$target_chain"_M119 -sel Chain_"\$target_chain"_W79
  gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_ALLOSTERY_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/ALLOSTERY/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_K118_Chain_"\$target_chain"_W79_"\$rep".xvg -ref Chain_"\$target_chain"_K118 -sel Chain_"\$target_chain"_W79
  gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_ALLOSTERY_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/ALLOSTERY/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_I76_Chain_"\$target_chain"_W79_"\$rep".xvg -ref Chain_"\$target_chain"_I76 -sel Chain_"\$target_chain"_W79
  gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_ALLOSTERY_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/ALLOSTERY/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_H88_Chain_"\$target_chain"_V54_"\$rep".xvg -ref Chain_"\$target_chain"_H88 -sel Chain_"\$target_chain"_V54
  gmx pairdist -n "\$mutant"_"\$Factin_or_Gactin"_ALLOSTERY_DISTANCE_INDEX.ndx -f "\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -o DISTANCES/ALLOSTERY/"\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_H88_Chain_"\$target_chain"_D56_"\$rep".xvg -ref Chain_"\$target_chain"_H88 -sel Chain_"\$target_chain"_D56

done


#########COMBINE TRAJECTORIES FOR OTHER ANALYSES

cd "\$directory_path"/"\$mutant"/

mkdir "\$directory_path"/"\$mutant"/COMBINED

cd "\$directory_path"/"\$mutant"/COMBINED/

#Combine trajectories for subsequent analyses, like network analysis
if [ "\$nucleotide_state" == "ADP-Pi" ]; then
  echo -e 'Protein_'\$cation'_ADP_PO4\nc\nc\nc\n' | gmx trjcat -f ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc ../Replicate2/"\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc ../Replicate3/"\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_index.ndx -o "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -settime
else
  echo -e 'Protein_'\$cation'_'\$nucleotide'\nc\nc\nc\n' | gmx trjcat -f ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc ../Replicate2/"\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc ../Replicate3/"\$mutant"_"\$Factin_or_Gactin"_md_trajec_fit.xtc -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_index.ndx -o "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -settime
fi

echo -e 'Protein\n Protein' | gmx trjconv -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -o "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.pdb -fit rot+trans -dump 0

echo -e 'Chain_'\$target_chain'\n Chain_'\$target_chain | gmx trjconv -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -o "\$mutant"_"\$Factin_or_Gactin"_"\$target_chain"_COMBINED.xtc -fit rot+trans

echo -e 'Chain_'\$target_chain'_protein\n Chain_'\$target_chain'_protein' | gmx trjconv -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -o "\$mutant"_"\$Factin_or_Gactin"_"\$target_chain"_protein_COMBINED.xtc -fit rot+trans
echo -e 'Chain_'\$target_chain'_protein\n Chain_'\$target_chain'_protein' | gmx trjconv -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -o "\$mutant"_"\$Factin_or_Gactin"_"\$target_chain"_protein_COMBINED.pdb -fit rot+trans -dump 0

#cluster analysis
echo -e 'C-alpha\n Chain_'\$target_chain | gmx cluster -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -method gromos -cutoff 0.15 -g "\$mutant"_"\$Factin_or_Gactin"_cluster_rep1.log -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -b 1 -e 200 -tu ns
echo -e 'C-alpha\n Chain_'\$target_chain | gmx cluster -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -method gromos -cutoff 0.15 -g "\$mutant"_"\$Factin_or_Gactin"_cluster_rep2.log -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -b 201 -e 400 -tu ns
echo -e 'C-alpha\n Chain_'\$target_chain | gmx cluster -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -method gromos -cutoff 0.15 -g "\$mutant"_"\$Factin_or_Gactin"_cluster_rep3.log -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -b 401 -e 600 -tu ns
echo -e 'C-alpha\n Chain_'\$target_chain | gmx cluster -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -method gromos -cutoff 0.15 -g "\$mutant"_"\$Factin_or_Gactin"_cluster_COMBINED.log -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -b 1 -e 600 -tu ns


#pull out the timepoints associated with each cluster, then separate clusters into separate variables
rep1_cluster_values=\$(awk -F '|' '/[[:space:]]*[0-9]+[[:space:]]*\|[[:space:]]*[0-9]+[[:space:]]*[0-9]+\.[0-9]+[[:space:]]*\|[[:space:]]*[0-9]+(\.[0-9]+)?[[:space:]]*\.[0-9]+[[:space:]]*\|/ {print \$3}' "\$mutant"_"\$Factin_or_Gactin"_cluster_rep1.log)
rep1_cluster1=\$(echo "\$rep1_cluster_values" | awk 'NR==1{print \$1}')
rep1_cluster2=\$(echo "\$rep1_cluster_values" | awk 'NR==2{print \$1}')
rep1_cluster3=\$(echo "\$rep1_cluster_values" | awk 'NR==3{print \$1}')
rep1_cluster4=\$(echo "\$rep1_cluster_values" | awk 'NR==4{print \$1}')
rep1_cluster5=\$(echo "\$rep1_cluster_values" | awk 'NR==5{print \$1}')
rep1_cluster6=\$(echo "\$rep1_cluster_values" | awk 'NR==6{print \$1}')
rep1_cluster7=\$(echo "\$rep1_cluster_values" | awk 'NR==7{print \$1}')
rep1_cluster8=\$(echo "\$rep1_cluster_values" | awk 'NR==8{print \$1}')
rep1_cluster9=\$(echo "\$rep1_cluster_values" | awk 'NR==9{print \$1}')
rep1_cluster10=\$(echo "\$rep1_cluster_values" | awk 'NR==10{print \$1}')
echo -e 'Chain_'\$target_chain'\n Chain_'\$target_chain | gmx trjconv -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -o "\$mutant"_rep1_cluster1.pdb -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -fit rot+trans -tu ns -dump "\$rep1_cluster1"
echo -e 'Chain_'\$target_chain'\n Chain_'\$target_chain | gmx trjconv -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -o "\$mutant"_rep1_cluster2.pdb -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -fit rot+trans -tu ns -dump "\$rep1_cluster2"
echo -e 'Chain_'\$target_chain'\n Chain_'\$target_chain | gmx trjconv -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -o "\$mutant"_rep1_cluster3.pdb -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -fit rot+trans -tu ns -dump "\$rep1_cluster3"
echo -e 'Chain_'\$target_chain'\n Chain_'\$target_chain | gmx trjconv -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -o "\$mutant"_rep1_cluster4.pdb -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -fit rot+trans -tu ns -dump "\$rep1_cluster4"
echo -e 'Chain_'\$target_chain'\n Chain_'\$target_chain | gmx trjconv -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -o "\$mutant"_rep1_cluster5.pdb -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -fit rot+trans -tu ns -dump "\$rep1_cluster5"
echo -e 'Chain_'\$target_chain'\n Chain_'\$target_chain | gmx trjconv -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -o "\$mutant"_rep1_cluster6.pdb -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -fit rot+trans -tu ns -dump "\$rep1_cluster6"
echo -e 'Chain_'\$target_chain'\n Chain_'\$target_chain | gmx trjconv -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -o "\$mutant"_rep1_cluster7.pdb -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -fit rot+trans -tu ns -dump "\$rep1_cluster7"
echo -e 'Chain_'\$target_chain'\n Chain_'\$target_chain | gmx trjconv -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -o "\$mutant"_rep1_cluster8.pdb -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -fit rot+trans -tu ns -dump "\$rep1_cluster8"
echo -e 'Chain_'\$target_chain'\n Chain_'\$target_chain | gmx trjconv -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -o "\$mutant"_rep1_cluster9.pdb -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -fit rot+trans -tu ns -dump "\$rep1_cluster9"
echo -e 'Chain_'\$target_chain'\n Chain_'\$target_chain | gmx trjconv -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -o "\$mutant"_rep1_cluster10.pdb -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -fit rot+trans -tu ns -dump "\$rep1_cluster10"

rep2_cluster_values=\$(awk -F '|' '/[[:space:]]*[0-9]+[[:space:]]*\|[[:space:]]*[0-9]+[[:space:]]*[0-9]+\.[0-9]+[[:space:]]*\|[[:space:]]*[0-9]+(\.[0-9]+)?[[:space:]]*\.[0-9]+[[:space:]]*\|/ {print \$3}' "\$mutant"_"\$Factin_or_Gactin"_cluster_rep2.log)
rep2_cluster1=\$(echo "\$rep2_cluster_values" | awk 'NR==1{print \$1}')
rep2_cluster2=\$(echo "\$rep2_cluster_values" | awk 'NR==2{print \$1}')
rep2_cluster3=\$(echo "\$rep2_cluster_values" | awk 'NR==3{print \$1}')
rep2_cluster4=\$(echo "\$rep2_cluster_values" | awk 'NR==4{print \$1}')
rep2_cluster5=\$(echo "\$rep2_cluster_values" | awk 'NR==5{print \$1}')
rep2_cluster6=\$(echo "\$rep2_cluster_values" | awk 'NR==6{print \$1}')
rep2_cluster7=\$(echo "\$rep2_cluster_values" | awk 'NR==7{print \$1}')
rep2_cluster8=\$(echo "\$rep2_cluster_values" | awk 'NR==8{print \$1}')
rep2_cluster9=\$(echo "\$rep2_cluster_values" | awk 'NR==9{print \$1}')
rep2_cluster10=\$(echo "\$rep2_cluster_values" | awk 'NR==10{print \$1}')
echo -e 'Chain_'\$target_chain'\n Chain_'\$target_chain | gmx trjconv -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -o "\$mutant"_rep2_cluster1.pdb -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -fit rot+trans -tu ns -dump "\$rep2_cluster1"
echo -e 'Chain_'\$target_chain'\n Chain_'\$target_chain | gmx trjconv -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -o "\$mutant"_rep2_cluster2.pdb -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -fit rot+trans -tu ns -dump "\$rep2_cluster2"
echo -e 'Chain_'\$target_chain'\n Chain_'\$target_chain | gmx trjconv -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -o "\$mutant"_rep2_cluster3.pdb -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -fit rot+trans -tu ns -dump "\$rep2_cluster3"
echo -e 'Chain_'\$target_chain'\n Chain_'\$target_chain | gmx trjconv -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -o "\$mutant"_rep2_cluster4.pdb -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -fit rot+trans -tu ns -dump "\$rep2_cluster4"
echo -e 'Chain_'\$target_chain'\n Chain_'\$target_chain | gmx trjconv -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -o "\$mutant"_rep2_cluster5.pdb -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -fit rot+trans -tu ns -dump "\$rep2_cluster5"
echo -e 'Chain_'\$target_chain'\n Chain_'\$target_chain | gmx trjconv -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -o "\$mutant"_rep2_cluster6.pdb -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -fit rot+trans -tu ns -dump "\$rep2_cluster6"
echo -e 'Chain_'\$target_chain'\n Chain_'\$target_chain | gmx trjconv -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -o "\$mutant"_rep2_cluster7.pdb -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -fit rot+trans -tu ns -dump "\$rep2_cluster7"
echo -e 'Chain_'\$target_chain'\n Chain_'\$target_chain | gmx trjconv -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -o "\$mutant"_rep2_cluster8.pdb -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -fit rot+trans -tu ns -dump "\$rep2_cluster8"
echo -e 'Chain_'\$target_chain'\n Chain_'\$target_chain | gmx trjconv -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -o "\$mutant"_rep2_cluster9.pdb -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -fit rot+trans -tu ns -dump "\$rep2_cluster9"
echo -e 'Chain_'\$target_chain'\n Chain_'\$target_chain | gmx trjconv -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -o "\$mutant"_rep2_cluster10.pdb -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -fit rot+trans -tu ns -dump "\$rep2_cluster10"

rep3_cluster_values=\$(awk -F '|' '/[[:space:]]*[0-9]+[[:space:]]*\|[[:space:]]*[0-9]+[[:space:]]*[0-9]+\.[0-9]+[[:space:]]*\|[[:space:]]*[0-9]+(\.[0-9]+)?[[:space:]]*\.[0-9]+[[:space:]]*\|/ {print \$3}' "\$mutant"_"\$Factin_or_Gactin"_cluster_rep3.log)
rep3_cluster1=\$(echo "\$rep3_cluster_values" | awk 'NR==1{print \$1}')
rep3_cluster2=\$(echo "\$rep3_cluster_values" | awk 'NR==2{print \$1}')
rep3_cluster3=\$(echo "\$rep3_cluster_values" | awk 'NR==3{print \$1}')
rep3_cluster4=\$(echo "\$rep3_cluster_values" | awk 'NR==4{print \$1}')
rep3_cluster5=\$(echo "\$rep3_cluster_values" | awk 'NR==5{print \$1}')
rep3_cluster6=\$(echo "\$rep3_cluster_values" | awk 'NR==6{print \$1}')
rep3_cluster7=\$(echo "\$rep3_cluster_values" | awk 'NR==7{print \$1}')
rep3_cluster8=\$(echo "\$rep3_cluster_values" | awk 'NR==8{print \$1}')
rep3_cluster9=\$(echo "\$rep3_cluster_values" | awk 'NR==9{print \$1}')
rep3_cluster10=\$(echo "\$rep3_cluster_values" | awk 'NR==10{print \$1}')
echo -e 'Chain_'\$target_chain'\n Chain_'\$target_chain | gmx trjconv -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -o "\$mutant"_rep3_cluster1.pdb -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -fit rot+trans -tu ns -dump "\$rep3_cluster1"
echo -e 'Chain_'\$target_chain'\n Chain_'\$target_chain | gmx trjconv -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -o "\$mutant"_rep3_cluster2.pdb -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -fit rot+trans -tu ns -dump "\$rep3_cluster2"
echo -e 'Chain_'\$target_chain'\n Chain_'\$target_chain | gmx trjconv -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -o "\$mutant"_rep3_cluster3.pdb -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -fit rot+trans -tu ns -dump "\$rep3_cluster3"
echo -e 'Chain_'\$target_chain'\n Chain_'\$target_chain | gmx trjconv -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -o "\$mutant"_rep3_cluster4.pdb -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -fit rot+trans -tu ns -dump "\$rep3_cluster4"
echo -e 'Chain_'\$target_chain'\n Chain_'\$target_chain | gmx trjconv -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -o "\$mutant"_rep3_cluster5.pdb -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -fit rot+trans -tu ns -dump "\$rep3_cluster5"
echo -e 'Chain_'\$target_chain'\n Chain_'\$target_chain | gmx trjconv -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -o "\$mutant"_rep3_cluster6.pdb -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -fit rot+trans -tu ns -dump "\$rep3_cluster6"
echo -e 'Chain_'\$target_chain'\n Chain_'\$target_chain | gmx trjconv -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -o "\$mutant"_rep3_cluster7.pdb -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -fit rot+trans -tu ns -dump "\$rep3_cluster7"
echo -e 'Chain_'\$target_chain'\n Chain_'\$target_chain | gmx trjconv -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -o "\$mutant"_rep3_cluster8.pdb -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -fit rot+trans -tu ns -dump "\$rep3_cluster8"
echo -e 'Chain_'\$target_chain'\n Chain_'\$target_chain | gmx trjconv -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -o "\$mutant"_rep3_cluster9.pdb -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -fit rot+trans -tu ns -dump "\$rep3_cluster9"
echo -e 'Chain_'\$target_chain'\n Chain_'\$target_chain | gmx trjconv -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -o "\$mutant"_rep3_cluster10.pdb -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -fit rot+trans -tu ns -dump "\$rep3_cluster10"


combined_cluster_values=\$(awk -F '|' '/[[:space:]]*[0-9]+[[:space:]]*\|[[:space:]]*[0-9]+[[:space:]]*[0-9]+\.[0-9]+[[:space:]]*\|[[:space:]]*[0-9]+(\.[0-9]+)?[[:space:]]*\.[0-9]+[[:space:]]*\|/ {print \$3}' "\$mutant"_"\$Factin_or_Gactin"_cluster_COMBINED.log)
combined_cluster1=\$(echo "\$combined_cluster_values" | awk 'NR==1{print \$1}')
combined_cluster2=\$(echo "\$combined_cluster_values" | awk 'NR==2{print \$1}')
combined_cluster3=\$(echo "\$combined_cluster_values" | awk 'NR==3{print \$1}')
combined_cluster4=\$(echo "\$combined_cluster_values" | awk 'NR==4{print \$1}')
combined_cluster5=\$(echo "\$combined_cluster_values" | awk 'NR==5{print \$1}')
combined_cluster6=\$(echo "\$combined_cluster_values" | awk 'NR==6{print \$1}')
combined_cluster7=\$(echo "\$combined_cluster_values" | awk 'NR==7{print \$1}')
combined_cluster8=\$(echo "\$combined_cluster_values" | awk 'NR==8{print \$1}')
combined_cluster9=\$(echo "\$combined_cluster_values" | awk 'NR==9{print \$1}')
combined_cluster10=\$(echo "\$combined_cluster_values" | awk 'NR==10{print \$1}')
echo -e 'Chain_'\$target_chain'\n Chain_'\$target_chain | gmx trjconv -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -o "\$mutant"_combined_cluster1.pdb -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -fit rot+trans -tu ns -dump "\$combined_cluster1"
echo -e 'Chain_'\$target_chain'\n Chain_'\$target_chain | gmx trjconv -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -o "\$mutant"_combined_cluster2.pdb -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -fit rot+trans -tu ns -dump "\$combined_cluster2"
echo -e 'Chain_'\$target_chain'\n Chain_'\$target_chain | gmx trjconv -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -o "\$mutant"_combined_cluster3.pdb -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -fit rot+trans -tu ns -dump "\$combined_cluster3"
echo -e 'Chain_'\$target_chain'\n Chain_'\$target_chain | gmx trjconv -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -o "\$mutant"_combined_cluster4.pdb -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -fit rot+trans -tu ns -dump "\$combined_cluster4"
echo -e 'Chain_'\$target_chain'\n Chain_'\$target_chain | gmx trjconv -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -o "\$mutant"_combined_cluster5.pdb -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -fit rot+trans -tu ns -dump "\$combined_cluster5"
echo -e 'Chain_'\$target_chain'\n Chain_'\$target_chain | gmx trjconv -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -o "\$mutant"_combined_cluster6.pdb -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -fit rot+trans -tu ns -dump "\$combined_cluster6"
echo -e 'Chain_'\$target_chain'\n Chain_'\$target_chain | gmx trjconv -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -o "\$mutant"_combined_cluster7.pdb -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -fit rot+trans -tu ns -dump "\$combined_cluster7"
echo -e 'Chain_'\$target_chain'\n Chain_'\$target_chain | gmx trjconv -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -o "\$mutant"_combined_cluster8.pdb -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -fit rot+trans -tu ns -dump "\$combined_cluster8"
echo -e 'Chain_'\$target_chain'\n Chain_'\$target_chain | gmx trjconv -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -o "\$mutant"_combined_cluster9.pdb -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -fit rot+trans -tu ns -dump "\$combined_cluster9"
echo -e 'Chain_'\$target_chain'\n Chain_'\$target_chain | gmx trjconv -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -o "\$mutant"_combined_cluster10.pdb -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -fit rot+trans -tu ns -dump "\$combined_cluster10"


#PCA
#PCA whole filament along WT trajectory, then along it's own trajectory
#for PCA, use common eigenvectors for gibbs free energy landscape (ie project along WT vectors)
#for PCA RMSF and projections of extremes, use individual trajectories
if [ \$mutant == 'WT' ]; then

  echo -e 'C-alpha\n C-alpha' | gmx covar -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -v "\$mutant"_"\$Factin_or_Gactin"_eigenval_landscape_calpha.trr -o "\$mutant"_"\$Factin_or_Gactin"_eigenval_landscape_calpha.xvg
  echo -e 'C-alpha\n C-alpha' | gmx anaeig -v "\$mutant"_"\$Factin_or_Gactin"_eigenval_landscape_calpha.trr -eig "\$mutant"_"\$Factin_or_Gactin"_eigenval_landscape_calpha.xvg -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -extr "\$mutant"_"\$Factin_or_Gactin"_extreme_calpha.pdb -2d "\$mutant"_"\$Factin_or_Gactin"_2d_landscape_calpha.xvg -rmsf "\$mutant"_"\$Factin_or_Gactin"_eigenrmsf_landscape_calpha.xvg -first 1 -last 2 -nframes 10

  gmx sham -f "\$mutant"_"\$Factin_or_Gactin"_2d_landscape_calpha.xvg -ls "\$mutant"_"\$Factin_or_Gactin"_gibbs_landscape_calpha.xpm -notime -gmax 10
  gmx xpm2ps -f "\$mutant"_"\$Factin_or_Gactin"_gibbs_landscape_calpha.xpm -o "\$mutant"_"\$Factin_or_Gactin"_landscape_plot_calpha.eps -rainbow red

  #for EQ time points
  echo -e 'C-alpha\n C-alpha' | gmx covar -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -v "\$mutant"_"\$Factin_or_Gactin"_eigenval_landscape_calpha_EQ.trr -o "\$mutant"_"\$Factin_or_Gactin"_eigenval_landscape_calpha_EQ.xvg -b \$((eq_start_time*1000)) -e \$((simulation_length*1000))
  echo -e 'C-alpha\n C-alpha' | gmx anaeig -v "\$mutant"_"\$Factin_or_Gactin"_eigenval_landscape_calpha_EQ.trr -eig "\$mutant"_"\$Factin_or_Gactin"_eigenval_landscape_calpha_EQ.xvg -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -extr "\$mutant"_"\$Factin_or_Gactin"_extreme_calpha_EQ.pdb -2d "\$mutant"_"\$Factin_or_Gactin"_2d_landscape_calpha_EQ.xvg -rmsf "\$mutant"_"\$Factin_or_Gactin"_eigenrmsf_landscape_calpha_EQ.xvg -first 1 -last 2 -nframes 10 -b \$((eq_start_time*1000)) -e \$((simulation_length*1000))

  gmx sham -f "\$mutant"_"\$Factin_or_Gactin"_2d_landscape_calpha_EQ.xvg -ls "\$mutant"_"\$Factin_or_Gactin"_gibbs_landscape_calpha_EQ.xpm -notime -gmax 10
  gmx xpm2ps -f "\$mutant"_"\$Factin_or_Gactin"_gibbs_landscape_calpha_EQ.xpm -o "\$mutant"_"\$Factin_or_Gactin"_landscape_plot_calpha_EQ.eps -rainbow red

else

  #along WT trajectory
  echo -e 'C-alpha\n C-alpha' | gmx covar -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -v "\$mutant"_"\$Factin_or_Gactin"_eigenval_landscape_calpha.trr -o "\$mutant"_"\$Factin_or_Gactin"_eigenval_landscape_calpha.xvg
  echo -e 'C-alpha\n C-alpha' | gmx anaeig -v "\$directory_path"/WT/COMBINED/WT_"\$Factin_or_Gactin"_eigenval_landscape_calpha.trr -eig "\$directory_path"/WT/COMBINED/WT_"\$Factin_or_Gactin"_eigenval_landscape_calpha.xvg -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -extr "\$mutant"_"\$Factin_or_Gactin"_extreme_calpha.pdb -2d "\$mutant"_"\$Factin_or_Gactin"_2d_landscape_calpha.xvg -rmsf "\$mutant"_"\$Factin_or_Gactin"_eigenrmsf_landscape_calpha.xvg -first 1 -last 2 -nframes 10

  gmx sham -f "\$mutant"_"\$Factin_or_Gactin"_2d_landscape_calpha.xvg -ls "\$mutant"_"\$Factin_or_Gactin"_gibbs_landscape_calpha.xpm -notime -gmax 10
  gmx xpm2ps -f "\$mutant"_"\$Factin_or_Gactin"_gibbs_landscape_calpha.xpm -o "\$mutant"_"\$Factin_or_Gactin"_landscape_plot_calpha.eps -rainbow red

  #along own trajectory
  echo -e 'C-alpha\n C-alpha' | gmx anaeig -v "\$mutant"_"\$Factin_or_Gactin"_eigenval_landscape_calpha.trr -eig "\$mutant"_"\$Factin_or_Gactin"_eigenval_landscape_calpha.xvg -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -extr "\$mutant"_"\$Factin_or_Gactin"_extreme_calpha_OWN_VECTORS.pdb -2d "\$mutant"_"\$Factin_or_Gactin"_2d_landscape_calpha_OWN_VECTORS.xvg -rmsf "\$mutant"_"\$Factin_or_Gactin"_eigenrmsf_landscape_calpha_OWN_VECTORS.xvg -first 1 -last 2 -nframes 10

  gmx sham -f "\$mutant"_"\$Factin_or_Gactin"_2d_landscape_calpha_OWN_VECTORS.xvg -ls "\$mutant"_"\$Factin_or_Gactin"_gibbs_landscape_calpha_OWN_VECTORS.xpm -notime -gmax 10
  gmx xpm2ps -f "\$mutant"_"\$Factin_or_Gactin"_gibbs_landscape_calpha_OWN_VECTORS.xpm -o "\$mutant"_"\$Factin_or_Gactin"_landscape_plot_calpha_OWN_VECTORS.eps -rainbow red

  #for EQ time points
  echo -e 'C-alpha\n C-alpha' | gmx covar -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -v "\$mutant"_"\$Factin_or_Gactin"_eigenval_landscape_calpha_EQ.trr -o "\$mutant"_"\$Factin_or_Gactin"_eigenval_landscape_calpha_EQ.xvg -b \$((eq_start_time*1000)) -e \$((simulation_length*1000))
  echo -e 'C-alpha\n C-alpha' | gmx anaeig -v "\$directory_path"/WT/COMBINED/WT_"\$Factin_or_Gactin"_eigenval_landscape_calpha_EQ.trr -eig "\$directory_path"/WT/COMBINED/WT_"\$Factin_or_Gactin"_eigenval_landscape_calpha_EQ.xvg -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -extr "\$mutant"_"\$Factin_or_Gactin"_extreme_calpha_EQ.pdb -2d "\$mutant"_"\$Factin_or_Gactin"_2d_landscape_calpha_EQ.xvg -rmsf "\$mutant"_"\$Factin_or_Gactin"_eigenrmsf_landscape_calpha_EQ.xvg -first 1 -last 2 -nframes 10 -b \$((eq_start_time*1000)) -e \$((simulation_length*1000))

  gmx sham -f "\$mutant"_"\$Factin_or_Gactin"_2d_landscape_calpha_EQ.xvg -ls "\$mutant"_"\$Factin_or_Gactin"_gibbs_landscape_calpha_EQ.xpm -notime -gmax 10
  gmx xpm2ps -f "\$mutant"_"\$Factin_or_Gactin"_gibbs_landscape_calpha_EQ.xpm -o "\$mutant"_"\$Factin_or_Gactin"_landscape_plot_calpha_EQ.eps -rainbow red

  echo -e 'C-alpha\n C-alpha' | gmx anaeig -v "\$mutant"_"\$Factin_or_Gactin"_eigenval_landscape_calpha_EQ.trr -eig "\$mutant"_"\$Factin_or_Gactin"_eigenval_landscape_calpha_EQ.xvg -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -extr "\$mutant"_"\$Factin_or_Gactin"_extreme_calpha_EQ_OWN_VECTORS.pdb -2d "\$mutant"_"\$Factin_or_Gactin"_2d_landscape_calpha_EQ_OWN_VECTORS.xvg -rmsf "\$mutant"_"\$Factin_or_Gactin"_eigenrmsf_landscape_calpha_EQ_OWN_VECTORS.xvg -first 1 -last 2 -nframes 10 -b \$((eq_start_time*1000)) -e \$((simulation_length*1000))

  gmx sham -f "\$mutant"_"\$Factin_or_Gactin"_2d_landscape_calpha_EQ_OWN_VECTORS.xvg -ls "\$mutant"_"\$Factin_or_Gactin"_gibbs_landscape_calpha_EQ_OWN_VECTORS.xpm -notime -gmax 10
  gmx xpm2ps -f "\$mutant"_"\$Factin_or_Gactin"_gibbs_landscape_calpha_EQ_OWN_VECTORS.xpm -o "\$mutant"_"\$Factin_or_Gactin"_landscape_plot_calpha_EQ_OWN_VECTORS.eps -rainbow red

fi

  #PCA for target chain only along WT then mutant trajectory

if [ \$mutant == 'WT' ]; then
  echo -e 'C-alpha_&_Chain_'\$target_chain'_protein\n C-alpha_&_Chain_'\$target_chain'_protein' | gmx covar -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -v "\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_eigenval_landscape_calpha.trr -o "\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_eigenval_landscape_calpha.xvg
  echo -e 'C-alpha_&_Chain_'\$target_chain'_protein\n C-alpha_&_Chain_'\$target_chain'_protein' | gmx anaeig -v "\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_eigenval_landscape_calpha.trr -eig "\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_eigenval_landscape_calpha.xvg -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -extr "\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_extreme_calpha.pdb -2d "\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_2d_landscape_calpha.xvg -rmsf "\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_eigenrmsf_landscape_calpha.xvg -first 1 -last 2 -nframes 10

  gmx sham -f "\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_2d_landscape_calpha.xvg -ls "\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_gibbs_landscape_calpha.xpm -notime -gmax 10
  gmx xpm2ps -f "\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_gibbs_landscape_calpha.xpm -o "\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_landscape_plot_calpha.eps -rainbow red

  #for EQ timepoints
  echo -e 'C-alpha_&_Chain_'\$target_chain'_protein\n C-alpha_&_Chain_'\$target_chain'_protein' | gmx covar -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -v "\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_eigenval_landscape_calpha_EQ.trr -o "\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_eigenval_landscape_calpha_EQ.xvg -b \$((eq_start_time*1000)) -e \$((simulation_length*1000))
  echo -e 'C-alpha_&_Chain_'\$target_chain'_protein\n C-alpha_&_Chain_'\$target_chain'_protein' | gmx anaeig -v "\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_eigenval_landscape_calpha_EQ.trr -eig "\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_eigenval_landscape_calpha_EQ.xvg -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -extr "\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_extreme_calpha_EQ.pdb -2d "\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_2d_landscape_calpha_EQ.xvg -rmsf "\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_eigenrmsf_landscape_calpha_EQ.xvg -first 1 -last 2 -nframes 10 -b \$((eq_start_time*1000)) -e \$((simulation_length*1000))

  gmx sham -f "\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_2d_landscape_calpha_EQ.xvg -ls "\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_gibbs_landscape_calpha_EQ.xpm -notime -gmax 10
  gmx xpm2ps -f "\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_gibbs_landscape_calpha_EQ.xpm -o "\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_landscape_plot_calpha_EQ.eps -rainbow red


else
  #along WT trajectory
  echo -e 'C-alpha_&_Chain_'\$target_chain'_protein\n C-alpha_&_Chain_'\$target_chain'_protein' | gmx covar -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -v "\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_eigenval_landscape_calpha.trr -o "\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_eigenval_landscape_calpha.xvg
  echo -e 'C-alpha_&_Chain_'\$target_chain'_protein\n C-alpha_&_Chain_'\$target_chain'_protein' | gmx anaeig -v "\$directory_path"/WT/COMBINED/WT_"\$Factin_or_Gactin"_Chain_"\$target_chain"_eigenval_landscape_calpha.trr -eig "\$directory_path"/WT/COMBINED/WT_"\$Factin_or_Gactin"_Chain_"\$target_chain"_eigenval_landscape_calpha.xvg -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -extr "\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_extreme_calpha.pdb -2d "\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_2d_landscape_calpha.xvg -rmsf "\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_eigenrmsf_landscape_calpha.xvg -first 1 -last 2 -nframes 10

  gmx sham -f "\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_2d_landscape_calpha.xvg -ls "\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_gibbs_landscape_calpha.xpm -notime -gmax 10
  gmx xpm2ps -f "\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_gibbs_landscape_calpha.xpm -o "\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_landscape_plot_calpha.eps -rainbow red

  #along own trajectory
  echo -e 'C-alpha_&_Chain_'\$target_chain'_protein\n C-alpha_&_Chain_'\$target_chain'_protein' | gmx anaeig -v "\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_eigenval_landscape_calpha.trr -eig "\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_eigenval_landscape_calpha.xvg -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -extr "\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_extreme_calpha_OWN_VECTORS.pdb -2d "\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_2d_landscape_calpha_OWN_VECTORS.xvg -rmsf "\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_eigenrmsf_landscape_calpha_OWN_VECTORS.xvg -first 1 -last 2 -nframes 10

  gmx sham -f "\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_2d_landscape_calpha_OWN_VECTORS.xvg -ls "\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_gibbs_landscape_calpha_OWN_VECTORS.xpm -notime -gmax 10
  gmx xpm2ps -f "\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_gibbs_landscape_calpha_OWN_VECTORS.xpm -o "\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_landscape_plot_calpha_OWN_VECTORS.eps -rainbow red

  #for EQ timepoints
  echo -e 'C-alpha_&_Chain_'\$target_chain'_protein\n C-alpha_&_Chain_'\$target_chain'_protein' | gmx covar -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -v "\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_eigenval_landscape_calpha_EQ.trr -o "\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_eigenval_landscape_calpha_EQ.xvg -b \$((eq_start_time*1000)) -e \$((simulation_length*1000))
  echo -e 'C-alpha_&_Chain_'\$target_chain'_protein\n C-alpha_&_Chain_'\$target_chain'_protein' | gmx anaeig -v "\$directory_path"/WT/COMBINED/WT_"\$Factin_or_Gactin"_Chain_"\$target_chain"_eigenval_landscape_calpha_EQ.trr -eig "\$directory_path"/WT/COMBINED/WT_"\$Factin_or_Gactin"_Chain_"\$target_chain"_eigenval_landscape_calpha_EQ.xvg -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -extr "\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_extreme_calpha_EQ.pdb -2d "\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_2d_landscape_calpha_EQ.xvg -rmsf "\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_eigenrmsf_landscape_calpha_EQ.xvg -first 1 -last 2 -nframes 10 -b \$((eq_start_time*1000)) -e \$((simulation_length*1000))

  gmx sham -f "\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_2d_landscape_calpha_EQ.xvg -ls "\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_gibbs_landscape_calpha_EQ.xpm -notime -gmax 10
  gmx xpm2ps -f "\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_gibbs_landscape_calpha_EQ.xpm -o "\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_landscape_plot_calpha_EQ.eps -rainbow red

  echo -e 'C-alpha_&_Chain_'\$target_chain'_protein\n C-alpha_&_Chain_'\$target_chain'_protein' | gmx anaeig -v "\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_eigenval_landscape_calpha_EQ.trr -eig "\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_eigenval_landscape_calpha_EQ.xvg -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -extr "\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_extreme_calpha_EQ_OWN_VECTORS.pdb -2d "\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_2d_landscape_calpha_EQ_OWN_VECTORS.xvg -rmsf "\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_eigenrmsf_landscape_calpha_EQ_OWN_VECTORS.xvg -first 1 -last 2 -nframes 10 -b \$((eq_start_time*1000)) -e \$((simulation_length*1000))

  gmx sham -f "\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_2d_landscape_calpha_EQ_OWN_VECTORS.xvg -ls "\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_gibbs_landscape_calpha_EQ_OWN_VECTORS.xpm -notime -gmax 10
  gmx xpm2ps -f "\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_gibbs_landscape_calpha_EQ_OWN_VECTORS.xpm -o "\$mutant"_"\$Factin_or_Gactin"_Chain_"\$target_chain"_landscape_plot_calpha_EQ_OWN_VECTORS.eps -rainbow red

fi

#ramachandran, run for top 3 clusters (most common conformations)

gmx rama -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -o "\$mutant"_"\$Factin_or_Gactin"_ramachandran1.xvg -b \$combined_cluster1 -e \$combined_cluster1
gmx rama -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -o "\$mutant"_"\$Factin_or_Gactin"_ramachandran2.xvg -b \$combined_cluster2 -e \$combined_cluster2
gmx rama -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -o "\$mutant"_"\$Factin_or_Gactin"_ramachandran3.xvg -b \$combined_cluster3 -e \$combined_cluster3



#gmx confrms to compute rmsd differences between structures
#first calculate differences of top clusters to starting structure

#second calculate differences of top clusters to WT top cluster

if [ \$mutant == 'WT' ]; then
  echo -e 'Backbone_&_Chain_'\$target_chain'_protein\n Backbone_&_Chain_'\$target_chain'_protein' | gmx trjconv -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -o "\$mutant"_backbone_combined_cluster1.pdb -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -fit rot+trans -tu ns -dump "\$combined_cluster1"
  echo -e 'Backbone_&_Chain_'\$target_chain'_protein\n Backbone_&_Chain_'\$target_chain'_protein' | gmx trjconv -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -o "\$mutant"_backbone_combined_cluster2.pdb -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -fit rot+trans -tu ns -dump "\$combined_cluster2"
  echo -e 'Backbone_&_Chain_'\$target_chain'_protein\n Backbone_&_Chain_'\$target_chain'_protein' | gmx trjconv -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -o "\$mutant"_backbone_combined_cluster3.pdb -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -fit rot+trans -tu ns -dump "\$combined_cluster3"

  echo -e 'Backbone_&_Chain_'\$target_chain'_protein\n Backbone_&_Chain_'\$target_chain'_protein' | gmx trjconv -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -o "\$mutant"_backbone_combined_firstframe.pdb -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -fit rot+trans -tu ns -dump 0

  echo -e 'Protein\n Protein' | gmx confrms -f1 "\$mutant"_backbone_combined_cluster1.pdb -f2 "\$mutant"_backbone_combined_firstframe.pdb -o cluster1_"\$mutant"_firstframe_confrms.pdb -label -bfac
  echo -e 'Protein\n Protein' | gmx confrms -f1 "\$mutant"_backbone_combined_cluster2.pdb -f2 "\$mutant"_backbone_combined_firstframe.pdb -o cluster2_"\$mutant"_firstframe_confrms.pdb -label -bfac
  echo -e 'Protein\n Protein' | gmx confrms -f1 "\$mutant"_backbone_combined_cluster3.pdb -f2 "\$mutant"_backbone_combined_firstframe.pdb -o cluster3_"\$mutant"_firstframe_confrms.pdb -label -bfac

else
  echo -e 'Backbone_&_Chain_'\$target_chain'_protein\n Backbone_&_Chain_'\$target_chain'_protein' | gmx trjconv -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -o "\$mutant"_backbone_combined_cluster1.pdb -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -fit rot+trans -tu ns -dump "\$combined_cluster1"
  echo -e 'Protein\n Protein' | gmx confrms -f1 "\$directory_path"/WT/COMBINED/WT_backbone_combined_cluster1.pdb -f2 "\$mutant"_backbone_combined_cluster1.pdb -o cluster1_"\$mutant"_WT_confrms.pdb -label -bfac
  echo -e 'Backbone_&_Chain_'\$target_chain'_protein\n Backbone_&_Chain_'\$target_chain'_protein' | gmx trjconv -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -o "\$mutant"_backbone_combined_cluster2.pdb -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -fit rot+trans -tu ns -dump "\$combined_cluster2"
  echo -e 'Protein\n Protein' | gmx confrms -f1 "\$directory_path"/WT/COMBINED/WT_backbone_combined_cluster2.pdb -f2 "\$mutant"_backbone_combined_cluster2.pdb -o cluster2_"\$mutant"_WT_confrms.pdb -label -bfac
  echo -e 'Backbone_&_Chain_'\$target_chain'_protein\n Backbone_&_Chain_'\$target_chain'_protein' | gmx trjconv -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -o "\$mutant"_backbone_combined_cluster3.pdb -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -fit rot+trans -tu ns -dump "\$combined_cluster3"
  echo -e 'Protein\n Protein' | gmx confrms -f1 "\$directory_path"/WT/COMBINED/WT_backbone_combined_cluster3.pdb -f2 "\$mutant"_backbone_combined_cluster3.pdb -o cluster3_"\$mutant"_WT_confrms.pdb -label -bfac

  echo -e 'Backbone_&_Chain_'\$target_chain'_protein\n Backbone_&_Chain_'\$target_chain'_protein' | gmx trjconv -s ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_md.tpr -f "\$mutant"_"\$Factin_or_Gactin"_protein_COMBINED.xtc -o "\$mutant"_backbone_combined_firstframe.pdb -n ../Replicate1/"\$mutant"_"\$Factin_or_Gactin"_analysis_index.ndx -fit rot+trans -tu ns -dump 0

  echo -e 'Protein\n Protein' | gmx confrms -f1 "\$mutant"_backbone_combined_cluster1.pdb -f2 "\$mutant"_backbone_combined_firstframe.pdb -o cluster1_"\$mutant"_firstframe_confrms.pdb -label -bfac
  echo -e 'Protein\n Protein' | gmx confrms -f1 "\$mutant"_backbone_combined_cluster2.pdb -f2 "\$mutant"_backbone_combined_firstframe.pdb -o cluster2_"\$mutant"_firstframe_confrms.pdb -label -bfac
  echo -e 'Protein\n Protein' | gmx confrms -f1 "\$mutant"_backbone_combined_cluster3.pdb -f2 "\$mutant"_backbone_combined_firstframe.pdb -o cluster3_"\$mutant"_firstframe_confrms.pdb -label -bfac
fi



#DYNAMIC NETWORK ANALYSIS
cd "\$directory_path"/"\$mutant"/COMBINED/

mkdir "\$directory_path"/"\$mutant"/COMBINED/NETWORK_ANALYSIS

cd "\$directory_path"/"\$mutant"/COMBINED/NETWORK_ANALYSIS

#copy files from network tutorial
cp "\$network_tutorial_path"/subopt "\$directory_path"/"\$mutant"/COMBINED/NETWORK_ANALYSIS
cp "\$network_tutorial_path"/gncommunities "\$directory_path"/"\$mutant"/COMBINED/NETWORK_ANALYSIS

#for network analysis, we need to specify that we want to use our own VMD install to run
#then we need to specify the paths for carma and catdcd which are modules used to generate networks
module load vmd
export PATH=\$PATH:/home/ksteffen/NETWORK_ANALYSIS/network-tutorial/carma-2.01/bin/linux
export PATH=\$PATH:/home/ksteffen/NETWORK_ANALYSIS/network-tutorial/catdcd-4.0b/LINUXAMD64/bin/catdcd4.0

cp \$directory_path/\$mutant/COMBINED/\${mutant}_\${Factin_or_Gactin}_\${target_chain}_protein_COMBINED.pdb \$directory_path/\$mutant/COMBINED/NETWORK_ANALYSIS
cp \$directory_path/\$mutant/COMBINED/\${mutant}_\${Factin_or_Gactin}_\${target_chain}_protein_COMBINED.xtc \$directory_path/\$mutant/COMBINED/NETWORK_ANALYSIS


#make scrip tthat will generat psf file and dcd file
make_psf="generate_psf.tcl"

cat <<EOL > "\$make_psf"

mol new \${mutant}_\${Factin_or_Gactin}_\${target_chain}_protein_COMBINED.pdb
mol addfile \${mutant}_\${Factin_or_Gactin}_\${target_chain}_protein_COMBINED.xtc type xtc waitfor all
animate write dcd \${mutant}_\${Factin_or_Gactin}_\${target_chain}_protein_COMBINED.dcd waitfor all
package require psfgen
resetpsf
topology /home/ksteffen/vmd_library/plugins/noarch/tcl/readcharmmtop1.2/top_all36_prot.rtf
pdbalias residue HIS HSD
pdbalias atom ILE CD1 CD
segment XP1 { pdb \${mutant}_\${Factin_or_Gactin}_\${target_chain}_protein_COMBINED.pdb }
patch CTER XP1:\$number_of_residues
patch NTER XP1:1
coordpdb \${mutant}_\${Factin_or_Gactin}_\${target_chain}_protein_COMBINED.pdb XP1
guesscoord
writepdb \${mutant}_\${target_chain}.pdb
writepsf \${mutant}_\${target_chain}.psf
EOL



#export variables to script that will generate psf for network
sed -i "1s/^/set target_chain \$target_chain\n/" \$directory_path/\$mutant/COMBINED/NETWORK_ANALYSIS/generate_psf.tcl
sed -i "1s/^/set Factin_or_Gactin \$Factin_or_Gactin\n/" \$directory_path/\$mutant/COMBINED/NETWORK_ANALYSIS/generate_psf.tcl
sed -i "1s/^/set mutant \$mutant\n/" \$directory_path/\$mutant/COMBINED/NETWORK_ANALYSIS/generate_psf.tcl

#run script to generate psf and dcd files
vmd -dispdev text -e \$directory_path/\$mutant/COMBINED/NETWORK_ANALYSIS/generate_psf.tcl

# Define the filename
network_file="network.config"

# Create the network.config file with the details for network analysis plugin to run
#the psf, dcd, systemselection, and node selection are REQUIRED
#clusters associated  specified set of atoms within a node. Default is atoms within a given residue are clustered together
#restrictions are constraints that specify edges to be removed from final network. Notsameresidue disallows edges between nodes in the same residue
cat <<EOL > "\$network_file"

>Psf
\${mutant}_\${target_chain}.psf

>Dcds
\${mutant}_\${Factin_or_Gactin}_\${target_chain}_protein_COMBINED.dcd

>SystemSelection
(chain X) and (not hydrogen)

>NodeSelection
(name CA)

>Clusters

>Restrictions
notSameResidue
notNeighboringCAlpha
EOL

#make script to run network analysis
tcl_file="network_analysis.tcl"

cat <<'EOL' > "\$tcl_file"

set path_to_network "/home/ksteffen/vmd-1.9.4a55/plugins/noarch/tcl/networkview1.41/"

source "\$path_to_network/adjacencyMatrix.tcl"
source "\$path_to_network/communityNetwork.tcl"
source "\$path_to_network/networkSetup.tcl"
source "\$path_to_network/networkView.tcl"
source "\$path_to_network/pkgIndex.tcl"
source "\$path_to_network/suboptimalPaths.tcl"

mol new \${mutant}_\${target_chain}.psf

mol addfile \${mutant}_\${Factin_or_Gactin}_\${target_chain}_protein_COMBINED.dcd type dcd waitfor all

networkSetup network.config

./gncommunities contact.dat communities.out

::NetworkView::networkView 0 contact.dat

set F375 [::NetworkView::getNodesFromSelection "chain X and resid 375 and name CA"]
set C374 [::NetworkView::getNodesFromSelection "chain X and resid 374 and name CA"]
set H371 [::NetworkView::getNodesFromSelection "chain X and resid 371 and name CA"]

set D1 [::NetworkView::getNodesFromSelection "chain X and resid 1 and name CA"]
set D2 [::NetworkView::getNodesFromSelection "chain X and resid 2 and name CA"]
set E3 [::NetworkView::getNodesFromSelection "chain X and resid 3 and name CA"]
set E4 [::NetworkView::getNodesFromSelection "chain X and resid 4 and name CA"]
set T5 [::NetworkView::getNodesFromSelection "chain X and resid 5 and name CA"]

set Q137 [::NetworkView::getNodesFromSelection "chain X and resid 137 and name CA"]
set S14 [::NetworkView::getNodesFromSelection "chain X and resid 14 and name CA"]
set G74 [::NetworkView::getNodesFromSelection "chain X and resid 74 and name CA"]

set R312 [::NetworkView::getNodesFromSelection "chain X and resid 312 and name CA"]
set H88 [::NetworkView::getNodesFromSelection "chain X and resid 88 and name CA"]

set T303 [::NetworkView::getNodesFromSelection "chain X and resid 303 and name CA"]

set I34 [::NetworkView::getNodesFromSelection "chain X and resid 34 and name CA"]
set V54 [::NetworkView::getNodesFromSelection "chain X and resid 54 and name CA"]
set V35 [::NetworkView::getNodesFromSelection "chain X and resid 35 and name CA"]
set V43 [::NetworkView::getNodesFromSelection "chain X and resid 43 and name CA"]

./subopt contact.dat F375-D1_subopt 20 \$F375 \$D1
./subopt contact.dat F375-D2_subopt 20 \$F375 \$D2
./subopt contact.dat F375-E3_subopt 20 \$F375 \$E3
./subopt contact.dat F375-E4_subopt 20 \$F375 \$E4
./subopt contact.dat F375-T5_subopt 20 \$F375 \$T5
./subopt contact.dat F375-Q137_subopt 20 \$F375 \$Q137
./subopt contact.dat F375-S14_subopt 20 \$F375 \$S14
./subopt contact.dat F375-G74_subopt 20 \$F375 \$G74
./subopt contact.dat F375-R312_subopt 20 \$F375 \$R312
./subopt contact.dat F375-H88_subopt 20 \$F375 \$H88
./subopt contact.dat F375-T303_subopt 20 \$F375 \$T303
./subopt contact.dat F375-I34_subopt 20 \$F375 \$I34
./subopt contact.dat F375-V54_subopt 20 \$F375 \$V54
./subopt contact.dat F375-V35_subopt 20 \$F375 \$V35
./subopt contact.dat F375-V43_subopt 20 \$F375 \$V43

./subopt contact.dat C374-D1_subopt 20 \$C374 \$D1
./subopt contact.dat C374-D2_subopt 20 \$C374 \$D2
./subopt contact.dat C374-E3_subopt 20 \$C374 \$E3
./subopt contact.dat C374-E4_subopt 20 \$C374 \$E4
./subopt contact.dat C374-T5_subopt 20 \$C374 \$T5
./subopt contact.dat C374-Q137_subopt 20 \$C374 \$Q137
./subopt contact.dat C374-S14_subopt 20 \$C374 \$S14
./subopt contact.dat C374-G74_subopt 20 \$C374 \$G74
./subopt contact.dat C374-R312_subopt 20 \$C374 \$R312
./subopt contact.dat C374-H88_subopt 20 \$C374 \$H88
./subopt contact.dat C374-T303_subopt 20 \$C374 \$T303
./subopt contact.dat C374-I34_subopt 20 \$C374 \$I34
./subopt contact.dat C374-V54_subopt 20 \$C374 \$V54
./subopt contact.dat C374-V35_subopt 20 \$C374 \$V35
./subopt contact.dat C374-V43_subopt 20 \$C374 \$V43

./subopt contact.dat H371-D1_subopt 20 \$H371 \$D1
./subopt contact.dat H371-D2_subopt 20 \$H371 \$D2
./subopt contact.dat H371-E3_subopt 20 \$H371 \$E3
./subopt contact.dat H371-E4_subopt 20 \$H371 \$E4
./subopt contact.dat H371-T5_subopt 20 \$H371 \$T5
./subopt contact.dat H371-Q137_subopt 20 \$H371 \$Q137
./subopt contact.dat H371-S14_subopt 20 \$H371 \$S14
./subopt contact.dat H371-G74_subopt 20 \$H371 \$G74
./subopt contact.dat H371-R312_subopt 20 \$H371 \$R312
./subopt contact.dat H371-H88_subopt 20 \$H371 \$H88
./subopt contact.dat H371-T303_subopt 20 \$H371 \$T303
./subopt contact.dat H371-I34_subopt 20 \$H371 \$I34
./subopt contact.dat H371-V54_subopt 20 \$H371 \$V54
./subopt contact.dat H371-V35_subopt 20 \$H371 \$V35
./subopt contact.dat H371-V43_subopt 20 \$H371 \$V43

EOL

# Add variables using sed
sed -i "1s/^/set target_chain \$target_chain\n/" "\$tcl_file"
sed -i "1s/^/set Factin_or_Gactin \$Factin_or_Gactin\n/" "\$tcl_file"
sed -i "1s/^/set mutant \$mutant\n/" "\$tcl_file"


#the lines defining the variables need to be added with sed so that they can accurately call the variables from earlier in this script
sed -i "1s/^/set target_chain \$target_chain\n/" \$directory_path/\$mutant/COMBINED/NETWORK_ANALYSIS/network_analysis.tcl
sed -i "1s/^/set Factin_or_Gactin \$Factin_or_Gactin\n/" \$directory_path/\$mutant/COMBINED/NETWORK_ANALYSIS/network_analysis.tcl
sed -i "1s/^/set mutant \$mutant\n/" \$directory_path/\$mutant/COMBINED/NETWORK_ANALYSIS/network_analysis.tcl
sed -i "1s/^/set network_view \$path_to_vmd\n/" \$directory_path/\$mutant/COMBINED/NETWORK_ANALYSIS/network_analysis.tcl

#run script to generate network analysis
vmd -dispdev text -e \$directory_path/\$mutant/COMBINED/NETWORK_ANALYSIS/network_analysis.tcl


EOF
)

job_script="${mutant}_job_script.sh"
echo "$script_content" > "$job_script"

