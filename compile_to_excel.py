import mutant_functions_for_excel as mut

##################################################################################################################
##THIS SECTION REQUIRES MANUAL ENTRY

number_of_subunits = 1
target_chain = 'A'
mutant = 'R312H'
F_or_G = 'Gactin'
#not nucleotide state, just ADP or ATP
nucleotide = 'ATP'
nucleotide_state = 'ATP'
#length of simulation in NANOSECONDS
simulation_length = 200
##NOTE: FOR path_for_compiled_data, important to have the final / at the end
path_for_compiled_data = f'/Users/karl/Desktop/Computational_Work/Mutant_Gactin/ANALYSES/{mutant}/'

#SET FILE PATHS. NOTE: important to have the final / at the end
rmsd_rep1_path = f'{path_for_compiled_data}Replicate1/RMSD/'
rmsd_rep2_path = f'{path_for_compiled_data}Replicate2/RMSD/'
rmsd_rep3_path = f'{path_for_compiled_data}Replicate3/RMSD/'

rmsf_rep1_path = f'{path_for_compiled_data}Replicate1/RMSF/'
rmsf_rep2_path = f'{path_for_compiled_data}Replicate2/RMSF/'
rmsf_rep3_path = f'{path_for_compiled_data}Replicate3/RMSF/'

rmsdev_rep1_path = f'{path_for_compiled_data}Replicate1/RMSDEV/'
rmsdev_rep2_path = f'{path_for_compiled_data}Replicate2/RMSDEV/'
rmsdev_rep3_path = f'{path_for_compiled_data}Replicate3/RMSDEV/'

PCA_path = f'{path_for_compiled_data}COMBINED/PCA/OWN_VECTORS/'

hbonds_rep1_path = f'{path_for_compiled_data}Replicate1/HBONDS/'
hbonds_rep2_path = f'{path_for_compiled_data}Replicate2/HBONDS/'
hbonds_rep3_path = f'{path_for_compiled_data}Replicate3/HBONDS/'

allostery_distance_folder_path_rep1 = f'{path_for_compiled_data}Replicate1/DISTANCES/ALLOSTERY'
allostery_distance_folder_path_rep2 = f'{path_for_compiled_data}Replicate2/DISTANCES/ALLOSTERY'
allostery_distance_folder_path_rep3 = f'{path_for_compiled_data}Replicate3/DISTANCES/ALLOSTERY'

interprotomer_distance_folder_path_rep1 = f'{path_for_compiled_data}Replicate1/DISTANCES/INTERPROTOMER'
interprotomer_distance_folder_path_rep2 = f'{path_for_compiled_data}Replicate2/DISTANCES/INTERPROTOMER'
interprotomer_distance_folder_path_rep3 = f'{path_for_compiled_data}Replicate3/DISTANCES/INTERPROTOMER'

interprotomer2_folder_path_rep1 = f'{path_for_compiled_data}Replicate1/DISTANCES/INTERPROTOMER2'
interprotomer2_folder_path_rep2 = f'{path_for_compiled_data}Replicate2/DISTANCES/INTERPROTOMER2'
interprotomer2_folder_path_rep3 = f'{path_for_compiled_data}Replicate3/DISTANCES/INTERPROTOMER2'

interstrand_distance_folder_path_rep1 = f'{path_for_compiled_data}Replicate1/DISTANCES/INTERSTRAND'
interstrand_distance_folder_path_rep2 = f'{path_for_compiled_data}Replicate2/DISTANCES/INTERSTRAND'
interstrand_distance_folder_path_rep3 = f'{path_for_compiled_data}Replicate3/DISTANCES/INTERSTRAND'

nucleotide_distance_folder_path_rep1 = f'{path_for_compiled_data}Replicate1/DISTANCES/NUCLEOTIDE'
nucleotide_distance_folder_path_rep2 = f'{path_for_compiled_data}Replicate2/DISTANCES/NUCLEOTIDE'
nucleotide_distance_folder_path_rep3 = f'{path_for_compiled_data}Replicate3/DISTANCES/NUCLEOTIDE'

###################################################################################################################
##THIS SECTION RUNS THE SCRIPT

#for calculations

number_of_timepoints = simulation_length*100
number_of_residues = number_of_subunits*375

if target_chain == 'A':
    target_chain_number = 0
if target_chain == 'B':
    target_chain_number = 1
if target_chain == 'C':
    target_chain_number = 2
if target_chain == 'D':
    target_chain_number = 3
if target_chain == 'E':
    target_chain_number = 4
if target_chain == 'F':
    target_chain_number = 5
if target_chain == 'G':
    target_chain_number = 6
if target_chain == 'H':
    target_chain_number = 7
if target_chain == 'I':
    target_chain_number = 8

chains = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L']

#RMSD
#run rmsd processing. File_counter is for helping to sort/organize files numerically
file_counter = 1
mut.process_files_rmsd(rmsd_rep1_path, file_counter, F_or_G, mutant, target_chain)
mut.process_files_rmsd(rmsd_rep2_path, file_counter, F_or_G, mutant, target_chain)
mut.process_files_rmsd(rmsd_rep3_path, file_counter, F_or_G, mutant, target_chain)

mut.make_compiled_excel(F_or_G, mutant, path_for_compiled_data, number_of_timepoints, target_chain)


#RMSF
mut.process_rmsf_files(rmsf_rep1_path, F_or_G, mutant, target_chain)
mut.process_rmsf_files(rmsf_rep2_path, F_or_G, mutant, target_chain)
mut.process_rmsf_files(rmsf_rep3_path, F_or_G, mutant, target_chain)

mut.compile_rmsf(rmsf_rep1_path, rmsf_rep2_path, rmsf_rep3_path, F_or_G, mutant, number_of_residues, path_for_compiled_data, target_chain_number)

#PCA

#mut.process_PCA(PCA_path, F_or_G, mutant)

#mut.compile_PCA(path_for_compiled_data, F_or_G, mutant, PCA_path, number_of_residues)

#HBONDS

mut.process_HBOND_files(hbonds_rep1_path, F_or_G, mutant, chains, target_chain_number, nucleotide, target_chain)
mut.process_HBOND_files(hbonds_rep2_path, F_or_G, mutant, chains, target_chain_number, nucleotide, target_chain)
mut.process_HBOND_files(hbonds_rep3_path, F_or_G, mutant, chains, target_chain_number, nucleotide, target_chain)

mut.compile_hbonds(F_or_G, mutant, path_for_compiled_data, hbonds_rep1_path, hbonds_rep2_path, hbonds_rep3_path,
                   chains, target_chain_number, nucleotide)

#distance_files

mut.process_distance_files(path_for_compiled_data, F_or_G, mutant, number_of_timepoints, allostery_distance_folder_path_rep1,
                           interprotomer_distance_folder_path_rep1, interstrand_distance_folder_path_rep1, nucleotide_distance_folder_path_rep1,
                           chains, target_chain_number, nucleotide, nucleotide_state)

mut.process_distance_files(path_for_compiled_data, F_or_G, mutant, number_of_timepoints, allostery_distance_folder_path_rep2,
                           interprotomer_distance_folder_path_rep2, interstrand_distance_folder_path_rep2, nucleotide_distance_folder_path_rep2,
                           chains, target_chain_number, nucleotide, nucleotide_state)

mut.process_distance_files(path_for_compiled_data, F_or_G, mutant, number_of_timepoints, allostery_distance_folder_path_rep3,
                           interprotomer_distance_folder_path_rep3, interstrand_distance_folder_path_rep3, nucleotide_distance_folder_path_rep3,
                           chains, target_chain_number, nucleotide, nucleotide_state)

mut.compile_distances(F_or_G, mutant, path_for_compiled_data, nucleotide_state,
                      allostery_distance_folder_path_rep1, allostery_distance_folder_path_rep2, allostery_distance_folder_path_rep3,
                      interprotomer_distance_folder_path_rep1, interprotomer_distance_folder_path_rep2, interprotomer_distance_folder_path_rep3,
                      interstrand_distance_folder_path_rep1, interstrand_distance_folder_path_rep2, interstrand_distance_folder_path_rep3,
                      nucleotide_distance_folder_path_rep1, nucleotide_distance_folder_path_rep2, nucleotide_distance_folder_path_rep3)


#RMSF
mut.process_rmsdeviation(rmsdev_rep1_path, F_or_G, mutant, target_chain)
mut.process_rmsdeviation(rmsdev_rep2_path, F_or_G, mutant, target_chain)
mut.process_rmsdeviation(rmsdev_rep3_path, F_or_G, mutant, target_chain)

mut.compile_rmsdeviation(rmsdev_rep1_path, rmsdev_rep2_path, rmsdev_rep3_path, F_or_G, mutant, number_of_residues, path_for_compiled_data, target_chain_number)





#INTERPROTOMER2 DISTANCES

#mut.process_interprotomer_distances2(path_for_compiled_data, F_or_G, mutant, number_of_timepoints, interprotomer2_folder_path_rep1, chains, target_chain_number)
#mut.process_interprotomer_distances2(path_for_compiled_data, F_or_G, mutant, number_of_timepoints, interprotomer2_folder_path_rep2, chains, target_chain_number)
#mut.process_interprotomer_distances2(path_for_compiled_data, F_or_G, mutant, number_of_timepoints, interprotomer2_folder_path_rep3, chains, target_chain_number)

#mut.compile_interprotomer_distances2(F_or_G, mutant, path_for_compiled_data,
                      #interprotomer2_folder_path_rep1, interprotomer2_folder_path_rep2, interprotomer2_folder_path_rep3)
