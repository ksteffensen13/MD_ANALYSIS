from pymol import cmd
from pymol import math
from pymol import stored
import numpy as np
import sys

# Extract command-line arguments. Here, we can set the start/end distance/other values for the docking function every time we call this script
#the line below will just call variables from within the UNIX script running this pythons cript
if len(sys.argv) != 8:
    print("Usage: python prepare_for_Tn_MD.py <mutant> <initial_file_directory> <path_for_output> <path_for_pymol_axes_script> <thin_filament_state> <initial_distance> <docking_distance>")
    sys.exit(1)

mutant = sys.argv[1]
initial_file_directory = sys.argv[2]
path_for_output = sys.argv[3]
path_for_pymol_axes_script = sys.argv[4]
state = sys.argv[5]
initial_distance = sys.argv[6]
docking_distance = sys.argv[7]

#first, add TPM chains to minimized actin filament
for i in range(1, 6):
    cmd.delete('all')
    cmd.load(f'{initial_file_directory}DOCKING_{initial_distance}A_{mutant}{i}_Tn_{state}_MINIMIZED_SD.pdb')
    # Pymol changes chain lettering, so take post-minimization structure and fix it
    cmd.alter(selection='chain M', expression="chain='U'")
    cmd.alter(selection='chain L', expression="chain='V'")
    cmd.alter(selection='chain K', expression="chain='Y'")
    cmd.alter(selection='resn ADP', expression="chain='R'")
    cmd.alter(selection='resn PO4', expression="chain='T'")
    cmd.alter(selection='resn MG', expression="chain='S'")
    #to make all TPM translations consistent across structures, we'll move the structure to the same point in Pymol
    #That point will be the origin
    #First define the origin, calculate relative position from camera view, then shift everything to the origin
    cmd.origin(position=[0, 0, 0])
    cmd.run(f'{path_for_pymol_axes_script}')
    cmd.orient(selection='all')
    # get the coordinates of the structure upon loading
    v = cmd.get_view()
    # first translate to origin. In this case, we set x,y,z to a formula
    # the formula says to subtract the x,y, and z values from the translation matrix from the x,y,z coordinates
    # aka, what are my current x/y/z in the matrix v (corresponding to matrix indices 12,13,14)? Subtract those from current x,y,z coordinates
    # the result is 0,0, and 0 because we want to move to the origin
    cmd.alter_state(1, f'DOCKING_{initial_distance}A_{mutant}{i}_Tn_{state}_MINIMIZED_SD', f"x,y,z=x-{v[12]},y-{v[13]},z-{v[14]}")
    # then we have to use the translation matrix again to apply a rotation as per https://pymolwiki.org/index.php/Transform_selection
    cmd.alter_state(1, f'DOCKING_{initial_distance}A_{mutant}{i}_Tn_{state}_MINIMIZED_SD',
                    f"x,y,z=x*{v[0]}+y*{v[3]}+z*{v[6]},x*{v[1]}+y*{v[4]}+z*{v[7]},x*{v[2]}+y*{v[5]}+z*{v[8]}")
    cmd.orient(selection='all')
    cmd.translate(vector=[0.0, 0.5, 0.0], selection=f'chain U+V+Y', camera=0)
    cmd.save(f'{path_for_output}DOCKING_{docking_distance}A_{mutant}{i}_Tn_{state}.pdb')
    cmd.delete('all')
