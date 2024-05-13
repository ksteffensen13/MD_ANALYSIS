from pymol import cmd
from pymol import math
import numpy as np
import math
import sys
import math

def initial_distance(path_for_output, mutant, Factin_or_Gactin, file_to_load, tpm_chains):
    cmd.delete('all')
    #load the initial filament
    cmd.load(f'{file_to_load}')
    # isolate the chain identifiers for tpm chains into separate variables (needed for calculating vectors)
    tpm_chain_list = tpm_chains.split('+')
    tpm_chain1 = tpm_chain_list[0]
    tpm_chain2 = tpm_chain_list[1]
    #calculate a vector through the central actin filament axis
    v1 = draw_axis('chain I', 'chain N')
    # calculate a vector through the center of the 2 tpm chains
    v2 = draw_axis(f'chain {tpm_chain1}', f'chain {tpm_chain2}')
    #calculate the vector from v1 through v2 by calculating the crossproduct
    v3 = np.cross(v1, v2)
    # calculate normalized vector (with a magnitude of 1)
    # need normalized vector so that it scales with our translationd istance
    v3_normalized = v3 / np.linalg.norm(v3)
    #create a transformation matrix for the calculated vector and the desired initial translation distance
    move_vector = [1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0, 0.0 , 0.0, 1.0, 0.0, 5.0 * v3_normalized[0], 5.0 * v3_normalized[1], 5.0 * v3_normalized[2], 1.0]
    # move the TPM chains along the desginated vector/distance
    cmd.transform_selection(f'chain {tpm_chains}', move_vector)
    cmd.select(f'chain {tpm_chains} & i. 254-284')
    cmd.remove('sele')
    cmd.save(f'{path_for_output}DOCKING_5A_{mutant}_{Factin_or_Gactin}0.00A0.00D.pdb')


#NOTE: THESE FUNCTIONS WERE MODIFIED FROM https://pymolwiki.org/index.php/RotationAxis
#  THESE FUNCTIONS ARE USED TO DRAW THE AXIS FROM THE DRAW_AXIS FUNCTION CALL ABOVE

def transf_matrix(chA, chB):
    '''
    DESCRIPTION

    Align two selections/chains, and returns the transformation matrix. I used super to carry out the alignment, likely is possible to use cmd.align and
    is going to be a bit faster, but I think is not going to work well with low-sequence-identity alignments.

    '''
    cmd.create('working', chA)
    cmd.super('working', chB)
    T = cmd.get_object_matrix('working')
    global cmW
    cmW = center_of_Mass('working')
    cmd.delete('working')
    return T


def center_of_Mass(selection):
    '''
    DESCRIPTION

    Calculates the center of mass of a given selection

    '''
    model= cmd.get_model(selection)
    x,y,z=0,0,0
    totmass = 0
    for a in model.atom:
        m = a.get_mass()
        x+= a.coord[0]*m
        y+= a.coord[1]*m
        z+= a.coord[2]*m
        totmass += m
    cM = np.array([x/totmass, y/totmass, z/totmass])
    return cM

def direction_cosines(chA, chB):
    '''
    DESCRIPTION

    Calculates the direction cosines of the rotation axis from the transformation matrix.

    '''
    t=transf_matrix(chA, chB)
    a1= (t[6]-t[9])/math.sqrt((t[6]-t[9])**2+(t[8]-t[2])**2+(t[1]-t[4])**2)
    b1= (t[8]-t[2])/math.sqrt((t[6]-t[9])**2+(t[8]-t[2])**2+(t[1]-t[4])**2)
    c1= (t[1]-t[4])/math.sqrt((t[6]-t[9])**2+(t[8]-t[2])**2+(t[1]-t[4])**2)
    axis = np.array([a1, b1, c1])
    return axis

def angle_axis(chA, chB):
    '''
    DESCRIPTION

    Calculates the rotation angle from the transformation matrix

    '''
    t=transf_matrix(chA, chB)
    angle_rad = math.acos((t[0]+t[5]+t[10]-1)/2)
    return angle_rad

def center_of_rot(chA, chB):
    '''
    DESCRIPTION

    Calculates the center of rotation of the axis

    '''
    cm_Working=center_of_Mass(chA)
    cm_Reference=cmW
    u=direction_cosines(chA, chB)[0]
    u2=u**2
    v=direction_cosines(chA, chB)[1]
    v2=v**2
    w=direction_cosines(chA, chB)[2]
    w2=w**2
    cos_theta=np.cos(angle_axis(chA, chB))
    sin_theta=np.sin(angle_axis(chA, chB))
    fx=cm_Working[0]
    fy=cm_Working[1]
    fz=cm_Working[2]
    x=cm_Reference[0]
    y=cm_Reference[1]
    z=cm_Reference[2]

    T11 = (v2 + w2)*(1-cos_theta)
    T12 = (w*sin_theta)-((u*v)*(1-cos_theta))
    T13 = -(v*sin_theta)-(u*w*(1-cos_theta))
    T14 =fx-((((u2*x)+((u*v)*y)+((u*w)*z))*(1-cos_theta))+(x*cos_theta)+((-(w*y)+(v*z))*sin_theta))
    T21 = -(w*sin_theta)-((u*v)*(1-cos_theta))
    T22 = (u2 + w2)*(1-cos_theta)
    T23 =  (u*sin_theta)-(w*v*(1-cos_theta))
    T24 =fy-((((v*u*x)+(v2*y)+(v*w*z))*(1-cos_theta))+(y*cos_theta)+(((w*x)-(u*z))*sin_theta))
    T31 = (v*sin_theta)-(w*u*(1-cos_theta))
    T32 = -(u*sin_theta)-(w*v*(1-cos_theta))
    T33 = (u2 + v2)*(1-cos_theta)
    T34 =fz-(((((u*x)*w)+((v*y)*w)+(w2*z))*(1-cos_theta))+(z*cos_theta)+((-(v*x)+(u*y))*sin_theta))

    term_lig = np.array([[T11, T12, T13], [T21, T22, T23], [T31, T32, T33]])
    term_ind = np.array([T14, T24, T34])

    sol_lstsq = np.linalg.lstsq(term_lig, term_ind, rcond=None)
    sol = sol_lstsq[0]

    return sol

def nearest_point_to_axis(chA, chB):
    '''
    DESCRIPTION

    Calculates the nearest point of the axis, I use it to create the cgo object.

    '''
    cmA=center_of_Mass(chA)
    cmB=cmW
    cmAver=(cmB+cmA)/2
    vector=np.array([(cmB[0]-cmA[0]), (cmB[1]-cmA[1]), (cmB[2]-cmA[2])])
    moduli_vector=np.linalg.norm(vector)
    vector_director=np.array([(cmB[0]-cmA[0])/moduli_vector, (cmB[1]-cmA[1])/moduli_vector, (cmB[2]-cmA[2])/moduli_vector])
    axis1= direction_cosines(chA, chB)
    sol=center_of_rot(chA, chB)
    term_lig2=np.array([[vector_director[0], vector_director[1], vector_director[2], 0], [1, 0, 0, -axis1[0]], [0, 1, 0, -axis1[1]], [0, 0, 1, -axis1[2]]])
    term_ind2=np.array([(cmAver[0]*(vector_director[0]))+(cmAver[1]*(vector_director[1]))+(cmAver[2]*(vector_director[2])), sol[0], sol[1], sol[2]])
    term_j=(cmAver[0]*vector_director[0])+(cmAver[1]*vector_director[1])+(cmAver[2]*vector_director[2])
    suma_vect_director=vector_director+axis1
    term_ji=(cmAver[0]*suma_vect_director[0])+(cmAver[1]*suma_vect_director[1])+(cmAver[2]*suma_vect_director[2])
    if np.dot(vector_director, axis1) != 0:
        t = ((-np.dot(vector_director, sol))+term_j)/np.dot(vector_director, axis1)
    else:
        t = ((-np.dot(suma_vect_director, sol))+term_ji)/np.dot(suma_vect_director, axis1)
    p = [sol[0]+axis1[0]*t, sol[1]+axis1[1]*t, sol[2]+axis1[2]*t]

    return p

def proyeccion_centroide(selection, chA, chB):
    '''
    DESCRIPTION

    Calculates the proyection of the mass center for the working molecule before being aligned with the reference molecule. For representation purpuses.

    '''
    axis1=np.array([direction_cosines(chA, chB)[0], direction_cosines(chA, chB)[1], direction_cosines(chA, chB)[2]])
    sol=center_of_rot(chA, chB)
    cmSel=center_of_Mass(selection)
    t_cen=np.dot(cmSel, axis1)-np.dot(sol, axis1)
    proy_cen= [sol[0]+(t_cen*axis1[0]), sol[1]+(t_cen*axis1[1]), sol[2]+(t_cen*axis1[2])]
    return proy_cen

def proyeccion_centroide_working(chA, chB):
    '''
    DESCRIPTION

    Calculates the proyection of the mass center for working molecule after being aligned with the reference molecule. For representation purpuses.

    '''
    axis1=np.array([direction_cosines(chA, chB)[0], direction_cosines(chA, chB)[1], direction_cosines(chA, chB)[2]])
    sol=center_of_rot(chA, chB)
    cmSel=cmW
    t_cen=np.dot(cmSel, axis1)-np.dot(sol, axis1)
    proy_cen= [sol[0]+(t_cen*axis1[0]), sol[1]+(t_cen*axis1[1]), sol[2]+(t_cen*axis1[2])]
    return proy_cen

def print_information(T, axis1, angle_degrees,  moduli_vector, obj, x1, y1, z1, x2, y2, z2, w, r1, g1, b1, r2, g2, b2, modu_tr=0):
    '''
    DESCRIPTION

    Print to basic information to the screen.
    '''
    print("#################################################################################################")
    print("Transformation (TTT) matrix")
    print("%8.2f, %8.2f, %8.2f, %8.2f" %(T[0], T[1], T[2], T[3]))
    print("%8.2f, %8.2f, %8.2f, %8.2f" %(T[4], T[5], T[6], T[7]))
    print("%8.2f, %8.2f, %8.2f, %8.2f" %(T[8], T[9], T[10], T[11]))
    print("%8.2f, %8.2f, %8.2f, %8.2f" %(T[12], T[13], T[14], T[15]))
    print(".................................................................................................")
    print("")
    print("The direction cosines of the rotation axis is: %3.2f, %3.2f, %3.2f" %(axis1[0], axis1[1], axis1[2]))
    print("The angle of rotation is %3.2f degrees" %(angle_degrees))
    print("The lenght of the translation vector along the rotation axis is %3.2f Angstroms" %(modu_tr))
    print("The distance between mass centers is %3.2f Angstroms" %(moduli_vector))
    print(".................................................................................................")
    print("")
    print("Lines to be used in a pml script to generate the axis")
    print("")
    print("CYLINDER, %3.2f, %3.2f, %3.2f, %3.2f, %3.2f, %3.2f, %3.2f, %3.2f, %3.2f, %3.2f, %3.2f, %3.2f, %3.2f, 0.0" %(x1, y1, z1, x2, y2, z2, w, r1, g1, b1, r2, g2, b2))
    print("cmd.load_cgo(obj, %3.2f)" %(angle_degrees))
    print("")
    print("#################################################################################################")

def draw_axis(chA, chB):
    T = transf_matrix(chA, chB)
    angle = angle_axis(chA, chB)
    angle_degrees = (angle * 180) / math.pi
    axis1 = [direction_cosines(chA, chB)[0], direction_cosines(chA, chB)[1], direction_cosines(chA, chB)[2]]

    cmA = center_of_Mass(chA)
    cmB = cmW
    cmAver = (cmB + cmA) / 2
    vector = np.array([(cmB[0] - cmA[0]), (cmB[1] - cmA[1]), (cmB[2] - cmA[2])])
    moduli_vector = np.linalg.norm(vector)
    vector_director = np.array([(cmB[0] - cmA[0]) / moduli_vector, (cmB[1] - cmA[1]) / moduli_vector, (cmB[2] - cmA[2]) / moduli_vector])
    pC_A = proyeccion_centroide(chA, chA, chB)
    pC_B = proyeccion_centroide_working(chA, chB)

    trans_vector = np.array([(pC_B[0] - pC_A[0]), (pC_B[1] - pC_A[1]), (pC_B[2] - pC_A[2])])
    modu_tr = np.linalg.norm(trans_vector)
    rota_centroid_rad = np.dot(vector_director, axis1)
    rota_centroid = (rota_centroid_rad * 180) / math.pi
    rota_centroid_absol_0 = np.absolute(rota_centroid)
    rota_centroid_absol = round(rota_centroid_absol_0, 2)

    return axis1

#create a function that creates a transformation matrix
def transmat(vector, dist=1):
    mat = [1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, dist * vector[0], dist * vector[1], dist * vector[2], 1.0]
    return(mat)

def calculate_rotation_matrix(actin_chain1, actin_chain2, angle_degrees):
#run the draw_axis function (above) where we tell it to return the axis coordinates
    #the function will draw an axis down the center of the actin filament between chains I and N
    axis1_values = draw_axis(f'chain {actin_chain1}', f'chain {actin_chain2}')
    #pull out the coordinate of the actin central axis
    x_coord = axis1_values[0]
    y_coord = axis1_values[1]
    z_coord = axis1_values[2]

    # set angle. In this case, 2.5 degrees and convert to radians for calculations
    angle_radians = math.radians(angle_degrees)
    sin_angle = math.sin(angle_radians)
    cos_angle = math.cos(angle_radians)
    #calculate each element of the transformation matrix
    matrix1 = cos_angle + (x_coord * x_coord) * (1 - cos_angle)
    matrix2 = (x_coord * y_coord) * (1 - cos_angle) - (z_coord * sin_angle)
    matrix3 = (y_coord * sin_angle) + (x_coord * z_coord) * (1 - cos_angle)
    matrix4 = 0.0
    matrix5 = (z_coord * sin_angle) + (x_coord * y_coord) * (1 - cos_angle)
    matrix6 = cos_angle + (y_coord * y_coord) * (1 - cos_angle)
    matrix7 = (-x_coord) * sin_angle + (y_coord * z_coord) * (1 - cos_angle)
    matrix8 = 0.0
    matrix9 = (-y_coord) * sin_angle + (x_coord * z_coord) * (1 - cos_angle)
    matrix10 = (x_coord * sin_angle) + (y_coord * z_coord) * (1 - cos_angle)
    matrix11 = cos_angle + (z_coord * z_coord) * (1 - cos_angle)
    matrix12 = 0.0
    matrix13 = 0.0
    matrix14 = 0.0
    matrix15 = 0.0
    matrix16 = 1.0
    rotation_matrix = [matrix1, matrix2, matrix3, matrix4, matrix5, matrix6, matrix7, matrix8, matrix9, matrix10,
                       matrix11, matrix12, matrix13, matrix14, matrix15, matrix16]
    return rotation_matrix


#####################################################END OF FUNCTIONS#############################################





# Extract command-line arguments. Here, we can set the start/end distance/other values for the docking function every time we call this script
#the line below will just call variables from within the UNIX script running this pythons cript
if len(sys.argv) != 14:
    print("Usage: python tpm_translations.py <mutant> <Factin_or_Gactin> <nucleotide_state> <cation> <reference_structure> <path_for_reference_structure> <chains_to_delete> <tpm_chains> <path_for_translation_output> <path_for_pymol_axes_script> <path_to_cluster_chain> <length_of_filament_to_build> <reference_starting_actin_chain>")
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

#index all possible letters for actin chain naming
actin_chains = ('A' 'B' 'C' 'D' 'E' 'F' 'G' 'H' 'I' 'J' 'K' 'L' 'M' 'N' 'O' 'P' 'Q' 'R' 'S' 'T' 'U' 'V' 'W' 'X' 'Y' 'Z'
                'a' 'b' 'c' 'd' 'e' 'f' 'g' 'h' 'i' 'j' 'k' 'l' 'm' 'n' 'o' 'p' 'q' 'r' 's' 't' 'u' 'v' 'w' 'x' 'y' 'z')

#grab the index number of the starting actin chain from the reference structure.
#we need this to calculate the letters of other chains based on where they are relative to starting chain
reference_starting_actin_chain_index = actin_chains.index(reference_starting_actin_chain)

#determine the nucleotide from the state. Need this to select nucleotide for naming/numbering
if nucleotide_state == 'ADP-Pi' or nucleotide_state == 'ADP':
    nucleotide = 'ADP'
else:
    nucleotide = 'ATP'

#first, add TPM chains to minimized actin filament
for i in range(1, 6):
    # load reference structure and remove excess chains (troponin and extra TPM)
    cmd.load(f'{path_for_reference_structure}{reference_structure}.pdb')
    cmd.select(f'c. {chains_to_delete}')
    cmd.remove('sele')
    #load minimized filament WITHOUT TPM
    translation_output_per_cluster_path = f'{path_for_translation_output}cluster{i}/'
    cmd.load(f'{translation_output_per_cluster_path}{mutant}_base_filament_cluster{i}_MINIMIZED.pdb')


    #NOTE: THESE LOOPS ARE BASED ON THE GROMACS OUTPUT FROM MINIMIZING THE BARE ACTIN FILAMENTS
    #Basically, the minimized pdb doesn't assign chains to cation or PO4, but does for the nucleotide
    #So, the total number of assigned chains is double the number of actin protomers (16 protomers --> 32 total chains in my unaltered script)
    #We need to loop through every chain and re-assign the lettering/numbering of the ligands
    #the reason for this is that GROMACS gets confused unless everything is unique

    #for this first loop, go through all of the nucleotides and change the chain/numbering
    #set the gromacs assigned residue number (377 in my unaltered script)
    j = 377
    #set the starting residue number for what you'll assign (500 in my unaltered script)
    #We set it to a number that is high enough to avoid selecting protein residues
    #because actin only has 375 residues, 500 is way above that so selecting 'residue 500' won't select an actin residue
    k = 500
    #start the loop at 1 because index 0 is chain A, and the first ADP is chain B
    for l in range(1, 32, 2):
        #select chain B, residue 377, and assign it residue number 500
        #these values will change with each iteration of the loop to hit ALL nucleotides
        cmd.alter(selection=f'chain {actin_chains[l]} & resi {j}', expression=f"resi='{k}'")
        #change the chain to R. In my case, chain lettering stops at P so I will start ligand chains at R (skipping Q for a buffer)
        #We will assign all nucleotides to be in the same chain, but because they will have unique residue numbers GROMACS won't be confused
        cmd.alter(selection=f'chain {actin_chains[l]}', expression="chain='R'")
        #increase the residue number to match the next nucleotide
        j += 3
        #increase the assigned residue number
        k += 1
    j = 1
    # then change all the actin chains
    #start at C because we don't need to change A
    for l in range(2, 31, 2):
        #here, we are shifting the chain lettering down so that it's sequential
        #Chain C (index 2 in actin_chains) becomes chain B (index 1 in actin chains)
        #chain E (index 4 in actin chains) becomes chain C (index 2 in actin_chains) etc.
        cmd.alter(selection=f'chain {actin_chains[l]}', expression=f"chain='{actin_chains[j]}'")
        j += 1

    # then change all the cations
    #again, set the residue numbers far above anything else to avoid selection
    j = 600
    #set the loop to go through all of the residue numbers for the cations, as found in the minimized actin PDB file
    #starting at 376 going to 421 in my unaltered script
    #increase loop by 3 because gromacs assigns cation 376, nucleotide 377, and PO4 378
    #this would increase by just 2 in ADP or ATP nucleotide states
    if nucleotide_state == 'ADP-Pi':
        for l in range(376, 422, 3):
            #select anything with residue name of the cation and iterate through residue numbers
            #assign them all to chain S since that chain letter is unused
            cmd.alter(selection=f'resn {cation} & resi {l}', expression="chain='S'")
            #give each cation a unique number
            cmd.alter(selection=f'resn {cation} & resi {l}', expression=f"resi='{j}'")
            j += 1
    else:
        for l in range(376, 409, 2):
            cmd.alter(selection=f'resn {cation} & resi {l}', expression="chain='S'")
            cmd.alter(selection=f'resn {cation} & resi {l}', expression=f"resi='{j}'")
            j += 1


    # then change all the PO4 if ADP-Pi nucleotide state. Same logic as changing cation
    if nucleotide_state == 'ADP-Pi':
        j = 700
        for l in range(378, 424, 3):
            cmd.alter(selection=f'resn PO4 & resi {l}', expression="chain='T'")
            cmd.alter(selection=f'resn PO4 & resi {l}', expression=f"resi='{j}'")
            j += 1


    #here, I am using a 16 protomer long filament instead of 18 like in reference PDB 7UTL because of problems with 7UTL pointed end protomers
    #as a result, to add TPM to proper location I need to superimpose my filament with 7UTL filament but shifted 'upwards' 2 protomers (since chains A/B are barbed end)

    #first, select the first 2 actin chains and delete them
    cmd.select(f'{reference_structure} & c. {actin_chains[reference_starting_actin_chain_index]}:{actin_chains[reference_starting_actin_chain_index + 1]}')
    cmd.remove('sele')
    #now, do structural superposition between my filament and the remaining actin chains
    #WE DO THIS BECAUSE OUR FILAMENT IS SHORTER THAN THE REFERENCE FILAMENT.
    #NOW OUR MINIMIZED FILAMENT IS PROPERLY ALIGNED WITH THE BOUND TPM.
    cmd.super(f'{mutant}_base_filament_cluster{i}_MINIMIZED', f'{reference_structure}')

    #create an object with just TPM chains
    cmd.create("tropomyosin", f'{reference_structure} & c. {tpm_chains}')
    #remove the reference structure, leaving just our minimized actin and TPM
    cmd.delete(f'{reference_structure}')
    #save it as a PDB
    cmd.save(f'{translation_output_per_cluster_path}{mutant}_filament_cluster{i}.pdb')
    #delete everything so that we can re-load our PDB which will be just 1 combined object
    cmd.delete('all')
    cmd.load(f'{translation_output_per_cluster_path}{mutant}_filament_cluster{i}.pdb')
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
    cmd.alter_state(1, f'{mutant}_filament_cluster{i}', f"x,y,z=x-{v[12]},y-{v[13]},z-{v[14]}")
    # then we have to use the translation matrix again to apply a rotation as per https://pymolwiki.org/index.php/Transform_selection
    cmd.alter_state(1, f'{mutant}_filament_cluster{i}',
                    f"x,y,z=x*{v[0]}+y*{v[3]}+z*{v[6]},x*{v[1]}+y*{v[4]}+z*{v[7]},x*{v[2]}+y*{v[5]}+z*{v[8]}")
    cmd.save(f'{translation_output_per_cluster_path}{mutant}_filament_cluster{i}.pdb')
    cmd.delete('all')

#pull TPM out to 5A distance
for i in range(1, 6):
    translation_output_per_cluster_path = f'{path_for_translation_output}cluster{i}/'
    file_to_load = f'{translation_output_per_cluster_path}{mutant}_filament_cluster{i}.pdb'
    initial_distance(translation_output_per_cluster_path, mutant, Factin_or_Gactin, file_to_load, tpm_chains)

#translate TPM across full grid. Grid points determined from Viswanathan et al. M305L paper (supplementary data excel file)
#now, for positive and negative azimuthal rotations, we want to translate up and down actin filament
for i in range(1, 6):
    translation_output_per_cluster_path = f'{path_for_translation_output}cluster{i}/'
    cmd.delete('all')
    cmd.load(f'{translation_output_per_cluster_path}DOCKING_5A_{mutant}_{Factin_or_Gactin}0.00A0.00D.pdb')

    # positiveTranslations up filament
    translation = 2.5
    for s in np.arange(float(translation), 32.5, float(translation)):
        #maximum translation is +30A
        #translate up X-axis 2.5Angstroms, save, translate, save etc.
        #Because we moved the filament to the origin earlier, all of these translations will be along the same vector/axis
        cmd.translate(vector=[float(translation), 0, 0], selection=f'DOCKING_5A_{mutant}_{Factin_or_Gactin}0.00A0.00D & c. {tpm_chains}', camera=0)
        cmd.save(f'{translation_output_per_cluster_path}DOCKING_5A_{mutant}_{Factin_or_Gactin}%2.2fA0.00D.pdb' % s)

    cmd.delete('all')
    cmd.load(f'{translation_output_per_cluster_path}DOCKING_5A_{mutant}_{Factin_or_Gactin}0.00A0.00D.pdb')
    # negativeTranslations down filament
    #grid starts at -10A
    translation = -2.5
    for s in np.arange(float(translation), -12.5, float(translation)):
        cmd.translate(vector=[float(translation), 0, 0], selection=f'DOCKING_5A_{mutant}_{Factin_or_Gactin}0.00A0.00D & c. {tpm_chains}', camera=0)
        cmd.save(f'{translation_output_per_cluster_path}DOCKING_5A_{mutant}_{Factin_or_Gactin}%2.2fA0.00D.pdb' % s)


    # translate tpm for the full range of rotations, starting at 2.5 degrees up to 20 in increments of 2.5
for i in range(1, 6):
    #define 2 actin chains that we'll use to calculate vector through central axis
    actin_chain1 = 'I'
    actin_chain2 = 'N'
    translation_output_per_cluster_path = f'{path_for_translation_output}cluster{i}/'

    rotation_angle_plus = 2.5
    max_angle_plus = 5

    rotation_angle_minus = -2.5
    max_angle_minus = -35

    for a in range(-10, 35, 5):
        #load the translation structure
        cmd.load(f'{translation_output_per_cluster_path}DOCKING_5A_{mutant}_{Factin_or_Gactin}%2.2fA0.00D.pdb' % a)
        #isolate TPM chains
        cmd.extract('tropomyosin', selection=f'DOCKING_5A_{mutant}_{Factin_or_Gactin}%2.2fA0.00D & chain {tpm_chains}' % a)
        #positive rotations
        rotation_matrix = calculate_rotation_matrix(actin_chain1, actin_chain2, rotation_angle_plus)
        for t in np.arange(float(rotation_angle_plus), float(max_angle_plus) + float(rotation_angle_plus), float(rotation_angle_plus)):
            cmd.transform_object('tropomyosin', rotation_matrix)
            cmd.save(f'{translation_output_per_cluster_path}DOCKING_5A_{mutant}_{Factin_or_Gactin}%2.2fA%2.2fD.pdb' % (a, t))
        cmd.delete('all')

        max_angle_plus += 2.50

        # load the translation structure
        cmd.load(f'{translation_output_per_cluster_path}DOCKING_5A_{mutant}_{Factin_or_Gactin}%2.2fA0.00D.pdb' % a)
        # isolate TPM chains
        cmd.extract('tropomyosin', selection=f'DOCKING_5A_{mutant}_{Factin_or_Gactin}%2.2fA0.00D & chain {tpm_chains}' % a)
        # negative rotations
        rotation_matrix = calculate_rotation_matrix(actin_chain1, actin_chain2, rotation_angle_minus)
        for t in np.arange(float(rotation_angle_minus), float(max_angle_minus) + float(rotation_angle_minus), float(rotation_angle_minus)):
            cmd.transform_object('tropomyosin', rotation_matrix)
            cmd.save(f'{translation_output_per_cluster_path}DOCKING_5A_{mutant}_{Factin_or_Gactin}%2.2fA%2.2fD.pdb' % (a, t))
        cmd.delete('all')

        max_angle_minus += 2.50


    starting_angle_plus = 1.25
    rotation_angle_plus = 2.50
    max_angle_plus = 6.25

    starting_angle_minus = -1.25
    rotation_angle_minus = -2.50
    max_angle_minus = -33.75

    for a in np.arange(-7.5, 32.5, 5):
        # load the translation structure
        cmd.load(f'{translation_output_per_cluster_path}DOCKING_5A_{mutant}_{Factin_or_Gactin}%2.2fA0.00D.pdb' % a)
        # isolate TPM chains
        cmd.extract('tropomyosin', selection=f'DOCKING_5A_{mutant}_{Factin_or_Gactin}%2.2fA0.00D & chain {tpm_chains}' % a)
        # define the start angle
        rotation_matrix = calculate_rotation_matrix(actin_chain1, actin_chain2, starting_angle_plus)
        # move to the starting angle
        cmd.transform_object('tropomyosin', rotation_matrix)
        cmd.save(f'{translation_output_per_cluster_path}DOCKING_5A_{mutant}_{Factin_or_Gactin}%2.2fA{starting_angle_plus}D.pdb' % a)
        # positive rotations
        rotation_matrix = calculate_rotation_matrix(actin_chain1, actin_chain2, rotation_angle_plus)
        for t in np.arange(float(starting_angle_plus) + float(rotation_angle_plus), float(max_angle_plus) + float(rotation_angle_plus), float(rotation_angle_plus)):
            cmd.transform_object('tropomyosin', rotation_matrix)
            cmd.save(f'{translation_output_per_cluster_path}DOCKING_5A_{mutant}_{Factin_or_Gactin}%2.2fA%2.2fD.pdb' % (a, t))
        cmd.delete('all')

        max_angle_plus += 2.50


        # load the translation structure
        cmd.load(f'{translation_output_per_cluster_path}DOCKING_5A_{mutant}_{Factin_or_Gactin}%2.2fA0.00D.pdb' % a)
        # isolate TPM chains
        cmd.extract('tropomyosin', selection=f'DOCKING_5A_{mutant}_{Factin_or_Gactin}%2.2fA0.00D & chain {tpm_chains}' % a)
        # define the start angle
        rotation_matrix = calculate_rotation_matrix(actin_chain1, actin_chain2, starting_angle_minus)
        # move to the starting angle
        cmd.transform_object('tropomyosin', rotation_matrix)
        cmd.save(f'{translation_output_per_cluster_path}DOCKING_5A_{mutant}_{Factin_or_Gactin}%2.2fA{starting_angle_minus}D.pdb' % a)
        # positive rotations
        rotation_matrix = calculate_rotation_matrix(actin_chain1, actin_chain2, rotation_angle_minus)
        for t in np.arange(float(starting_angle_minus) + float(rotation_angle_minus), float(max_angle_minus) + float(rotation_angle_minus), float(rotation_angle_minus)):
            cmd.transform_object('tropomyosin', rotation_matrix)
            cmd.save(f'{translation_output_per_cluster_path}DOCKING_5A_{mutant}_{Factin_or_Gactin}%2.2fA%2.2fD.pdb' % (a, t))
        cmd.delete('all')

        max_angle_minus += 2.50
