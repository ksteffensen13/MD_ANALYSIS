from pymol import cmd
from pymol import math
from pymol import stored
import numpy as np
import sys

#define the function for docking
def docking(path_for_translation_output, mutant, Factin_or_Gactin, tpm_chains, start_distance, end_distance, docking_increment, file_to_dock, docked_output_file):
    cmd.delete('all')
    #load the reference structure so that we can use the full TPM for calculating translation vectors
    cmd.load(f'{path_for_reference_structure}{reference_structure}.pdb')
    cmd.select(f'c. {chains_to_delete}')
    cmd.remove('sele')
    cmd.select(f'{reference_structure} & c. {reference_starting_actin_chain}:{reference_ending_actin_chain}')
    cmd.remove('sele')
    #load the minimized PDB file from the previous TPM distance
    #the file is loaded from folder containing all of the minimized PDBs
    cmd.load(f'{file_to_dock}.pdb')
    cmd.super(f'{reference_structure}', f'{file_to_dock}')
    #isolate the chain identifiers for tpm chains into separate variables (needed for calculating vectors)
    tpm_chain_list = tpm_chains.split('+')
    tpm_chain1 = tpm_chain_list[0]
    tpm_chain2 = tpm_chain_list[1]
    #calculate a vector through the central actin filament axis
    v1 = draw_axis('chain I', 'chain N')
    #calculate a vector through the center of the 2 COMPLETE tpm chains
    #Note: calculating on the trimmed TPM chains creates a less accurate vector
    v2 = draw_axis(f'{reference_structure} & chain {tpm_chain1}', f'{reference_structure} & chain {tpm_chain2}')
    #calculate the vector from v1 through v2 by calculating the crossproduct
    v3 = np.cross(v1, v2)
    #calculate normalized vector (with a magnitude of 1)
    #need normalized vector so that it scales with our translationd istance
    v3_normalized = v3 / np.linalg.norm(v3)
    # create a transformation matrix for the calculated vector and the desired initial translation distance
    move_vector = [1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, docking_increment * v3_normalized[0], docking_increment * v3_normalized[1], docking_increment * v3_normalized[2], 1.0]
    cmd.remove(f'{reference_structure}')
    #update the chains after minimization
    cmd.alter(selection='chain P', expression="chain='i'")
    cmd.alter(selection='chain Q', expression="chain='j'")
    cmd.alter(selection='resn ADP', expression="chain='R'")
    cmd.alter(selection='resn MG', expression="chain='S'")
    cmd.alter(selection='resn PO4', expression="chain='T'")
    #move the TPM chains along the desginated vector/distance
    cmd.transform_selection(f'chain {tpm_chains}', move_vector)
    #save the file
    cmd.save(f'{docked_output_file}.pdb')


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



#####################################################END OF FUNCTIONS#############################################




# Extract command-line arguments. Here, we can set the start/end distance/other values for the docking function every time we call this script
#the line below will just call variables from within the UNIX script running this pythons cript
if len(sys.argv) != 14:
    print("Usage: python tpm_docking.py <MUTANT> <Factin_or_Gactin> <REFERENCE_STRUCTURE> <CHAINS_TO_DELETE> <TPM_CHAINS> <PATH_FOR_TRANSLATION_OUTPUT> <PATH_FOR_PYMOL_AXES_SCRIPT> <PATH_FOR_ROTATION_AXES_SCRIPT> <START_DISTANCE> <END_DISTANCE> <DOCKING_INCREMENT> <reference_starting_actin_chain>")
    sys.exit(1)

#number 0 is the script name itself, all the variables are 'system argument #1, 2, 3, etc."
mutant = sys.argv[1]
Factin_or_Gactin = sys.argv[2]
reference_structure = sys.argv[3]
path_for_reference_structure = sys.argv[4]
chains_to_delete = sys.argv[5]
tpm_chains = sys.argv[6]
path_for_translation_output = sys.argv[7]
path_for_pymol_axes_script = sys.argv[8]
path_for_rotation_axes_script = sys.argv[9]
start_distance = float(sys.argv[10])
end_distance = float(sys.argv[11])
docking_increment = float(sys.argv[12])
reference_starting_actin_chain = sys.argv[13]

if reference_structure == '7UTL':
    state = "BLOCKED"
    reference_ending_actin_chain = 'R'
if reference_structure == '7UTI':
    state = "CLOSED"
    reference_ending_actin_chain = 'R'
if reference_structure == '8EFI':
    state = "OPEN"
    reference_ending_actin_chain = 'R'

for i in range(1, 6):
    super f'WT{i}', 7UTL_blocked
#for clusters 1 to 5 (this loop doesn't include 6)
for i in range(1, 6):
    translation_output_per_cluster_path = f'{path_for_translation_output}cluster{i}/'
    
    
    #First, define the azimuthal rotation coordinates to help setup the loop
    #the starting angle (for the whole number translations, it's 0 degrees)
    starting_angle_plus = 0
    #how much each azimuthal coordinate is separated by
    rotation_angle_plus = 2.5
    #the max angle of the first translation point we'll start with (for -10A, it's 5 degrees)
    max_angle_plus = 5
    
    #same thing as above but for negative azimuthal rotations
    rotation_angle_minus = -2.5
    max_angle_minus = -35
    for s in range(-10, 35, 5):
        #set loop for + rotations
        for t in np.arange(float(starting_angle_plus), float(max_angle_plus) + float(rotation_angle_plus), float(rotation_angle_plus)):
            file_to_dock = f'{path_for_translation_output}MINIMIZED/{start_distance}A/MINIMIZED_DOCKING_{start_distance}A_{mutant}_{Factin_or_Gactin}%2.2fA%2.2fD.pdb' % (s, t)
            docked_output_file = f'{path_for_translation_output}DOCKING_{end_distance}A_{mutant}_{Factin_or_Gactin}%2.2fA%2.2fD.pdb' % (s, t)
            #run the function to dock
            docking(translation_output_per_cluster_path, mutant, Factin_or_Gactin, tpm_chains, start_distance, end_distance, docking_increment, file_to_dock, docked_output_file)
        #update the max angle for the next loop iteration
        max_angle_plus += 2.50
        
        for t in np.arange(float(rotation_angle_minus), float(max_angle_minus) + float(rotation_angle_minus), float(rotation_angle_minus)):
            file_to_dock = f'{path_for_translation_output}MINIMIZED/{start_distance}A/MINIMIZED_DOCKING_{start_distance}A_{mutant}_{Factin_or_Gactin}%2.2fA%2.2fD.pdb' % (s, t)
            docked_output_file = f'{path_for_translation_output}DOCKING_{end_distance}A_{mutant}_{Factin_or_Gactin}%2.2fA%2.2fD.pdb' % (s, t)
            docking(translation_output_per_cluster_path, mutant, Factin_or_Gactin, tpm_chains, start_distance, end_distance, docking_increment, file_to_dock, docked_output_file)
        max_angle_minus += 2.50
    
    
    
    #repeat the process for the decimal translation points which start at 1.25 and -1.25 degrees instead of 0
    starting_angle_plus = 1.25
    rotation_angle_plus = 2.50
    max_angle_plus = 6.25

    starting_angle_minus = -1.25
    rotation_angle_minus = -2.50
    max_angle_minus = -33.75
    for s in np.arange(-7.5, 32.5, 5):
        for t in np.arange(float(starting_angle_plus), float(max_angle_plus) + float(rotation_angle_plus), float(rotation_angle_plus)):
            file_to_dock = f'{path_for_translation_output}MINIMIZED/{start_distance}A/MINIMIZED_DOCKING_{start_distance}A_{mutant}_{Factin_or_Gactin}%2.2fA%2.2fD.pdb' % (s, t)
            docked_output_file = f'{path_for_translation_output}DOCKING_{end_distance}A_{mutant}_{Factin_or_Gactin}%2.2fA%2.2fD.pdb' % (s, t)
            docking(translation_output_per_cluster_path, mutant, Factin_or_Gactin, tpm_chains, start_distance, end_distance, docking_increment, file_to_dock, docked_output_file)
        max_angle_plus += 2.50
        
        for t in np.arange(float(starting_angle_minus), float(max_angle_minus) + float(rotation_angle_minus), float(rotation_angle_minus)):
            file_to_dock = f'{path_for_translation_output}MINIMIZED/{start_distance}A/MINIMIZED_DOCKING_{start_distance}A_{mutant}_{Factin_or_Gactin}%2.2fA%2.2fD.pdb' % (s, t)
            docked_output_file = f'{path_for_translation_output}DOCKING_{end_distance}A_{mutant}_{Factin_or_Gactin}%2.2fA%2.2fD.pdb' % (s, t)
            docking(translation_output_per_cluster_path, mutant, Factin_or_Gactin, tpm_chains, start_distance, end_distance, docking_increment, file_to_dock, docked_output_file)
        max_angle_minus += 2.50
    
