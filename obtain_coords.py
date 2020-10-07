# LAST VERSION OF THE CODE FOR OBTAINING XYZ COORDINATES AND Z-MATRIX 07.10.2020

from chemml.chem import Molecule
import chemcoord as cc
import pandas as pd
import numpy as np
import time

# ------------------
# XYZ COORDINTES
# ------------------

test_case = pd.read_csv('test_case.csv')
all_molecules = []

for smiles, formula in zip(test_case['SMILES'], test_case['Formula']):
    
    # These three lines read the SMILES input, add the hydrogens to structure and
    # optimise the geomtry to obtain the XYZ coordintes.
    mol = Molecule(smiles, input_type = 'smiles')
    mol.hydrogens('add')
    mol.to_xyz('UFF')
    
    # This final array merges and organises the atomic symbols and geometries properly
    # as the initial versions are separeted.
    final_array = np.column_stack((mol.xyz.atomic_symbols, mol.xyz.geometry))
    
    # Anna's contribution.
    # Prints the geometries into XYZ files named after the molecular formula for each compound, considering isomers.
    # https://stackoverflow.com/questions/5914627/prepend-line-to-beginning-of-a-file
    # https://stackoverflow.com/questions/12309976/how-do-i-convert-a-list-into-a-string-with-spaces-in-python/12309982
    
    # If the molecular formula is unique, the XYZ file will be naed after the formula (e.g. PH3.xyz)
    if formula not in all_molecules:
        isomers = 0
        all_molecules.append(formula)
        with open(formula+'.xyz', 'a') as inp:
            inp.write(str('atom'+'\t x'+'\t\t\t y'+'\t\t\t z')+'\n\n')
            for row in final_array:
                print('\t'.join(map(str, row)), file = inp)

    # If there are multiple entries with the same formula, i.e. isomers, the XYZ file will be named according with
    # the number of isomers present, e.g. C2H5P.xyz, C2H5P_1.xyz, C2H5P_2.xyz, etc.
    else:
        isomers = isomers + 1
        all_molecules.append(formula+'_'+str(isomers))
        with open(formula+'_'+str(isomers)+'.xyz', 'a') as inp:
            inp.write(str('atom'+'\t x'+'\t\t\t y'+'\t\t\t z')+'\n\n')
            for row in final_array:
                print('\t'.join(map(str, row)), file = inp)

# Writes the molecular_reference.txt file that contains the file names along with the SMILES code, in order to
# have track of the molecular identities.
for file_name, smiles_code in zip(all_molecules, test_case['SMILES']):
    print(file_name+','+smiles_code, file = open('molecular_reference.txt', 'a'))
                
# ---------
# Z-MATRIX
# ---------                

# The list all_molecules, created along with the XYZ coordintes, contains the names of all .XYZ files, including
# the the counting for isomers. All these names are used to read the XYZ coordinates and create the
# appropriate Z-matrix.

for name in all_molecules:
    
    atoms = []
    bonds = ['']
    number_bonds = ['']
    angles = ['','']
    number_angles = ['','']
    dihedrals = ['','','']
    number_dihedrals = ['','','']
    
    # Read XYZ cooridinates and transforms them into z-matrix representation.
    # The z-matrix is saved as a DataFrame
    name_zmat = cc.Cartesian.read_xyz(name+'.xyz', start_index=1).get_zmat()
    # Anna's contribution to list the data frame properly
    name_zmat_new = name_zmat.change_numbering(new_index=list(range(1, len(name_zmat) + 1)))
    
    # The name_zmat contains some extra information.
    # These for loops are meant to store the important information from the z-matrix representation from above.
    for element in name_zmat_new['atom']:
        atoms.append(element)
    
    for b, bond in zip(name_zmat_new['b'],name_zmat_new['bond']):
        if isinstance(b, str) == False:
            number_bonds.append(b)
            bonds.append(bond)
            
    for a, angle in zip(name_zmat_new['a'],name_zmat_new['angle']):
        if isinstance(a, str) == False:
            number_angles.append(a)
            angles.append(angle)

    for c, dihedral in zip(name_zmat_new['d'],name_zmat_new['dihedral']):
        if isinstance(c, str) == False:
            number_dihedrals.append(c)
            dihedrals.append(dihedral)
    
    # The z-matrix is printed in files named after the molecular formula for each compound.
    # The charge and multiplicity are included at the top of the file.
    with open(name+'.zmax', 'a') as f:
        f.write(str('0'+' '+'1')+'\n')
        for row in zip(atoms, number_bonds, bonds, number_angles, angles, number_dihedrals, dihedrals):
            print('\t'.join(map(str, row)), file = f)

# Writes a molecules.txt file listing all molecules for which geometries have been generated.
# This file is needed for running the input files templates from the run_gaussian.sh code.
for formula in all_molecules:
    print(formula, file = open('molecules.txt', 'a'))
