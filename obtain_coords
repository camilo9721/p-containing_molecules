# CODE FOR OBTAINING XYZ COORDINATES AND Z-MATRIX

from chemml.chem import Molecule
import chemcoord as cc
import pandas as pd
import numpy as np
import time

# ------------------
# XYZ COORDINTES
# ------------------

# Loads the csv file
test_case = pd.read_csv('test_case.csv')

for smiles, formula in zip(test_case['SMILES'], test_case['Formula']):
    
    # These three lines read the SMILES input, add the hydrogens to the structure and
    # optimise the geomtry to obtain the xyz coordintes.
    mol = Molecule(smiles, input_type = 'smiles')
    mol.hydrogens('add')
    mol.to_xyz('UFF')
    
    # This final array merges and organises the atomic symbols and geometries properly
    # as the initial versions are separeted
    final_array = np.column_stack((mol.xyz.atomic_symbols, mol.xyz.geometry))
    
    # Prints the geometries into xyz files named after the molecular formula for each compound
    # Anna's contribution
    # https://stackoverflow.com/questions/5914627/prepend-line-to-beginning-of-a-file
    # https://stackoverflow.com/questions/12309976/how-do-i-convert-a-list-into-a-string-with-spaces-in-python/12309982
    with open(formula+'.xyz', 'a') as inp:
        inp.write(str('atom'+'\t x'+'\t\t\t y'+'\t\t\t z')+'\n\n')
        for row in final_array:
            print('\t'.join(map(str, row)), file = inp)

# ---------
# Z-MATRIX
# ---------

for name in test_case['Formula']:
    
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
for formula in test_case['Formula']:
    print(formula, file = open('molecules.txt', 'a'))
