# LAST VERSION AS OF 09.10.2020

# ---------
# Z-MATRIX
# ---------   

from chemml.chem import Molecule
import chemcoord as cc
import pandas as pd
import numpy as np
import time
import os

# The list all_molecules, created along with the XYZ coordintes, contains the names of all .XYZ files, including
# the the counting for isomers. All these names are used to read the XYZ coordinates and create the
# appropriate Z-matrix.

# Creates the folder to store the Z-matrix representation of the geometry.
os.mkdir('geom_zmax')

failed_zmax = {'SMILES'  : [], 'Formula' : []}
z_matrix_all = []

for name, smiles in zip(all_molecules, all_smiles):
    
    try:
        # Gets into the Z-matrix folder.
        os.chdir('geom_zmax')
        
        atoms = []
        bonds = ['']
        number_bonds = ['']
        angles = ['','']
        number_angles = ['','']
        dihedrals = ['','','']
        number_dihedrals = ['','','']
    
        # Read XYZ cooridinates and transforms them into z-matrix representation.
        # The z-matrix is saved as a DataFrame
        name_zmat = cc.Cartesian.read_xyz('../geom_xyz/'+name+'.xyz', start_index=1).get_zmat()
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
                
        z_matrix_all.append(name)

    except:
        failed_zmax['SMILES'].append(smiles)
        failed_zmax['Formula'].append(name)
    
    os.chdir('../')

# Writes a molecules.txt file listing all molecules for which geometries have been generated.
# This file is needed for running the input files templates from the run_gaussian.sh code.
for formula in z_matrix_all:
    print(formula, file = open('molecules.txt', 'a'))
