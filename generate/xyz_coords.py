# LAST VERSION AS OF 09.10.2020

# ----------------
# GET XYZ COORDS:
# ----------------

from chemml.chem import Molecule
import chemcoord as cc
import pandas as pd
import numpy as np
import time
import os

# Creates folder for storing the XYZ geometries.
os.mkdir('geom_xyz')

test_case = pd.read_csv('se_mol.csv')

# Dictionary storing the smiles codes and formulas of those molecules for which the Python script doesn't work.
failed_jobs = {'SMILES'  : [], 'Formula' : []}

# The first list will store the name of the molecules that face problems with this program, the second sotres
# molecules with no problem and the last one the smiles codes.
#all_molecules_prob = []
all_molecules = []
all_smiles = []

for smiles, formula in zip(test_case['SMILES'], test_case['Formula']):
    #print(smiles,formula)
    
    # Tries to run the jobs and skips those for which the script doesn't work.
    try:
        
        # Gets into the XYZ folder.
        os.chdir('geom_xyz')
        
        # These three lines read the SMILES input, add the hydrogens to structure and
        # optimise the geomtry to obtain the XYZ coordintes.
        mol = Molecule(smiles, input_type = 'smiles')
        mol.hydrogens('add')
        mol.to_xyz('UFF')
        
        # This final array merges and organises the atomic symbols and geometries properly
        # as the initial versions are separeted.
        final_array = np.column_stack((mol.xyz.atomic_symbols, mol.xyz.geometry))
        #print(final_array)
        
        # Anna's contribution.
        # Prints the geometries into XYZ files named after the molecular formula for each compound, considering isomers.
        # https://stackoverflow.com/questions/5914627/prepend-line-to-beginning-of-a-file
        # https://stackoverflow.com/questions/12309976/how-do-i-convert-a-list-into-a-string-with-spaces-in-python/12309982
    
        # If the molecular formula is unique, the XYZ file will be naed after the formula (e.g. PH3.xyz)
        if formula not in all_molecules:
            isomers = 0
            all_molecules.append(formula)
            all_smiles.append(smiles)
            with open(formula+'.xyz', 'a') as inp:
                inp.write(str('atom'+'\t x'+'\t\t\t y'+'\t\t\t z')+'\n\n')
                for row in final_array:
                    print('\t'.join(map(str, row)), file = inp)
        
        # If there are multiple entries with the same formula, i.e. isomers, the XYZ file will be named according
        # with the number of isomers present, e.g. C2H5P.xyz, C2H5P_1.xyz, C2H5P_2.xyz, etc.
        else:
            isomers = isomers + 1
            all_molecules.append(formula+'_'+str(isomers))
            all_smiles.append(smiles)
            with open(formula+'_'+str(isomers)+'.xyz', 'a') as inp_iso:
                inp_iso.write(str('atom'+'\t x'+'\t\t\t y'+'\t\t\t z')+'\n\n')
                for row in final_array:
                    print('\t'.join(map(str, row)), file = inp_iso)

    except:
        failed_jobs['SMILES'].append(smiles)
        failed_jobs['Formula'].append(formula)

    os.chdir('../')

# Writes the molecular_reference.txt file that contains the file names along with the SMILES code, in order to
# have track of the molecular identities.
for file_name, smiles_code in zip(all_molecules, all_smiles):
    print(file_name+','+smiles_code, file = open('molecular_reference.txt', 'a'))
