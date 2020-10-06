# P-containing Molecules Project

This repository contains different Python scripts for different tasks needed calcualting vibrational spectra of p-containing molecules.
The scrips included are the follwonig:
  - obtain_coords.py: This script reads the SMILES input for each molecule and obtains its respective XYZ and Z-matrix coordintes separetly. The outputs are named                         accordingly with the molecule's name, e.g. PH3.xyz and PH3.zmax.
  - create_input.py:  This script allows to generate Gaussian input files systematically. The script needs the files: molecules.txt (which lists all molecules                             considered), harmonic_def2.gjf (input file template), subtem.pbs (submission file template) and a folder containing all geometries in Z-matrix                       coordinates.
