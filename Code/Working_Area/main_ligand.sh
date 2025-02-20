#!/bin/bash

COMPLEX_PDB_FILE_1="$1"  # Archivo PDB del complejo
gmx solvate -cp ../Output/"$COMPLEX_PDB_FILE_1" -cs spc216.gro -o ../Output/complex_water_box.pdb -p ../Output/topol.top
gmx grompp -f ../Input_files/MDP_FILES/ions.mdp -c ../Output/complex_water_box.gro -p ../Output/topol.top -o ../Output/ions.tpr -maxwarn 1