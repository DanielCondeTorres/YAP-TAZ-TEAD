#!/bin/bash

COMPLEX_PDB_FILE_1="$1"  # Archivo PDB del complejo
gmx solvate -cp ../Output/"$COMPLEX_PDB_FILE_1" -cs spc216.gro -o ../Output/complex_water_box.pdb -p ../Output/topol.top
# Add ions
gmx grompp -f ../Input_files/MDP_FILES/ions.mdp -c ../Output/complex_water_box.pdb -p ../Output/topol.top -o ../Output/ions.tpr -maxwarn 1
echo "SOL" | gmx genion -s ../Output/ions.tpr -o ../Output/complex_solv_ions.gro -p ../Output/topol.top -pname NA -nname CL -neutral
# Minimization of energy
gmx grompp -f ../Input_files/MDP_FILES/min.mdp -c ../Output/complex_solv_ions.gro -p ../Output/topol.top -o ../Output/min.tpr -maxwarn 2
gmx mdrun -v -deffnm ../Output/min -ntmpi 4


# FEP equilibration NVT
gmx grompp -f ../Input_files/MDP_FILES/eq_nvt_fep.mdp -c ../Output/min.gro  -r ../Output/min.gro -p ../Output/topol.top -o ../Output/eq_fes.tpr -maxwarn 2
gmx mdrun -v -deffnm ../Output/eq_fes -ntmpi 4

# Equilibración NPT
gmx grompp -f ../Input_files/MDP_FILES/eq_npt.mdp -c ../Output/eq_fes.gro -p ../Output/topol.top -o ../Output/eq_2.tpr -r ../Output/eq_fes.gro -maxwarn 2
gmx mdrun -v -deffnm ../Output/eq_2 -ntmpi 4