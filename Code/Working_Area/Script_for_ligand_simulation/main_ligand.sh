#!/bin/bash

COMPLEX_PDB_FILE_1="$1"  # Archivo PDB del complejo
OUTPUT_SAVE="$2"
gmx solvate -cp "$COMPLEX_PDB_FILE_1" -cs spc216.gro -o "$OUTPUT_SAVE"complex_water_box.pdb -p "$OUTPUT_SAVE"topol_mod2.top
# Add ions


gmx grompp -f ../Input_files/MDP_FILES/ions.mdp -c "$OUTPUT_SAVE"complex_water_box.pdb -p "$OUTPUT_SAVE"topol_mod2.top -o "$OUTPUT_SAVE"ions.tpr -maxwarn 1
echo "SOL" | gmx genion -s "$OUTPUT_SAVE"ions.tpr -o "$OUTPUT_SAVE"complex_solv_ions.gro -p "$OUTPUT_SAVE"topol_mod2.top -pname NA -nname CL -neutral
# Minimization of energy
gmx grompp -f ../Input_files/MDP_FILES/min.mdp -c "$OUTPUT_SAVE"complex_solv_ions.gro -p "$OUTPUT_SAVE"topol_mod2.top -o "$OUTPUT_SAVE"min.tpr -maxwarn 2
gmx mdrun -v -deffnm "$OUTPUT_SAVE"min -ntmpi 4


# FEP equilibration NVT
#gmx grompp -f ../Input_files/MDP_FILES/eq_nvt_fep.mdp -c ../Output/min.gro  -r ../Output/min.gro -p ../Output/topol_mod2.top -o ../Output/eq_fes.tpr -maxwarn 2
#gmx mdrun -v -deffnm ../Output/eq_fes -ntmpi 4
# New minimization
#gmx grompp -f ../Input_files/MDP_FILES/min.mdp -c ../Output/eq_fes.gro -p ../Output/topol_mod2.top -o ../Output/min_fes.tpr -maxwarn 2
#gmx mdrun -v -deffnm ../Output/min_fes -ntmpi 4
# Equilibración NPT 1
#gmx grompp -f ../Input_files/MDP_FILES/eq_npt_fep_1.mdp -c ../Output/min_fes.gro -p ../Output/topol_mod2.top -o ../Output/eq_npt_fep_1.tpr -r ../Output/min_fes.gro -maxwarn 2
#gmx mdrun -v -deffnm ../Output/eq_npt_fep_1 -ntmpi 4
# Equilibración NPT 2
#gmx grompp -f ../Input_files/MDP_FILES/eq_npt_fep_2.mdp -c ../Output/eq_npt_fep_1.gro -p ../Output/topol_mod2.top -o ../Output/eq_npt_fep_2.tpr -maxwarn 2
#gmx mdrun -v -deffnm ../Output/eq_npt_fep_2 -ntmpi 4
#Create prod.tpr
#gmx grompp -f ../Input_files/MDP_FILES/prod.mdp -c ../Output/eq_npt_fep_2.gro -p ../Output/topol_mod2.top -o ../Output/prod.tpr -maxwarn 2
