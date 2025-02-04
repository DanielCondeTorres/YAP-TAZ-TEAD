#!/bin/bash

# Transform .pdb from Protein Data Bank to .itp GROMACS format

COMPLEX_PDB_FILE="$1"  # Recibe el argumento del Makefile
FORCEFIELD_FILE="$2"  # Recibe el argumento del Makefile



GMXLIB=../Input_files/FORCEFIELDS gmx pdb2gmx -f  $COMPLEX_PDB_FILE -o ../Output/complex.gro -ignh -p ../Output/ -i ../Output/
mv ../Output/.top ../Output/topol.top
# Create water box
gmx editconf -f ../Output/complex.gro -o ../Output/complex_water_box.gro -c -d 1.8 -bt cubic
gmx solvate -cp ../Output/complex_water_box.gro -cs spc216.gro -o ../Output/complex_water_box.gro -p ../Output/topol.top
rm ../Output/#* 
# Add ions
gmx grompp -f ../Input_files/MDP_FILES/ions.mdp -c ../Output/complex_water_box.gro -p ../Output/topol.top -o ../Output/ions.tpr -maxwarn 1
gmx genion -s ../Output/ions.tpr -o ../Output/complex_solv_ions.gro -p  ../Output/topol.top -pname NA -nname CL -neutral << EOF
SOL
EOF

# Energy minimization
gmx grompp -f ../Input_files/MDP_FILES/min.mdp -c ../Output/complex_solv_ions.gro -p ../Output/topol.top -o ../Output/min.tpr -maxwarn 2
gmx mdrun -v -deffnm  ../Output/min -ntmpi 4

# NPT equilibration
gmx grompp -f ../Input_files/MDP_FILES/eq.mdp -c ../Output/min.gro -p ../Output/topol.top -o ../Output/eq.tpr -maxwarn 2
gmx mdrun -v -deffnm  ../Output/eq -ntmpi 4

# Create .tpr to production
gmx grompp -f ../Input_files/MDP_FILES/prod.mdp -c ../Output/eq.gro -p ../Output/topol.top -o ../Output/prod.tpr -maxwarn 2
