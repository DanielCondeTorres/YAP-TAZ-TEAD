#!/bin/bash

# Transform .pdb from Protein Data Bank to .itp GROMACS format

COMPLEX_PDB_FILE="$1"  # Recibe el argumento del Makefile

gmx pdb2gmx -f  $COMPLEX_PDB_FILE -o complex.gro -water spce



gmx make_ndx -f prod_final.tpr -o mem.ndx << EOF
1|19|18
q
EOF
