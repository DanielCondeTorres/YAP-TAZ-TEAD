#!/bin/bash

# Transform .pdb from Protein Data Bank to .itp GROMACS format




gmx make_ndx -f prod_final.tpr -o mem.ndx << EOF
1|19|18
q
EOF
