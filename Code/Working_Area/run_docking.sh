#!/bin/bash

# Lee el valor de MODE desde mode_output.txt (asegúrate de que NO haya espacios)
MODE_VALUE=$(sed '1q' mode_output.txt)
echo "Valor de MODE_VALUE: '${MODE_VALUE}'"

# Extrae la mejor conformación
grep -v BRANCH all.pdbqt | grep -v ROOT | grep -v TORSDOF > docking_results.pdb
awk "/^MODEL ${MODE_VALUE}/,/^ENDMDL/" docking_results.pdb > best_conformation.pdb

# Ejecuta el script de formación del complejo
python complex_ligand_formation.py -receptor receptor_docking.pdb -ligand best_conformation.pdb
