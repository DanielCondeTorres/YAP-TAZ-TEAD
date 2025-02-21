#!/bin/bash
LIGAND="$1" 
PDB_INPUT="$2"
# Lee el valor de MODE desde mode_output.txt (asegúrate de que NO haya espacios)
MODE_VALUE=$(sed '1q' mode_output.txt)
echo "Valor de MODE_VALUE: '${MODE_VALUE}'"

# Extrae la mejor conformación
grep -v BRANCH all.pdbqt | grep -v ROOT | grep -v TORSDOF > docking_results.pdb
awk "/^MODEL ${MODE_VALUE}/,/^ENDMDL/" docking_results.pdb > best_conformation.pdb

# Ejecuta el script de formación del complejo
gmx editconf -f best_conformation.pdb -o best_conformation2.pdb

#Necesito pornerle variables!!!! que es la de verteporfin_gmx22.pdb 
pymol -cq -d "load "$LIGAND", ligand; load best_conformation2.pdb, template; align ligand, template; save best_conformation2.pdb, ligand; quit;"
python complex_ligand_formation.py -receptor receptor_docking.pdb -ligand best_conformation2.pdb
# Creation of posrte.itp file for ligand
echo "Other" | gmx genrestr -f best_conformation2.pdb -o posre.itp -fc 1000 1000 1000
head -n 5 e"$PDB_INPUT" > complex_ligand.pdb
cat best_conformation2.pdb >> complex_ligand.pdb

