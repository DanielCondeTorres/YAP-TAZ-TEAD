#!/bin/bash


echo "All inputs received: $@"
# Recibe los argumentos del Makefile
COMPLEX_PDB_FILE_1="$1"  # Archivo PDB del complejo
COMPLEX_PDB_FILE_2="$2"  # Nombre del archivo PDB del complejo


# Verifica si Order es igual a 1
# Crear directorio de salida
mkdir -p ../Output/

# Convertir PDB a GRO y generar topología
grep -v HOH "$COMPLEX_PDB_FILE_1" > ../Output/clean_1.pdb
GMXLIB=../Input_files/FORCEFIELDS gmx pdb2gmx -f ../Output/clean_1.pdb -o ../Output/complex.gro -ignh -p ../Output/topol.top -i ../Output/posre.itp 

# Crear caja de agua
gmx editconf -f ../Output/complex.gro -o ../Output/complex.gro -c -d 2.8 -bt cubic


#Añadir segundo complejo
if [ -n "$COMPLEX_PDB_FILE_2" ]; then
    echo "Processing second PDB file: $COMPLEX_PDB_FILE_2"
    sed -n '/^\[ moleculetype \]/,/^; Include water topology/p' ../Output/topol.top | sed '/^; Include water topology/d' > ../Output/molecule.itp
    grep -v HOH "$COMPLEX_PDB_FILE_2" > ../Output/clean_2.pdb
    GMXLIB=../Input_files/FORCEFIELDS gmx pdb2gmx -f ../Output/clean_2.pdb -o ../Output/complex_2.gro -ignh -p ../Output/topol_2.top -i ../Output/posre_2.itp
    sed -i '' 's/^Protein_chain_A/Protein_chain_B/' ../Output/topol_2.top
    sed -n '/^\[ moleculetype \]/,/^; Include water topology/p' ../Output/topol_2.top | sed '/^; Include water topology/d' > ../Output/molecule_2.itp

    cat << 'EOF' > ../Output/topol.top
; Include forcefield parameters
#include "../Input_files/FORCEFIELDS/charmm36-mar2019.ff/forcefield.itp"
; Include topology for ions
#include "../Input_files/FORCEFIELDS/charmm36-mar2019.ff/ions.itp"
; Include water topology
#include "../Input_files/FORCEFIELDS/charmm36-mar2019.ff/tip3p.itp"

#ifdef POSRES_WATER
; Position restraint for each water oxygen
[ position_restraints ]
;  i funct       fcx        fcy        fcz
   1    1       1000       1000       1000
#endif

EOF

# Ahora, recorrer todos los archivos molecule.itp en el directorio actual y agregarlos:
for file in ../Output/molecule*.itp; do
    # Comprueba que existe al menos uno
    [ -e "$file" ] || continue
    echo "#include \"$file\"" >> ../Output/topol.top
done




    cat << 'EOF' >> ../Output/topol.top
[ system ]
; Name
Protein

[ molecules ]
; Compound        #mols
EOF

# Recorrer cada archivo itp y extraer el nombre de la molécula
for itp in ../Output/molecule*.itp; do
    if [[ -f "$itp" ]]; then
        molecule_name=$(awk '/^\[ moleculetype \]/ {getline; getline; print $1}' "$itp")
        
        if [[ -n "$molecule_name" ]]; then
            echo "Añadiendo molécula: $molecule_name"
            echo "$molecule_name     1" >> ../Output/topol.top
        fi
    else
        echo "Advertencia: No se encontró ningún archivo en ../Output/molecule*.itp"
    fi
done

gmx insert-molecules -f ../Output/complex.gro -ci ../Output/complex_2.gro -o ../Output/complex.gro -nmol 1
fi

gmx solvate -cp ../Output/complex.gro -cs spc216.gro -o ../Output/complex_water_box.gro  #-p ../Output/topol.top
gmx editconf -f ../Output/complex_water_box.gro -o ../Output/complex_water_box.pdb
grep -wc OW ../Output/complex_water_box.pdb 
water_count=$(grep -wc OW ../Output/complex_water_box.pdb)
echo -e "SOL     $water_count" >> ../Output/topol.top
# Eliminar archivos temporales
rm -f ../Output/#*


# Añadir iones
gmx grompp -f ../Input_files/MDP_FILES/ions.mdp -c ../Output/complex_water_box.gro -p ../Output/topol.top -o ../Output/ions.tpr -maxwarn 1
echo "SOL" | gmx genion -s ../Output/ions.tpr -o ../Output/complex_solv_ions.gro -p ../Output/topol.top -pname NA -nname CL -neutral

# Minimización de energía
gmx grompp -f ../Input_files/MDP_FILES/min.mdp -c ../Output/complex_solv_ions.gro -p ../Output/topol.top -o ../Output/min.tpr -maxwarn 2
gmx mdrun -v -deffnm ../Output/min -ntmpi 4

# Equilibración NVT
gmx grompp -f ../Input_files/MDP_FILES/eq_nvt.mdp -c ../Output/min.gro -p ../Output/topol.top -o ../Output/eq.tpr -r ../Output/min.gro -maxwarn 2
gmx mdrun -v -deffnm ../Output/eq -ntmpi 4

# Equilibración NPT
gmx grompp -f ../Input_files/MDP_FILES/eq_npt.mdp -c ../Output/eq.gro -p ../Output/topol.top -o ../Output/eq_2.tpr -r ../Output/eq.gro -maxwarn 2
gmx mdrun -v -deffnm ../Output/eq_2 -ntmpi 4

# Eliminar archivos temporales
rm -f ../Output/#* *mdp *pdb

# Crear archivo .tpr para producción
gmx grompp -f ../Input_files/MDP_FILES/prod.mdp -c ../Output/eq_2.gro -p ../Output/topol.top -o ../Output/prod.tpr -maxwarn 2

echo "Proceso completado exitosamente."
