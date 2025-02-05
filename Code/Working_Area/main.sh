#!/bin/bash


echo "All inputs received: $@"
# Recibe los argumentos del Makefile
COMPLEX_PDB_FILE="$1"  # Archivo PDB del complejo
Order="$2"  # Orden de ejecución (1 para ejecutar, otro valor para no hacer nada)
echo "Executing with order: $Order"
# Verifica si Order es igual a 1
if [ "$Order" -eq 1 ]; then
    # Crear directorio de salida
    mkdir -p ../Output/

    # Convertir PDB a GRO y generar topología
    GMXLIB=../Input_files/FORCEFIELDS gmx pdb2gmx -f "$COMPLEX_PDB_FILE" -o ../Output/complex.gro -ignh -p ../Output/topol.top -i ../Output/posre.itp

    # Crear caja de agua
    gmx editconf -f ../Output/complex.gro -o ../Output/complex_water_box.gro -c -d 1.2 -bt cubic
    gmx solvate -cp ../Output/complex_water_box.gro -cs spc216.gro -o ../Output/complex_water_box.gro -p ../Output/topol.top

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
else
    echo "Order no es igual a 1. No se ejecutará el proceso."
fi
