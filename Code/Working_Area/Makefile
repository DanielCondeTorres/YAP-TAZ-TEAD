# Define variables
SCRIPT = main.sh

# Default rule
all: execute

# Rule to execute the script with required input files
# make execute_complex_simultation pdb_complex=complex.pdb
execute_complex_simultation:
	@if [ -z "$(pdb_complex_1)" ]; then \
		echo "Error: Missing required input: pdb_complex_1"; \
		echo "Usage: make execute_complex_simultation pdb_complex_1=file1.pdb pdb_complex_2=file2.pdb"; \
		exit 1; \
	fi
	# Ejecutar el script pasando los parámetros
	bash $(SCRIPT) $(pdb_complex_1) $(pdb_complex_2)


# make cif_to_pdb cif_file=file.cif output_file=file.pdb
cif_to_pdb:
# Check if both files are provided
	@if [ ! -z "$(cif_file)" ] && [ ! -z "$(output_file)" ]; then \
		python cif_to_pdb_converter.py -cif $(cif_file) -o $(output_file); \
	else \
		echo "Error: Missing required input: cif_file or output_file"; \
		echo "Usage: make execute cif_file=file.cif output_file=output.pdb"; \
		exit 1; \
	fi
# Help rule
production:
	sbatch production.sh

# Rule usage:  make execute_complex_ligand_formation PDB_INPUT=../HUMAN_3KYS/YAP_TEAD_YAP_FROM_PDB/receptor_docking.pdb  RECEPTOR_CHAIN=B LIGAND=verteporfin_gmx22.pdb PATH_TO=/Users/danielcondetorres/Applications/mgltools_1.5.7_MacOS-X/MGLToolsPckgs/AutoDockTools/Utilities24 PATH_TO_pythonsh=/Users/danielcondetorres/Applications/mgltools_1.5.7_MacOS-X/bin/pythonsh
execute_complex_ligand_formation:
# Convert .pdb complex last frame to a complex without water.
        echo "Protein" | gmx trjconv -f $(PDB_INPUT) -o complex_protein.pdb -s $(PDB_INPUT)
# Obtain the CHAIN that it is going to be used to the docking
        python obtain_specific_chains.py -pdb complex_protein.pdb -c $(RECEPTOR_CHAIN) -o receptor.pdb
		cp receptor.pdb ../Output/
# Create pdbqt files
        $(PATH_TO_pythonsh) $(PATH_TO)/prepare_receptor4.py -r receptor.pdb -o receptor.pdbqt
        $(PATH_TO_pythonsh) $(PATH_TO)/prepare_ligand4.py -l $(LIGAND) -o ligand.pdbqt
# Create config.txt file
        python Script_for_ligand_protein_formation/generate_vina_config.py -r receptor.pdbqt -l ligand.pdbqt -p $(PDB_INPUT) -o config.txt
        vina --config config.txt --cpu 1 --num_modes 15 > log.txt
        python Script_for_ligand_protein_formation/obtain_best_affinity.py -log log.txt
		cp log.txt ../Output/
# Extract the best conformation
#run_docking.sh
		bash Script_for_ligand_protein_formation/run_docking.sh $(LIGAND) $(PDB_INPUT)
		cp complex_ligand.pdb ../Output/
		cp posre.itp ../Output/

# Rule usage:  make execute_complex_ligand_simulation PDB_INPUT_COMPLEX_LIGAND=../Output/complex_ligand.pdb
execute_complex_ligand_simulation:
	python Script_for_ligand_simulation/mod_topology_for_ligand.py -top $(TOPOLOGY) -o $(OUTPUT_TOPOLOGY_PATH)
	bash main_ligand.sh $(PDB_INPUT_COMPLEX_LIGAND)

# Help rule
help:
	@echo "Required input files must be provided as follows:"
	@echo "  pdb_complex=file.pdb                          - Complex's input file (required)"
	@echo "  forcefield=forcefield                         - Forcefield to use (required)"
	@echo "  config=config.txt                             - Configuration file (optional)"
	@echo "  ligand=file.mol2                              - Ligand input file (optional)"
	@echo ""
	@echo "Makefile Usage:"
	@echo "  make execute pdb_complex=file.pdb config=config.txt [ligand=file.mol2]  - Runs $(SCRIPT) to generate production.tpr file"
	@echo "  make production                                                         - Runs production.sh to submit the simulation to a supercomputer using an SLURM job scheduler"
	@echo "    (NOTE: Modify production.sh according to your supercomputer's job scheduler)"
	@echo "    (The current script is set up for the following SLURM settings):"
	@echo "  make help                                                               - Displays this help message"


# Regla para limpiar (si fuera necesario)
clean:
	# No es necesario en este caso, pero se podría agregar para otros archivos
	rm *mdp *pdb *txt 

cleanOutput:
	# No es necesario en este caso, pero se podría agregar para otros archivos
	rm ../Output/* 
	rm *mdp *pdb ../Output/.itp