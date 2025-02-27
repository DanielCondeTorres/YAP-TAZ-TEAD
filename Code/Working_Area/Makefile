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
	# Execute the script passing the parameters
	bash $(SCRIPT) $(pdb_complex_1) $(pdb_complex_2)

# Rule to convert CIF to PDB
# make cif_to_pdb cif_file=file.cif output_file=file.pdb
cif_to_pdb:
	# Check if both files are provided
	@if [ ! -z "$(cif_file)" ] && [ ! -z "$(output_file)" ]; then \
		python cif_to_pdb_converter.py -cif $(cif_file) -o $(output_file); \
	else \
		echo "Error: Missing required input: cif_file or output_file"; \
		echo "Usage: make cif_to_pdb cif_file=file.cif output_file=output.pdb"; \
		exit 1; \
	fi

# Help rule to provide general information about the Makefile commands
production:
	sbatch production.sh

# Rule for executing complex ligand formation from PDB file
# Usage example:
# make execute_complex_ligand_formation PDB_INPUT=../HUMAN_3KYS/YAP_TEAD_YAP_FROM_PDB/receptor_docking.pdb  RECEPTOR_CHAIN=B LIGAND=verteporfin_gmx22.pdb PATH_TO=/Users/danielcondetorres/Applications/mgltools_1.5.7_MacOS-X/MGLToolsPckgs/AutoDockTools/Utilities24 PATH_TO_pythonsh=/Users/danielcondetorres/Applications/mgltools_1.5.7_MacOS-X/bin/pythonsh OUTPUT_DIR=./output/
execute_complex_ligand_formation:
	# Convert the last frame of the complex PDB to a complex without water
	echo "Protein" | gmx trjconv -f $(PDB_INPUT) -o complex_protein.pdb -s $(PDB_INPUT)
	# Obtain the chain for docking
	python obtain_specific_chains.py -pdb complex_protein.pdb -c $(RECEPTOR_CHAIN) -o receptor.pdb
	cp receptor.pdb $(OUTPUT_DIR)
	cp $(LIGAND) .
	# Prepare receptor and ligand in pdbqt format
	$(PATH_TO_pythonsh) $(PATH_TO)/prepare_receptor4.py -r receptor.pdb -o receptor.pdbqt
	$(PATH_TO_pythonsh) $(PATH_TO)/prepare_ligand4.py -l $(LIGAND) -o ligand.pdbqt
	# Create config.txt file for docking
	python Script_for_ligand_protein_formation/generate_vina_config.py -r receptor.pdbqt -l ligand.pdbqt -p $(PDB_INPUT) -o config.txt
	# Run Vina for docking
	vina --config config.txt --cpu 1 --num_modes 50 --exhaustiveness=24 > log.txt
	# Process the log file to get best affinity
	python Script_for_ligand_protein_formation/obtain_best_affinity.py -log log.txt
	cp log.txt $(OUTPUT_DIR)
	# Extract the best conformation and run docking
	bash Script_for_ligand_protein_formation/run_docking.sh $(LIGAND) $(PDB_INPUT)
	cp posre.itp $(OUTPUT_DIR)
	head -n 5 "$(PDB_INPUT)" > complex_ligand.pdb
	cat estructura_combinada.pdb >> complex_ligand.pdb
	cp complex_ligand.pdb $(OUTPUT_DIR)

# Rule for executing complex ligand simulation
# Usage example:
# make execute_complex_ligand_simulation TOPOLOGY=../HUMAN_3KYS/YAP_TEAD_YAP_FROM_ALPHA_FOLD/topol.top OUTPUT_DIR=../HUMAN_3KYS/YAP_TEAD_YAP_FROM_ALPHA_FOLD/ LIGAND_TOP=../Input_files/LIGANDS/For_YAP_TAZ/verteporfin_gromacs/verteporfin_gmx.top    PDB_INPUT_COMPLEX_LIGAND=../HUMAN_3KYS/YAP_TEAD_YAP_FROM_ALPHA_FOLD/complex_ligand.pdb FF=../Input_files/LIGANDS/For_YAP_TAZ/verteporfin_gromacs/charmm36.ff
execute_complex_ligand_simulation:
	# Modify topology for ligand
	python Script_for_ligand_simulation/mod_topology_for_ligand.py -top $(TOPOLOGY) -o topol_mod.top -ligand_topology $(LIGAND_TOP)
	# Modify the ligand topology further
	python Script_for_ligand_simulation/mod_topology_for_ligand_2.py -top topol_mod.top -o topol_mod2.top
	cp topol_mod2.top $(OUTPUT_DIR)
	cp -r $(FF) $(OUTPUT_DIR)
	# Convert ligand topology to itp format
	python Script_for_ligand_simulation/topol_to_itp.py -top $(LIGAND_TOP)
	cp ligand.itp $(OUTPUT_DIR)
	# Run the main ligand simulation
	bash Script_for_ligand_simulation/main_ligand.sh $(PDB_INPUT_COMPLEX_LIGAND) $(OUTPUT_DIR)

# Help rule to display usage information
help:
	@echo "Required input files must be provided as follows:"
	@echo "  pdb_complex_1=file.pdb                    - First complex input PDB file (required)"
	@echo "  pdb_complex_2=file.pdb                    - Second complex input PDB file (required)"
	@echo "  cif_file=file.cif                        - Input CIF file (required for cif_to_pdb rule)"
	@echo "  output_file=file.pdb                     - Output PDB file (required for cif_to_pdb rule)"
	@echo "  PDB_INPUT=file.pdb                       - PDB input file for ligand formation (required)"
	@echo "  RECEPTOR_CHAIN=A                          - Receptor chain to be used in ligand formation (required)"
	@echo "  LIGAND=file.pdb                          - Ligand file (required for ligand formation)"
	@echo "  PATH_TO=/path/to/AutoDockTools           - Path to AutoDockTools utilities (required)"
	@echo "  PATH_TO_pythonsh=/path/to/pythonsh       - Path to pythonsh (required)"
	@echo "  OUTPUT_DIR=/path/to/output               - Directory for output files (required)"
	@echo "  TOPOLOGY=file.top                        - Topology file for ligand simulation (required)"
	@echo "  LIGAND_TOP=file.top                      - Ligand topology file for simulation (required)"
	@echo "  FF=file.ff                               - Forcefield files directory (required)"
	@echo "  PDB_INPUT_COMPLEX_LIGAND=file.pdb        - PDB complex file for ligand simulation (required)"
	@echo ""
	@echo "Makefile Usage:"
	@echo "  make execute_complex_simultation pdb_complex_1=file1.pdb pdb_complex_2=file2.pdb   - Runs $(SCRIPT) to simulate the complex"
	@echo "  make cif_to_pdb cif_file=file.cif output_file=output.pdb                            - Converts CIF to PDB format"
	@echo "  make execute_complex_ligand_formation PDB_INPUT=complex.pdb RECEPTOR_CHAIN=A LIGAND=ligand.pdb PATH_TO=/path/to/tools PATH_TO_pythonsh=/path/to/pythonsh OUTPUT_DIR=output/   - Prepares receptor and ligand for docking"
	@echo "  make execute_complex_ligand_simulation TOPOLOGY=topol.top OUTPUT_DIR=output/ LIGAND_TOP=ligand.top PDB_INPUT_COMPLEX_LIGAND=complex_ligand.pdb FF=forcefield/   - Runs the ligand simulation"
	@echo "  make production                                                                - Submits the simulation job to a supercomputer using SLURM job scheduler"
	@echo "    (NOTE: Modify production.sh according to your supercomputer's job scheduler)"
	@echo "  make help                                                                      - Displays this help message"

# Rule to clean up generated files (if necessary)
clean:
	# Remove unnecessary files
	rm *mdp *pdb *txt *itp *pdbqt *top  

# Rule to clean up output directory
cleanOutput:
	# Remove output files
	rm ../Output/*
	rm *mdp *pdb ../Output/.itp
