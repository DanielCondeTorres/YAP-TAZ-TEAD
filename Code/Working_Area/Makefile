# Define variables
SCRIPT = main.sh

# Default rule
all: execute

# Rule to execute the script with required input files
# make execute_complex_simultation pdb_complex=complex.pdb
execute_complex_simultation:
	@order=1; \
	echo "Passing order: $$order to the script"; \

	if [ -z "$(pdb_complex)" ]; then \
		echo "Error: Missing required input: pdb_complex"; \
		echo "Usage: make execute pdb_complex=file.pdb [forcefield=forcefield] [ligand=file.mol2] [config=config.txt]"; \
		exit 1; \
	fi
	# Execute the script with optional forcefield, ligand, and config inputs
	echo "Passing order: $$order to the script"
	bash $(SCRIPT) $(pdb_complex) $$order


execute_complex_formation:



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

help:
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
	rm ../Output/* 
	rm *mdp *pdb ../Output/.itp