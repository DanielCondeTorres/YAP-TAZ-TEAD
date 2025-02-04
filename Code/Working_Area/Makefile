# Define variables
SCRIPT = main.sh

# Default rule
all: execute

# Rule to execute the script with required input files
execute:
	@if [ -z "$(pdb_complex)" ]; then \
		echo "Error: Missing required input: pdb_complex"; \
		echo "Usage: make execute pdb_complex=file.pdb [forcefield=forcefield] [ligand=file.mol2] [config=config.txt]"; \
		exit 1; \
	fi
	# Execute the script with optional forcefield, ligand, and config inputs
	bash $(SCRIPT) $(forcefield) $(pdb_complex) $(ligand) $(config)


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