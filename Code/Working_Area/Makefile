# Define variables
SCRIPT = main.sh

# Default rule
all: execute

# Rule to execute the script with an input file
execute:
	@if [ -z "$(pdb_complex)" ]; then \
		echo "Error: You must provide an input file. Use: make INPUT=file.txt"; \
		exit 1; \
	fi
	bash $(SCRIPT) $(pdb_complex)

# Help rule
help:
	@echo "Input file must be provided as follows:"
	@echo "  pdb_complex=file.pdb                          - Complex's input file"
	@echo ""
	@echo "Makefile Usage:"
	@echo "  make execute pdb_complex=file.pdb             - Runs $(SCRIPT)"
	@echo "  make help                                     - Displays this help message"

# Regla para limpiar (si fuera necesario)
clean:
	# No es necesario en este caso, pero se podría agregar para otros archivos
	rm -f *.pyc