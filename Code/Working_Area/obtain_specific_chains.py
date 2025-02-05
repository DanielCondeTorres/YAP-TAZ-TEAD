import MDAnalysis as mda
from MDAnalysis.coordinates.PDB import PDBWriter
import argparse

def filter_chains(input_pdb, output_pdb, chains):
    # Load the PDB file
    u = mda.Universe(input_pdb)
    
    # Select atoms from the specified chains
    selection = ' or '.join([f'segid {c}' for c in chains])
    atoms = u.select_atoms(selection)
    
    # Save the new PDB file
    with PDBWriter(output_pdb) as writer:
        writer.write(atoms)
    
    print(f"File saved: {output_pdb}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter specific chains from a PDB file")
    parser.add_argument("-pdb", required=True, help="Path to the input PDB file")
    parser.add_argument("-c", nargs='+', required=True, help="Chains to filter, separated by spaces (e.g., -c A B)")
    parser.add_argument("-o", required=True, help="Output file")
    args = parser.parse_args()
    
    input_pdb = args.pdb
    output_pdb = args.o
    filter_chains(input_pdb, output_pdb, args.c)


