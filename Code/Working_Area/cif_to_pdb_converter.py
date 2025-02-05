from Bio import PDB
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBIO import PDBIO
import argparse

def convert_cif_to_pdb(input_cif, output_pdb):
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("structure", input_cif)
    
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb)
    print(f"Conversion completed: {output_pdb}")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter specific chains from a PDB file")
    parser.add_argument("-cif", required=True, help="Path to the input CIF file")
    parser.add_argument("-o", required=True, help="Output file in pdb format")
    args = parser.parse_args()
    
    input_cif = args.cif
    output_pdb = args.o
    convert_cif_to_pdb(input_cif, output_pdb)


