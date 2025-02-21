from Bio import PDB
import argparse
import string

def get_next_chain_id(chain_ids_in_use):
    """Devuelve el siguiente identificador de cadena disponible."""
    all_chain_ids = string.ascii_uppercase + string.ascii_lowercase
    for chain_id in all_chain_ids:
        if chain_id not in chain_ids_in_use:
            return chain_id
    raise ValueError("No hay identificadores de cadena disponibles.")

def rename_chain(structure, new_chain_id):
    """Renombra el identificador de cadena de una estructura."""
    for model in structure:
        for chain in model:
            chain.id = new_chain_id

def merge_structures(receptor_pdb, ligand_pdb, output_pdb):
    parser = PDB.PDBParser(QUIET=True)
    receptor_structure = parser.get_structure('receptor', receptor_pdb)
    ligand_structure = parser.get_structure('ligand', ligand_pdb)

    # Obtener los identificadores de cadena existentes en el receptor
    receptor_chain_ids = {chain.id for model in receptor_structure for chain in model}

    # Asignar un nuevo identificador de cadena al ligando
    new_chain_id = get_next_chain_id(receptor_chain_ids)
    rename_chain(ligand_structure, new_chain_id)

    # Combinar las estructuras
    combined_structure = receptor_structure
    for model in ligand_structure:
        for chain in model:
            combined_structure[0].add(chain)

    # Guardar la estructura combinada
    io = PDB.PDBIO()
    io.set_structure(combined_structure)
    io.save(output_pdb)

    print(f"Estructura combinada guardada en {output_pdb}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Combina dos archivos PDB asegurando identificadores de cadena consistentes.")
    parser.add_argument("-receptor", required=True, help="Archivo PDB del receptor")
    parser.add_argument("-ligand", required=True, help="Archivo PDB del ligando")
    parser.add_argument("-output", default="estructura_combinada.pdb", help="Nombre del archivo PDB combinado")

    args = parser.parse_args()
    merge_structures(args.receptor, args.ligand, args.output)

