import MDAnalysis as mda
import numpy as np
import argparse

def calculate_box_dimensions(pdb_file):
    u = mda.Universe(pdb_file)
    protein = u.select_atoms('protein')
    
    # Calcular el centro de masa de la proteína
    center_of_mass = protein.center_of_mass()
    
    # Obtener las coordenadas de todos los átomos de la proteína
    coords = protein.positions
    
    # Calcular las dimensiones de la caja
    min_coords = np.min(coords, axis=0)
    max_coords = np.max(coords, axis=0)
    box_size = max_coords - min_coords
    
    return center_of_mass, box_size

def generate_config(receptor, ligand, output, center, size):
    config_content = f"""receptor = {receptor}
ligand = {ligand}
out = {output}

center_x = {center[0]:.3f}
center_y = {center[1]:.3f}
center_z = {center[2]:.3f}

size_x = {size[0]:.3f}
size_y = {size[1]:.3f}
size_z = {size[2]:.3f}

"""
    with open('config.txt', 'w') as config_file:
        config_file.write(config_content)
    print("Archivo 'config.txt' generado con éxito.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Genera un archivo config.txt para AutoDock Vina.")
    parser.add_argument("-r", "--receptor", required=True, help="Archivo PDBQT del receptor.")
    parser.add_argument("-l", "--ligand", required=True, help="Archivo PDBQT del ligando.")
    parser.add_argument("-p", "--pdb", required=True, help="Archivo PDB para calcular dimensiones de la caja.")
    parser.add_argument("-o", "--output", default="all.pdbqt", help="Nombre del archivo de salida.")
    
    args = parser.parse_args()
    
    center, size = calculate_box_dimensions(args.pdb)
    generate_config(args.receptor, args.ligand, args.output, center, size)

