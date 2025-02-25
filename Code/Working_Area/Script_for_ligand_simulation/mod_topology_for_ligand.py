import argparse

def leer_ligand_topology(archivo_ligand):
    with open(archivo_ligand, 'r') as archivo:
        lineas = archivo.readlines()
    
    parametros_ligand = []
    for i, linea in enumerate(lineas):
        if '; Include ligand specific parameters' in linea:
            parametros_ligand.append(linea)
            if i + 1 < len(lineas):
                parametros_ligand.append(lineas[i + 1])
            break
    return parametros_ligand

def modificar_topologia(archivo_original, archivo_modificado, archivo_ligand):
    # Leer los parámetros del archivo de topología del ligando
    parametros_ligand = leer_ligand_topology(archivo_ligand)
    
    # Leer el contenido del archivo original
    with open(archivo_original, 'r') as archivo:
        lineas = archivo.readlines()

    nueva_topologia = []
    ligand_included = False
    en_seccion_moleculas = False

    for linea in lineas:
        # Insertar los parámetros del ligando justo después del primer #include encontrado
        if not ligand_included and linea.lstrip().startswith("#include"):
            nueva_topologia.append(linea)
            nueva_topologia.append('\n; Include ligand specific parameters\n')
            nueva_topologia.extend(parametros_ligand)
            nueva_topologia.append('\n; Include ligand topology\n#include "ligand.itp"\n')
            ligand_included = True
            continue

        # Procesar la sección [ molecules ]: omitir líneas que contengan SOL o NA
        if '[ molecules ]' in linea:
            en_seccion_moleculas = True
            nueva_topologia.append(linea)
            continue

        if en_seccion_moleculas:
            if 'SOL' in linea or 'NA' in linea:
                continue
            if linea.strip() == '':
                en_seccion_moleculas = False

        nueva_topologia.append(linea)

    nueva_topologia.append('Other           1\n')

    with open(archivo_modificado, 'w') as archivo:
        archivo.writelines(nueva_topologia)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Modificar archivo de topología de GROMACS.")
    parser.add_argument("-top", required=True, help="Archivo de topología de entrada.")
    parser.add_argument("-o", required=True, help="Archivo de topología de salida.")
    parser.add_argument("-ligand_topology", required=True, help="Archivo de parámetros del ligando.")
    
    args = parser.parse_args()
    
    modificar_topologia(args.top, args.o, args.ligand_topology)

