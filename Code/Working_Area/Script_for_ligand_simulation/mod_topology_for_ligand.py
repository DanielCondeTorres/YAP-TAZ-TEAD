def modificar_topologia(archivo_original, archivo_modificado):
    # Leer el contenido del archivo original
    with open(archivo_original, 'r') as archivo:
        lineas = archivo.readlines()

    # Inicializar variables
    nueva_topologia = []
    incluido_ligand = False
    en_seccion_moleculas = False

    # Procesar cada línea del archivo original
    for linea in lineas:
        # Incluir ligand.itp antes de las secciones [ system ] y [ molecules ]
        if not incluido_ligand and ('[ system ]' in linea or '[ molecules ]' in linea):
            nueva_topologia.append('\n; Include ligand topology\n#include "ligand.itp"\n')
            incluido_ligand = True

        # Verificar si estamos en la sección [ molecules ]
        if '[ molecules ]' in linea:
            en_seccion_moleculas = True
            nueva_topologia.append(linea)
            continue

        # Omitir las líneas que contienen SOL o NA en la sección [ molecules ]
        if en_seccion_moleculas:
            if 'SOL' in linea or 'NA' in linea:
                continue
            # Si encontramos una línea vacía, significa el fin de la sección [ molecules ]
            if linea.strip() == '':
                en_seccion_moleculas = False

        # Agregar la línea actual a la nueva topología
        nueva_topologia.append(linea)

    # Añadir 'Other           1' al final del archivo
    nueva_topologia.append('Other           1\n')

    # Escribir el contenido modificado en un nuevo archivo
    with open(archivo_modificado, 'w') as archivo:
        archivo.writelines(nueva_topologia)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Modificar archivo de topología de GROMACS.")
    parser.add_argument("-top", required=True, help="Archivo de topología de entrada.")
    parser.add_argument("-o", required=True, help="Archivo de topología de salida.")
    
    args = parser.parse_args()
    
    modificar_topologia(args.top, args.o)


