import argparse

def comentar_includes(archivo_topologia):
    # Leer el contenido del archivo original
    with open(archivo_topologia, 'r') as archivo:
        lineas = archivo.readlines()

    nueva_topologia = []
    omitir_restante = False

    for i, linea in enumerate(lineas):
        # Si encontramos el bloque de restricciones de posición, lo agregamos y dejamos de procesar el archivo
        if "; Include Position restraint file" in linea:
            nueva_topologia.append(linea)
            nueva_topologia.append(lineas[i + 1])  # #ifdef POSRES
            nueva_topologia.append(lineas[i + 2])  # #include "posre.itp"
            nueva_topologia.append(lineas[i + 3])  # #endif
            break  # Salimos del bucle

        # Comentar líneas que empiezan con #include, excepto la de posre.itp
        if linea.strip().startswith("#include") or linea.strip().startswith("# include") and 'posre.itp' not in linea:
            nueva_topologia.append(f";{linea}")  # Comentamos la línea
        else:
            nueva_topologia.append(linea)

    # Guardar en ligand.itp
    with open("ligand.itp", 'w') as archivo:
        archivo.writelines(nueva_topologia)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Comentar includes en un archivo .top excepto posre.itp y guardar en ligand.itp.")
    parser.add_argument("-top", required=True, help="Archivo de topología de entrada.")
    
    args = parser.parse_args()
    
    comentar_includes(args.top)

