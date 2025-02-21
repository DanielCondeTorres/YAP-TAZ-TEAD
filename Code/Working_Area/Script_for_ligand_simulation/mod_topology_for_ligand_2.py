import os
import argparse

def modificar_topologia(archivo_entrada, archivo_salida):
    # Verificar si el archivo de entrada existe
    if not os.path.isfile(archivo_entrada):
        print(f"El archivo '{archivo_entrada}' no existe.")
        return

    # Definir las rutas a reemplazar
    ruta_antigua = '#include "../../Input_files/FORCEFIELDS/charmm36-mar2019.ff'
    ruta_nueva = '#include "./charmm36.ff'

    try:
        # Leer el contenido del archivo de entrada
        with open(archivo_entrada, 'r') as infile:
            contenido = infile.read()

        # Reemplazar todas las ocurrencias de la ruta antigua por la nueva
        contenido_modificado = contenido.replace(ruta_antigua, ruta_nueva)

        # Escribir el contenido modificado en el archivo de salida
        with open(archivo_salida, 'w') as outfile:
            outfile.write(contenido_modificado)

        print(f"El archivo ha sido procesado y guardado como '{archivo_salida}'.")

    except Exception as e:
        print(f"Ocurrió un error al procesar el archivo: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Modificar rutas de inclusión en un archivo de topología de GROMACS.")
    parser.add_argument("-top", required=True, help="Archivo de topología de entrada (e.g., topol.top)")
    parser.add_argument("-o", required=True, help="Archivo de topología de salida con las rutas modificadas")

    args = parser.parse_args()

    modificar_topologia(args.top, args.o)

