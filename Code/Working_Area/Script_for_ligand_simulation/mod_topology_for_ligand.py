# Nombres de los archivos
archivo_original = 'topol.top'
archivo_modificado = 'topol_modificado.top'

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

