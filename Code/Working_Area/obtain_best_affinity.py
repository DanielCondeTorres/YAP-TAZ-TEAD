import re
import argparse

def parse_vina_log(log_file):
    """
    Analiza el archivo log de AutoDock Vina para encontrar el modo con la menor afinidad.

    :param log_file: Ruta al archivo log de Vina.
    :return: Diccionario con información del modo de menor afinidad.
    """
    mode_data = []
    mode_pattern = re.compile(r'^\s*(\d+)\s+(-\d+\.\d+)\s+([\d\.]+)\s+([\d\.]+)')

    with open(log_file, 'r') as file:
        for line in file:
            match = mode_pattern.match(line)
            if match:
                mode = int(match.group(1))
                affinity = float(match.group(2))
                rmsd_lb = float(match.group(3))
                rmsd_ub = float(match.group(4))
                mode_data.append({
                    'mode': mode,
                    'affinity': affinity,
                    'rmsd_lb': rmsd_lb,
                    'rmsd_ub': rmsd_ub
                })

    if not mode_data:
        raise ValueError("No se encontraron datos de modos en el archivo log.")

    # Encontrar el modo con la menor afinidad
    best_mode = min(mode_data, key=lambda x: x['affinity'])
    return best_mode

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analiza un archivo de registro de Vina para encontrar el modo con la menor afinidad.")
    parser.add_argument('-log', '--log_file', type=str, required=True, help="Ruta al archivo de registro de Vina")
    args = parser.parse_args()

    log_file_path = args.log_file
    output_file_path = 'mode_output.txt'  # Archivo donde se guardará el resultado

    try:
        best_mode = parse_vina_log(log_file_path)
        mode = best_mode['mode']
        affinity = best_mode['affinity']
        print(f"El modo con la menor afinidad es el modo {mode} "
              f"con una afinidad de {affinity} kcal/mol.")
        # Escribir el valor de 'mode' en el archivo de salida
        with open(output_file_path, 'w') as f:
            f.write(str(mode))
    except Exception as e:
        print(f"Error al analizar el archivo log: {e}")
