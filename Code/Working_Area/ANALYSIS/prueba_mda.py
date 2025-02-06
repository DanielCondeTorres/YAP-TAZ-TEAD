import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
import numpy as np
import matplotlib.pyplot as plt
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA

# Cargar la trayectoria
u = mda.Universe('prod.tpr', 'out.xtc')

# Definir proteínas (ejemplo por cadenas)
prot1 = u.select_atoms('segid A or segid seg_0 or segid seg_0_Protein_chain_A or segid Protein_chain_A')
prot2 = u.select_atoms('segid B or segid seg_1 or segid seg_1_Protein_chain_B or segid Protein_chain_B')
residues_names1 = [f"{res.resname}{i+1}" for i, res in enumerate(prot1.residues)]
residues_names2 = [f"{res.resname}{i+1}" for i, res in enumerate(prot2.residues)]

# Obtener la lista de residuos de cada proteína (asegúrate de que sean únicos y ordenados)
residues1 = sorted(set(atom.resid for atom in prot1))
residues2 = sorted(set(atom.resid for atom in prot2))

# Inicializar la matriz (filas: residuos de prot1, columnas: residuos de prot2)
matriz = np.zeros((len(residues1), len(residues2)))

# Ejecutar el análisis de enlaces de hidrógeno (se consideran todos los átomos de proteína)
hbonds = HBA(universe=u, donors_sel="protein", acceptors_sel="protein")
hbonds.run()

# Iterar sobre los resultados usando hbonds.results.hbonds (en MDAnalysis 2.x)
for hbond in hbonds.results.hbonds:
    # Cada 'hbond' es un array con [frame, donor_index, hydrogen_index, acceptor_index, distance, angle]
    _, donor_idx, _, acceptor_idx, _, _ = hbond
    # Convertir a enteros
    donor_idx = int(donor_idx)
    acceptor_idx = int(acceptor_idx)
    
    # Obtener los átomos correspondientes
    donor_atom = u.atoms[donor_idx]
    acceptor_atom = u.atoms[acceptor_idx]

    # Considerar enlaces que conectan dos proteínas distintas
    if (donor_atom in prot1 and acceptor_atom in prot2):
        r1 = donor_atom.resid
        r2 = acceptor_atom.resid
    elif (donor_atom in prot2 and acceptor_atom in prot1):
        r1 = acceptor_atom.resid
        r2 = donor_atom.resid
    else:
        continue

    try:
        i = residues1.index(r1)
        j = residues2.index(r2)
    except ValueError:
        continue

    matriz[i, j] += 1;print(matriz[48, 27])
# (Opcional) Normalizar la matriz por el número total de frames para obtener la ocupación
matriz /= 2*(u.trajectory.n_frames+1)
print(u.trajectory.n_frames,'FRAMES')
print('MATRIZ: ',matriz)
# Obtener los índices donde la matriz no es cero
indices = np.nonzero(matriz)

# Recorrer e imprimir cada elemento no cero junto con sus índices
for i, j in zip(indices[0], indices[1]):
    print(f"matriz[{i}, {j}] = {matriz[i, j]}")
# Representación visual de la matriz usando un mapa de calor



plt.figure(figsize=(18, 10))


#im = plt.imshow(matriz, cmap='viridis', interpolation='nearest')
#plt.colorbar(im, label='Fracción de tiempo con enlace')
plt.pcolormesh(matriz.T, cmap="Blues", edgecolors="black", linewidth=0.5)
plt.colorbar(label="Hydrogen Bond Frequency")
plt.xlabel('Residue (Protein B)')
plt.ylabel('Residue (Protein A)')
plt.title('Matriz de Enlaces de Hidrógeno entre Proteínas')

plt.xticks(ticks=np.arange(0, len(residues_names1), 2), labels=residues_names1[::2], rotation=90, fontsize=8)  # Show every 5th label
plt.yticks(ticks=np.arange(0, len(residues_names2), 2), labels=residues_names2[::2], fontsize=8)  # Show every 5th label
plt.gca().set_xticks(np.arange(len(residues_names1)))  # Ensure every tick has a label
plt.gca().set_yticks(np.arange(len(residues_names2)))
plt.tight_layout()
plt.show()
