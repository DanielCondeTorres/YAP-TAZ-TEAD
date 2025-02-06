import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
import numpy as np
import matplotlib.pyplot as plt
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
import argparse



def hydrogen_bond_matrix(tpr, xtc):
    """
    This function calculates and visualizes the hydrogen bond matrix between two protein chains from a molecular dynamics trajectory.

    Parameters:
    - tpr (str): The TPR file of the simulation (topology).
    - xtc (str): The XTC file of the simulation (trajectory).

    The function performs the following tasks:
    1. Loads the molecular dynamics trajectory.
    2. Defines two protein chains (by segment identifiers).
    3. Extracts residue names for both protein chains.
    4. Calculates hydrogen bonds using the HBA analysis.
    5. Builds a matrix where each entry represents the number of hydrogen bonds between residues of different chains.
    6. Normalizes the matrix based on the total number of frames.
    7. Visualizes the hydrogen bond matrix as a heatmap.
    """

    # Load the trajectory
    u = mda.Universe(tpr, xtc)

    # Define proteins (example using chain identifiers)
    prot1 = u.select_atoms('segid A or segid seg_0 or segid seg_0_Protein_chain_A or segid Protein_chain_A')
    prot2 = u.select_atoms('segid B or segid seg_1 or segid seg_1_Protein_chain_B or segid Protein_chain_B')
    
    # Extract residue names for each protein chain
    residues_names1 = [f"{res.resname}{i+1}" for i, res in enumerate(prot1.residues)]
    residues_names2 = [f"{res.resname}{i+1}" for i, res in enumerate(prot2.residues)]

    # Get the unique and sorted residue list for each protein chain
    residues1 = sorted(set(atom.resid for atom in prot1))
    residues2 = sorted(set(atom.resid for atom in prot2))

    # Initialize the hydrogen bond matrix (rows: residues of prot1, columns: residues of prot2)
    matrix = np.zeros((len(residues1), len(residues2)))

    # Run the hydrogen bond analysis (considering all protein atoms)
    hbonds = HBA(universe=u, donors_sel="protein", acceptors_sel="protein")
    hbonds.run()

    # Iterate over the results using hbonds.results.hbonds (in MDAnalysis 2.x)
    for hbond in hbonds.results.hbonds:
        # Each 'hbond' is an array with [frame, donor_index, hydrogen_index, acceptor_index, distance, angle]
        _, donor_idx, _, acceptor_idx, _, _ = hbond
        
        # Convert to integers
        donor_idx = int(donor_idx)
        acceptor_idx = int(acceptor_idx)
        
        # Get the corresponding atoms
        donor_atom = u.atoms[donor_idx]
        acceptor_atom = u.atoms[acceptor_idx]

        # Consider hydrogen bonds that connect two different protein chains
        if (donor_atom in prot1 and acceptor_atom in prot2):
            r1 = donor_atom.resid
            r2 = acceptor_atom.resid
        elif (donor_atom in prot2 and acceptor_atom in prot1):
            r1 = acceptor_atom.resid
            r2 = donor_atom.resid
        else:
            continue

        # Try to find the residue indices for each protein
        try:
            i = residues1.index(r1)
            j = residues2.index(r2)
        except ValueError:
            continue

        # Increment the corresponding matrix element
        matrix[i, j] += 1

    # (Optional) Normalize the matrix by the total number of frames to get the occupancy
    matrix /= 2 * (u.trajectory.n_frames + 1)
    print(f"{u.trajectory.n_frames} FRAMES")
    print('MATRIX: ', matrix)

    # Get indices where the matrix is non-zero
    indices = np.nonzero(matrix)

    # Iterate over and print each non-zero element with its indices
    for i, j in zip(indices[0], indices[1]):
        print(f"matrix[{i}, {j}] = {matrix[i, j]}")

    # Visualize the matrix as a heatmap
    plt.figure(figsize=(18, 10))

    # Plot the heatmap
    plt.pcolormesh(matrix.T, cmap="Blues", edgecolors="black", linewidth=0.5)
    plt.colorbar(label="Hydrogen Bond Frequency")
    plt.xlabel('Residue (Protein B)')
    plt.ylabel('Residue (Protein A)')
    plt.title('Hydrogen Bond Matrix Between Proteins')

    # Customize tick labels for residues
    plt.xticks(ticks=np.arange(0, len(residues_names1), 2), labels=residues_names1[::2], rotation=90, fontsize=8)  # Show every 2nd label
    plt.yticks(ticks=np.arange(0, len(residues_names2), 2), labels=residues_names2[::2], fontsize=8)  # Show every 2nd label
    plt.gca().set_xticks(np.arange(len(residues_names1)))  # Ensure every tick has a label
    plt.gca().set_yticks(np.arange(len(residues_names2)))
    plt.tight_layout()
    plt.show()










if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter specific chains from a PDB file")
    parser.add_argument("-tpr", required=True, help="Path to the input tpr file")
    parser.add_argument("-xtc", required=True, help="Trajectory file")
    args = parser.parse_args()
    hydrogen_bond_matrix(args.tpr, args.xtc)

