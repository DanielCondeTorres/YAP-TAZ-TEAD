import numpy as np
import matplotlib.pyplot as plt
import argparse
import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
import mdtraj as md
import matplotlib.patches as mpatches

def load_trajectory(pdb, xtc):
    """
    Load the system (structure and trajectory) into memory for faster subsequent analyses.

    Parameters:
    pdb (str): Path to the PDB file containing the protein structure.
    xtc (str): Path to the XTC file containing the trajectory.

    Returns:
    mda.Universe: Loaded MDAnalysis universe object containing both structure and trajectory.
    """
    u = mda.Universe(pdb, xtc)
    return u

def contact_map(u, cutoff: float = 5, frames_to_analyze: int = 1):
    """
    This function calculates and visualizes the contact map between two protein subunits 
    based on their Cα atoms from a trajectory.

    Parameters:
    u (mda.Universe): Loaded MDAnalysis universe object containing both structure and trajectory.
    cutoff (float): Distance threshold (in Angstroms) for considering a contact. Default is 12 Å.
    frames_to_analyze (int): Number of frames to analyze from the trajectory. Default is 1 frame.

    Returns:
    None: Displays a contact map as a heatmap.
    """

    # Print information about segments (molecules/chains)
    for seg in u.segments:
        print(f"Segment (Unit): {seg.segid}, Name: {seg.residues.resnames}")

    # Access the trajectory frames
    for ts in u.trajectory[:10]:  # Iterate over the first 10 frames
        print(f"Frame {ts.frame}: Time {ts.time} ps")

    # Select subunits (adjust the segment IDs based on your system)
    subunit1 = u.select_atoms("protein and (segid A or segid seg_0 or segid seg_0_Protein_chain_A or segid Protein_chain_A) and name CA")  # Cα atoms from subunit 1
    print('subunit1:', subunit1)
    subunit2 = u.select_atoms("protein and (segid B or segid seg_1 or segid seg_1_Protein_chain_B or segid Protein_chain_B) and name CA")  # Cα atoms from subunit 2

    # Get the residue names and renumber them starting from 1
    residues1 = [f"{res.resname}{i+1}" for i, res in enumerate(subunit1.residues)]
    residues2 = [f"{res.resname}{i+1}" for i, res in enumerate(subunit2.residues)]

    # Initialize the contact map (matrix)
    contact_map = np.zeros((len(subunit2.residues), len(subunit1.residues)))

    # Calculate contacts for each frame
    for ts in u.trajectory[:frames_to_analyze]:
        # Calculate the distance array between residues of subunit 1 and subunit 2
        distances = distance_array(subunit2.positions, subunit1.positions)
        
        # Define contacts as those with distance less than the cutoff
        contacts = distances < cutoff
        
        # Update the contact map (sum contacts in each frame)
        contact_map += contacts.astype(int)

    # Normalize the contact map
    contact_map /= frames_to_analyze

    # Create the contact map figure
    plt.figure(figsize=(18, 10))
    plt.pcolormesh(contact_map, cmap="Blues", edgecolors="black", linewidth=0.5)

    # Add color bar for contact frequency
    plt.colorbar(label="Contact Frequency")

    # Force all residues to be labeled on the axes
    plt.xticks(ticks=np.arange(len(residues1)), labels=residues1, rotation=90, fontsize=8)
    plt.yticks(ticks=np.arange(len(residues2)), labels=residues2, fontsize=8)
    plt.gca().set_xticks(np.arange(len(residues1)))  # Ensure every tick has a label
    plt.gca().set_yticks(np.arange(len(residues2)))

    # Labels and title
    plt.xlabel("Residue in Subunit 2")
    plt.ylabel("Residue in Subunit 1")
    plt.title("Contact Map Between Protein Subunits")
        # Adjust xticks and yticks spacing to make the plot more readable
    plt.xticks(ticks=np.arange(0, len(residues1), 2), labels=residues1[::2], rotation=90, fontsize=8)  # Show every 5th label
    plt.yticks(ticks=np.arange(0, len(residues2), 2), labels=residues2[::2], fontsize=8)  # Show every 5th label
    plt.gca().set_xticks(np.arange(len(residues1)))  # Ensure every tick has a label
    plt.gca().set_yticks(np.arange(len(residues2)))
    plt.tight_layout()
    # Show the plot
    plt.show()
    return residues1, residues2











def dssp_analysis(chain_id, trajectory_file, topology_file, residues_names):
    """
    Analyze the secondary structure of a protein chain over time using DSSP and visualize the dominant structure per residue.
    
    Parameters:
    chain_id (int): The chain ID of the protein to analyze.
    trajectory_file (str): Path to the trajectory file.
    topology_file (str): Path to the topology file.
    residues_names (list of str): List of residue names for labeling the x-axis.
    
    Returns:
    None: Displays a bar plot of the most frequent secondary structure per residue.
    """
    # Load trajectory and select only protein atoms from the given chain
    traj = md.load(trajectory_file, top=topology_file)
    protein_atoms_indices = traj.topology.select(f'chainid {chain_id}')
    traj_protein_slice = traj.atom_slice(protein_atoms_indices) 
    
    # Compute secondary structure using DSSP
    dssp = md.compute_dssp(traj_protein_slice)
    
    # Classify secondary structures
    helix_codes = {'H', 'G', 'I'}  # Helices
    beta_codes = {'E', 'B'}        # Beta sheets
    coil_codes = {'T', 'S', 'C'}   # Random coil
    
    num_residues = dssp.shape[1]  # Number of residues
    
    # Count frequencies of each structure per residue
    helix_freq = np.zeros(num_residues)
    beta_freq = np.zeros(num_residues)
    coil_freq = np.zeros(num_residues)
    
    for i in range(num_residues):
        res_dssp = dssp[:, i]  # Secondary structure of residue i over time
        helix_freq[i] = np.sum(np.isin(res_dssp, list(helix_codes))) / len(res_dssp)
        beta_freq[i] = np.sum(np.isin(res_dssp, list(beta_codes))) / len(res_dssp)
        coil_freq[i] = np.sum(np.isin(res_dssp, list(coil_codes))) / len(res_dssp)
    
    # Determine the most frequent conformation
    max_freq = np.maximum(np.maximum(helix_freq, beta_freq), coil_freq)
    colors = [
        'blue' if h == mf else 'orange' if b == mf else 'red'
        for h, b, mf in zip(helix_freq, beta_freq, max_freq)
    ]
    
    # Plot
    plt.figure(figsize=(16, 8))
    plt.bar(range(num_residues), max_freq * 100, color=colors)
    plt.xlabel('Residues', fontsize=20)
    plt.ylabel('Percentage of dominant conformation', fontsize=20)
    plt.title('Most frequent secondary structure per residue')
    plt.xticks(ticks=np.arange(0, len(residues_names), 3), labels=residues_names[::3], rotation=90, fontsize=10)
    plt.yticks(fontsize=10)
    plt.gca().set_xticks(np.arange(len(residues_names)))
    
    # Add legend
    legend_handles = [
        mpatches.Patch(color='blue', label='Helix'),
        mpatches.Patch(color='orange', label='Beta sheet'),
        mpatches.Patch(color='red', label='Random coil')
    ]
    plt.legend(handles=legend_handles, loc='upper right', fontsize=18)
    
    plt.show()








if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter specific chains from a PDB file")
    parser.add_argument("-pdb", required=True, help="Path to the input PDB file")
    parser.add_argument("-xtc", required=True, help="Trajectory file")
    args = parser.parse_args()

    # Load trajectory once
    u = load_trajectory(args.pdb, args.xtc)

    # Perform analysis on the loaded trajectory
    residues_1_names, residues_2_names = contact_map(u, cutoff=8, frames_to_analyze=10)
    dssp_analysis(0, args.xtc, args.pdb,residues_1_names)
    dssp_analysis(1, args.xtc, args.pdb,residues_2_names)
    # Calculate and visualize the average hydrogen bonds