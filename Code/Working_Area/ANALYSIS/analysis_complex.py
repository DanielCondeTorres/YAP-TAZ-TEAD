import numpy as np
import matplotlib.pyplot as plt
import argparse
import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
import mdtraj as md

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

def contact_map(u, cutoff: float = 12, frames_to_analyze: int = 1):
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








if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter specific chains from a PDB file")
    parser.add_argument("-pdb", required=True, help="Path to the input PDB file")
    parser.add_argument("-xtc", required=True, help="Trajectory file")
    args = parser.parse_args()

    # Load trajectory once
    u = load_trajectory(args.pdb, args.xtc)

    # Perform analysis on the loaded trajectory
    contact_map(u, cutoff=12, frames_to_analyze=10)
    # Calculate and visualize the average hydrogen bonds
