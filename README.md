# YAP/TAZ-TEAD
This repository provides code that utilizes GROMACS to prepare simulations of protein-protein complexes in a water box. Additionally, it allows for docking of ligands to these proteins and prepares structures that can be simulated using molecular dynamics (MD).

## Expected outcomes

-	Detailed understanding of the protein-protein complex and evaluation of potential next steps. Molecular dynamics simulations evaluating the effectiveness of proposed drugs.
-	An automated GROMACS protocol for protein-ligand simulations.
- An integrated workflow combining docking and molecular dynamics simulations.

![taz_tead](https://github.com/user-attachments/assets/28c65242-0d69-45f2-8466-75fda9355c3d)


[![GitHub forks](https://img.shields.io/github/forks/DanielCondeTorres/YAP-TAZ-TEAD)](https://github.com/DanielCondeTorres/YAP-TAZ-TEAD/network)
[![GitHub issues](https://img.shields.io/github/issues/DanielCondeTorres/YAP-TAZ-TEAD)](https://github.com/DanielCondeTorres/YAP-TAZ-TEAD)
![GitHub pull requests](https://img.shields.io/github/issues-pr/DanielCondeTorres/YAP-TAZ-TEAD)
[![GitHub stars](https://img.shields.io/github/stars/DanielCondeTorres/YAP-TAZ-TEAD)](YAP-TAZ-TEAD/stargazers)
![GitHub repo size](https://img.shields.io/github/repo-size/DanielCondeTorres/YAP-TAZ-TEAD)
[![GitHub license](https://img.shields.io/github/license/DanielCondeTorres/YAP-TAZ-TEAD)](https://github.com/DanielCondeTorres/YAP-TAZ-TEAD)
![GitHub language count](https://img.shields.io/github/languages/count/DanielCondeTorres/YAP-TAZ-TEAD)

<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li><a href="#introduction-movie_camera">Introduction</a></li>
    <li><a href="#general-pre-requirements-computer">Pre-requirements</a></li>
    <li><a href="#required-files-">Required-files</a></li>
    <li><a href="#usage-%EF%B8%8F">Usage</a></li>
    <li><a href="#input-files-">Input-files</a></li>
    <li><a href="#output-files-">Output-files</a></li>
    <li><a href="#analysis-">Analysis</a></li>
    <li><a href="#wiki-">Wiki</a></li>
    <li><a href="#distribution-of-tasks-%EF%B8%8F">Distribution of tasks</a></li>
    <li><a href="#faqs-interrobang">FAQs</a></li>
    <li><a href="#contributing-%EF%B8%8F">Contributing</a></li>
    <li><a href="#cite-mortar_board">Cite</a></li>
    <li><a href="#acknowledgements">Acknowledgements</a></li>
    <li><a href="#license-green_book">License</a></li>
  </ol>
</details>

<!-- Introduction --> 
## Introduction :movie_camera:

### YAP/TAZ and Hippo Signaling Pathway

#### Overview

The YAP/TAZ-TEAD pathway is a key regulatory mechanism in cell growth, proliferation, and apoptosis. It plays a crucial role in organ development, tissue regeneration, and cancer progression.

#### Hippo Pathway Regulation

The Hippo signaling pathway controls the activity of YAP (Yes-associated protein) and TAZ (Transcriptional coactivator with PDZ-binding motif) by phosphorylating them through a kinase cascade:

Active Hippo Pathway: MST1/2 kinases activate LATS1/2 kinases, which phosphorylate YAP/TAZ, leading to their cytoplasmic retention or degradation.

Inactive Hippo Pathway: When Hippo signaling is suppressed, YAP/TAZ translocate into the nucleus, bind to TEAD transcription factors, and activate genes that promote cell proliferation and survival.

#### Key Functions

- Development & Regeneration: YAP/TAZ drive cell proliferation and organ growth.

- Cancer: Hyperactivation of YAP/TAZ leads to uncontrolled cell division, metastasis, and therapy resistance.

- Fibrosis: Persistent YAP/TAZ activation contributes to excessive tissue scarring.

#### Therapeutic Potential

- YAP/TAZ-TEAD inhibitors are being explored for cancer treatment.

- Hippo pathway modulators could aid in tissue repair and fibrosis management.

- Targeting mechanical signals may help regulate YAP/TAZ activity in disease contexts.


### Interaction from literature

#### Summary of Interface Amino Acids YAP-TEAD

|               | **YAP**             | **TEAD**                                |
|--------------|--------------------|-----------------------------------------|
| **Interface 1 $\beta$-Helix** | residues (52–58)     | residues (318–324)                        |
| **Interface 2 $\alpha$-Helix** | residues (Most important: L4,L7,F8 )(61–73)     | residues (345–369)                        |
| **Interface 3 $\alpha$-Helix** | residues (86–100)   (Most important: M25,L30,F34 )  | I(247), V(242), L(272), V(391, Y406, D249, E240, W276, H404) |

From [Li, Ze, et al. "Structural insights into the YAP and TEAD complex." Genes & development 24.3 (2010): 235-240.](https://genesdev.cshlp.org/content/24/3/235.short)
 #### Summary of Interface Amino Acids TAZ-TEAD

|               | **TAZ**                             | **TEAD**                                      |
|--------------|-----------------------------------|----------------------------------------------|
| **Interface 1** | Leu(28), Leu(31), Phe(32), Val(35), Met(36) | Tyr(362), Phe(366), Lys(369), Leu(370), Leu(373), Met(378), Val(382), Phe(386) |
| **Interface 2** | Phe(52), Phe(53), Pro(56)             | Leu(288), Lys(290), Trp(292), Val(407), Gln(418)       |

In () we put the number of the paper the other one is the number sequence.
From [Kaan, Hung Yi Kristal, et al. "Crystal structure of TAZ-TEAD complex reveals a distinct interaction mode from that of YAP-TEAD complex." Scientific reports 7.1 (2017): 2035.](https://www.nature.com/articles/s41598-017-02219-9)

### Goal
<p align="justify"> 
The goal of this repository is to conduct molecular simulations to study a protein-protein complex implicated in cancer and to identify potential drugs to disrupt its formation. 
</p>

#### About this repository
<p align="justify"> 
Given a .pdb file containing two units (or two separate .pdb files), the final .tpr file can be created to run a simulation. Additionally, various analyses can be performed afterward
</p>

<!-- Pre-requirements -->
## General Pre-requirements :computer:

In order to run the program it is necessary to have the following requirements installed:

The program is written in Python language so Python version 3 or higher is required. Also, for python programs to work properly the following libraries are needed:
- [Python](https://www.python.org)
- [PyMOL](https://pymol.org)
- [Gromacs](https://gromacs.bioexcel.eu)
- [Matplotlib](https://matplotlib.org)
- [MDAnalysis](https://www.mdanalysis.org)
- [Mdtraj](https://www.mdtraj.org/1.9.8.dev0/index.html)
- [Numpy](https://numpy.org)
- [Avogadro](https://sourceforge.net/projects/avogadro/files/latest/download)
- [CGenFF](https://app.cgenff.com/login) or [acpype](https://github.com/alanwilter/acpype) 
- [Autodock Vina](https://vina.scripps.edu/downloads/)
- [MGLTools](https://ccsb.scripps.edu/mgltools/downloads/)
- [VMD](https://www.ks.uiuc.edu/Research/vmd/)
- [Re](https://docs.python.org/3/library/re.html)
- [Os](https://docs.python.org/3/library/os.html)
- [Glob](https://docs.python.org/3/library/glob.html)
- [Fileinput](https://docs.python.org/es/3/library/fileinput.html)

### Doing the docking (provisional).
 You can use [Autodock Vina](https://vina.scripps.edu/downloads/) or [AutoDock](https://autodock.scripps.edu) and follow the [tutorial](https://autodock-vina.readthedocs.io/en/latest/docking_basic.html), also install [meeko](https://meeko.readthedocs.io/en/release-doc/installation.html) and [rdkit](https://www.rdkit.org).
 Use [MGLTools](https://ccsb.scripps.edu/mgltools/downloads/) to obtain the .pdbqt file


To be able to use these libraries, it is recommended to work with anaconda environment.


If you are working in an anaconda environment and you need to install a specific module:

```
conda install module
```

or, if you have pip installed:

```
conda install pip
```
```
pip install module
```


<!-- Required-files -->
## Required files and installation📋
Clone this repository. For example,

```
git clone https://github.com/DanielCondeTorres/YAP-TAZ-TEAD.git
```
- Forcefield: **charmm 36**
- Docking software: **Autodock Vina**
- Input files:
    * YAP–TEAD1 complex (PDB: [3KYS](https://www.rcsb.org/structure/3KYS))
    * mYAP-TEAD4 (PDB: [3JUA](https://www.rcsb.org/structure/3JUA))
    * Vgll1-TEAD4  (PDB [5Z2Q](https://www.rcsb.org/structure/5Z2Q))
    * TAZ-TEAD complex (PDB [5GN0](https://www.rcsb.org/structure/5GN0))
- Ligands topology can be obtained from: [CGenFF](https://cgenff.com)

Currently, some PDB files encounter issues when using pdb2gmx, so we perform simulations of the complex using the first two chains of [3JUA](https://www.rcsb.org/structure/3JUA). From this, we extract the TEAD .pdb structure. Meanwhile, the structures of YAP and TAZ are obtained from their amino acid sequences (as detailed in the next section) and predicted using [AlphaFold](https://alphafoldserver.com).


### Input files:
In this directory, you will find different folders, such as:

- **FORCEFIELD**: Contains **charmm36m**, which we use for simulations.  
- **CIF_FILES_FROM_ALPHA_FOLD**: Contains `.cif` files obtained from AlphaFold for **YAP/TAZ** with the short sequence exposed on GitHub, as well as for **human YAP and TAZ**. 
  Use (in **Working_Area**):
   ```
  python cif_to_pdb_converter.py -cif input_file.cif -o output_file.pdb
   ```
- **TEMPLATE_FOR_ALIGNMENT**: Contains **PDB** structures of **YAP-TEAD and TAZ-TEAD** complexes, which can be used to align the structures predicted by AlphaFold.  
- **COMPLEXES**: This folder is subdivided into:  
  - **FROM_PDB**: Contains **PDB** files of complete structures (various subunits of the **YAP/TAZ-TEAD** complex) obtained from **PDB**.  
  - **TWO_SUBUNITS_FORMED**: Contains two previously formed subunits:
    - Files with the prefix **aligned** were aligned using **PyMOL**.  
    - Other files were generated by selecting subunits from different **PDBs**.  
  - **WHOLE_SEQUENCES_FROM_ALPHA_FOLD:** Contains the whole human aminoacid sequence in humans: [YAP](https://www.rcsb.org/structure/AF_AFP46937F1),[TAZ](https://www.alphafold.ebi.ac.uk/entry/Q9GZV5),[TEAD](https://www.uniprot.org/uniprotkb/P28347/entry#structure) : 
  - **YAP_TAZ_TEAD_UNITS**: Contains individual `.pdb` files obtained from **YAP/TAZ** from AlphaFold  and **TEAD** from `3JUA.pdb`, allowing you to create your own free simulations.  
    - These can be aligned with **PyMOL** to a reference structure. 
      -**HUMAN_SEQUENCE_FROM_ALPHA_FOLD:** Contains human YAP/TAZ PDBs predicted by AlphaFold (sequence from [Li, Ze, et al. "Structural insights into the YAP and TEAD complex." Genes & development 24.3 (2010): 235-240.](https://genesdev.cshlp.org/content/24/3/235.short)).  
- **MDP**: Includes .mdp files to run the simulations.
- **LIGANDS**: Includes several *.pdb* and *.top* files that can be used in GROMACS (obtained from [CGenFF](https://app.cgenff.com/login))
  - **For_TEAD**: Ligands that are reported that interacts with TEAD
  - **For_YAP_TAZ**:  Ligands that are reported that interacts with YAP/TAZ

<!-- Usage -->
## Usage ⚙️
In order to run this program, the following command has to be used in the **Working_Area**:

### Only Protein-Complex Simulation.
```
make execute_complex_simultation pdb_complex_1=../Input_files/COMPLEXS/file.pdb 
```
For two units
```
 make execute_complex_simultation pdb_complex_1=file.pdb pdb_complex_2=file_2.pdb
```

This files can be filtered to select specific chains:

```
python obtain_specific_chains.py -pdb file.pdb -c A B -o ../Output/filtered.pdb
```
For example in the case of TAZ-TEAD complex we use:

```
python obtain_specific_chains.py -pdb ../Input_files/COMPLEXS/FORMED_FROM_PDB/5gn0.pdb -c B G -o TAZ_TEAD_2_SUBUNITS_COMPLEX.pdb
```

### Some Useful scripts:

#### PyMOL to template

In [PyMOL](https://pymol.org) terminal.

```
load template.pdb, template
```
 Then we can load our molecules, for example TAZ and TEAD for different .pdb such as:

```
load TAZ.pdb, YAP
```
```
load TEAD.pdb, TEAD
```
If they are in a same .pdb file (complex_not_aligned.pdb) but they are not aligned:

```
load complex_not_aligned.pdb, complex
```
We need to identify the chains of each element

```
create TEAD, complex and chain A
create TAZ, complex and chain B
```
Now we can do the alignment
If TEAD is in our template in chain A and TAZ is in the template in chain B

```
align TEAD, template and chain A
align TAZ, template and chain B
```
Instead of **align** you can use **super** to not break the secondary structure of your protein.

Finally we create a new object to save the created complex! and save the coordinates.

```
# create new_complex_aligned, TEAD or TAZ # Jump to the next, this one can create errors!
save aligned_TAZ_TEAD.pdb,  TEAD or TAZ
```

In one line:

```
pymol -cq -d "load archivo.pdb; align cadena_A, plantilla; align cadena_B, plantilla; save cadenas_alineadas.pdb"
```
#### .cif to .pdb    
From AlphaFold a .cif file will be obtained, this can be easy converted to a .pdb file with:

```
make cif_to_pdb cif_file=input_file.cif  output_file=output_file.pdb
```
or 
```
 python cif_to_pdb_converter.py -cif input_file.cif -o output_file.pdb
 ```


### Obtain the ligand topology
Steps:
- Search the [Smiles Code](https://www.raybiotech.com/) or similar
- Convert it to [.mol](https://www.aatbio.com/tools/smiles-to-mol-converter). If you have access to [ChemDraw3D](https://revvitysignals.com/products/research/chemdraw)
- Use [Avogadro](https://sourceforge.net/projects/avogadro/files/latest/download) to add hydrogens and obtain a good structure and save the file as a .mol2 file of the structure
- Obtain the topolgy from [CGenFF](https://app.cgenff.com/login) or [acpype](https://github.com/alanwilter/acpype) 
> [!NOTE]
Remember to convert the .top file to a valid .itp file

> [!IMPORTANT]
Need to copy in Output, the files from [CGenFF](https://app.cgenff.com/login) in GROMACS (Forcefield, .pdb and .top -which should be change to an .itp format-)

All the steps together with [Autodock Vina](https://vina.scripps.edu/downloads/) and [MGLTools](https://ccsb.scripps.edu/mgltools/downloads/):

**Form out.pdb** from the last frame of your traj

```
make execute_complex_ligand_formation PATH_TO=/Users/danielcondetorres/Applications/mgltools_1.5.7_MacOS-X/MGLToolsPckgs/AutoDockTools/Utilities24 PATH_TO_pythonsh=/Users/danielcondetorres/Applications/mgltools_1.5.7_MacOS-X/bin/pythonsh PDB_INPUT=../HUMAN_3KYS/YAP_TEAD_YAP_FROM_ALPHA_FOLD/out.pdb  RECEPTOR_CHAIN=A LIGAND=../../../YAP-TAZ-TEAD/Code/Input_files/LIGANDS/For_YAP_TAZ/verteporfin_gromacs/verteporfin_gmx22.pdb OUTPUT_DIR=../HUMAN_3KYS/YAP_TEAD_YAP_FROM_ALPHA_FOLD/
```
> [!IMPORTANT]
Inputs:
- **PDB_INPUT**: PDB of your protein complex
- **RECEPTOR_CHAIN**: Chain of PDB input where the ligand will be attach (i.e A, B...)
- **LIGAND**: PDB of you ligand (you can obtain it from [CGenFF](https://app.cgenff.com/login))
- **PATH_TO_pythonsh**: Path to your pythonsh installation (obtained from [MGLTools](https://ccsb.scripps.edu/mgltools/downloads/))
- **PATH_TO**: Path to Utilities24 of [MGLTools](https://ccsb.scripps.edu/mgltools/downloads/), it should be something like ../MGLToolsPckgs/AutoDockTools/Utilities24
- **OUTPUT_DIR**: Directory where files will be saved
Output:
- **complex_ligand.pdb**: PDB with the protein complex and the ligand after the Docking.

The next step is to run an [FEP simulation](https://tutorials.gromacs.org/docs/free-energy-of-solvation.html) in [Gromacs](https://gromacs.bioexcel.eu). This simulation is performed to allow the ligand to gradually incorporate into the protein complex rather than doing so abruptly 

```
make  make execute_complex_ligand_simulation TOPOLOGY=previous_complex_simulation_top/topol.top OUTPUT_DIR=Output_dir_path LIGAND_TOP=path_to_ligand_top/verteporfin_gmx.top    PDB_INPUT_COMPLEX_LIGAND= path_to_complex_ligand.pdb/complex_ligand.pdb FF=path_to_ligand_ff_generated/verteporfin_gromacs/charmm36.ff
```
OJO MODIFICAR EL TOPOL.TOP!!!!!!! puede ser modificando el charmm36 bonded y añadirlo en el mod_topol.py primero :) vamosss
Inputs:
- **PDB_INPUT_COMPLEX_LIGAND**:  PDB with the protein complex and the ligand after the Docking (complex_ligand.pdb in the previous step).

<!-- Output files -->
## Output files 📋
* **Output.txt:**  Contains the ....



> [!NOTE]  
These previous files are stored in the folder output.

<!-- Analysis -->
## Analysis 📊

This section covers the analysis of the complex (go to Working_Area/ANALYSIS)


We recommend focusing on the trajectory; for the initial analyses, water is not necessary. Therefore:
```
gmx trjconv -s prod.tpr -f prod.xtc -o complex.xtc
Select Protein
gmx trjconv -s prod.tpr -f prod.xtc -o complex.pdb -b 0 -e 0
Select Protein
```
MDtraj can present some problems with .pdb files so:
```
pdb_reres complex.pdb > complex_clean.pdb
 ```

 > [!NOTE] 
Installed 
```
pip install pdb-tools or conda install pdb-tools
```
Use this files for analysis_complex.py :)

So far with analysis complex performs the contact matrix between two subunits of the complex and the secondary structure per residue of each subunit:
The inputs are a .pdb of your system and an associated trajectory in .xtc format

```
 python analysis_complex.py -pdb complex.pdb -xtc trajectory.xtc
```
Similarlly, but it needs a .tpr file (in some way the topology), the hydrogen bond matrix between the two subunits can be performed as follows:

```
python hydrogend_bond_matrix.py -tpr ../../../../OUTPUT_CESGA/YAP_TEAD_COMPLEX/prod.tpr -xtc ../../../../OUTPUT_CESGA/YAP_TEAD_COMPLEX/out2.xtc
```

<!-- WIKI -->
## Wiki 📖

Additional information about this modifications can be found in our [Wiki](https://github.com/DanielCondeTorres/YAP-TAZ-TEAD/wiki) **(...Work in progress...)**.


<!-- DISTRIBUTION OF TASKS -->
## Distribution of tasks ✒️ 
Project coordinator: Hugo A. L. Filipe, Ángel Piñeiro and Rebeca García-Fandino

- Main program (Python program): Daniel Conde-Torres,Alejandro Seco, Hugo A. L. Filipe,  Ángel Piñeiro and Rebeca García-Fandino
- GitHub Designer: Daniel Conde-Torres


<!-- FAQs -->
## FAQs :interrobang:



### What ...? 
Answer

<!-- CONTRIBUTING -->
## Contributing 🖇️
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.


# References
[Kaan, Hung Yi Kristal, et al. "Crystal structure of TAZ-TEAD complex reveals a distinct interaction mode from that of YAP-TEAD complex." Scientific reports 7.1 (2017): 2035.](https://www.nature.com/articles/s41598-017-02219-9)

[Li, Ze, et al. "Structural insights into the YAP and TEAD complex." Genes & development 24.3 (2010): 235-240.](https://genesdev.cshlp.org/content/24/3/235.short)

[Crawford, James J., Sarah M. Bronner, and Jason R. Zbieg. "Hippo pathway inhibition by blocking the YAP/TAZ–TEAD interface: a patent review." Expert opinion on therapeutic patents 28.12 (2018): 867-873.](https://www.tandfonline.com/doi/full/10.1080/13543776.2018.1549226?casa_token=IfUAL-qtOCcAAAAA%3Aej1TIFMo8DoScoSl_N3p-RrKfbfgXQawOoM8bVANtZoGxc7gUyqCEdwMrz8Mtz7wST7-wP13Q_bRc9g)

[Pearson, Joel D., et al. "Binary pan-cancer classes with distinct vulnerabilities defined by pro-or anti-cancer YAP/TEAD activity." Cancer Cell 39.8 (2021): 1115-1134.](https://www.cell.com/cancer-cell/fulltext/S1535-6108(21)00338-X)

## References to start with MD simulations

[Braun E, Gilmer J, Mayes HB, Mobley DL, Monroe JI, Prasad S, Zuckerman DM. Best Practices for Foundations in Molecular Simulations [Article v1.0]. Living J Comput Mol Sci. 2019;1(1):5957. doi: 10.33011/livecoms.1.1.5957. Epub 2018 Nov 29. PMID: 31788666; PMCID: PMC6884151](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6884151/)

[Gromacs Tutorial](http://www.mdtutorials.com/gmx/lysozyme/01_pdb2gmx.html)

## Importance of colors to make plots for everybody

[ColorBrewer](https://colorbrewer2.org/#type=sequential&scheme=GnBu&n=5)

## Acknowledgments


To **Cesga**, for allowing us the use of their facilities to hold this seminar.

This work has received financial support from the Spanish Agencia Estatal de Investigación (AEI) and the European Regional Development Fund - ERDF (RTI2018-098795-A-I00, PID2019-111327GB-I00, PDC2022-133402-I00 and RYC-2016-20335), by the (FPU22/00636) and the RePo-SUDOE project.

## Team

Daniel Conde-Torres, Alejandro Seco, Carlos Lozano, António Costa,Hugo A. L. Filipe, Ángel Piñeiro y Rebeca García Fandiño 


## Contact
danielconde.torres@usc.es

alejandro.seco.gonzalez@usc.es

## Social

### LinkedIn

[Daniel Conde-Torres](https://www.linkedin.com/in/daniel-conde-torres-4683b521a)


# STEPS

- Reconocer aminoacidos que forman las interacciones entre el complex que coincida con los articulos
- Intentar alinar el TAZ para que cuadre bien al principio (fijarse en [5GN0](https://www.rcsb.org/structure/5GN0))
- Probar steered molecular dynamics para alinear (poco prioritario)
