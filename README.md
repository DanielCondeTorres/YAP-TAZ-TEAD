# YAP/TAZ-TEAD


## Expected outcomes

-	Detailed understanding of the protein-protein complex and evaluation of potential next steps. Molecular dynamics simulations evaluating the effectiveness of proposed drugs.
-	An automated GROMACS protocol for protein-ligand simulations.
- An integrated workflow combining docking and molecular dynamics simulations.



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
| **Interface 2 $\alpha$-Helix** | residues 1-12 (Most important: L4,L7,F8 )(61–73)     | residues (345–369)                        |
| **Interface 3 $\alpha$-Helix** | residues 25-39(86–100)   (Most important: M25,L30,F34 )  | I(247), V(242), L(272), V(391, Y406, D249, E240, W276, H404) |

From [Li, Ze, et al. "Structural insights into the YAP and TEAD complex." Genes & development 24.3 (2010): 235-240.](https://genesdev.cshlp.org/content/24/3/235.short)
 #### Summary of Interface Amino Acids TAZ-TEAD

|               | **TAZ**                             | **TEAD**                                      |
|--------------|-----------------------------------|----------------------------------------------|
| **Interface 1** | Leu6(28), Leu9(31), Phe10(32), Val13(35), Met14(36) | Tyr152(362), Phe156(366), Lys159(369), Leu160(370), Leu163(373), Met168(378), Val172(382), Phe176(386) |
| **Interface 2** | Phe30(52), Phe31(53), Pro34(56)             | Leu81(288), Lys83(290), Trp85(292), Val407, Gln418       |

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
- [Re](https://docs.python.org/3/library/re.html)
- [Os](https://docs.python.org/3/library/os.html)
- [Glob](https://docs.python.org/3/library/glob.html)
- [Fileinput](https://docs.python.org/es/3/library/fileinput.html)


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
- Docking software: **  **
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
- **COMPLEX**: This folder is subdivided into:  
  - **FROM_PDB**: Contains **PDB** files of complete structures (various subunits of the **YAP/TAZ-TEAD** complex) obtained from **PDB**.  
  - **TWO_SUBUNITS_FORMED**: Contains two previously formed subunits:
    - Files with the prefix **aligned** were aligned using **PyMOL**.  
    - Other files were generated by selecting subunits from different **PDBs**.  
  - **YAP_TAZ_TEAD_UNITS**: Contains individual `.pdb` files obtained from **YAP/TAZ** from AlphaFold and **TEAD** from `3JUA.pdb`, allowing you to create your own free simulations.  
    - These can be aligned with **PyMOL** to a reference structure. 
      -**HUMAN_SEQUENCE_FROM_ALPHA_FOLD:** Contains human YAP/TAZ PDBs predicted by AlphaFold.  


<!-- Usage -->
## Usage ⚙️
In order to run this program, the following command has to be used in the **Working_Area**:
```
make execute pdb_complex_1=../Input_files/COMPLEXS/file.pdb 
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

### PyMOL to template

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

```
create TEAD, complex and chain A
create TAZ, complex and chain B
```

If TEAD is in our template in chain A and TAZ is in the template in chain B
```
align TEAD, template and chain A
align TAZ, template and chain B
```
Instead of **align** you can use **super** to not break the secondary structure of your protein.

Finally we create a new object to save the created complex! and save the coordinates.

```
create new_complex_aligned, TAZ or TEAD
save aligned_TAZ_TEAD.pdb, new_complex_aligned
```

In one line:

```
pymol -cq -d "load archivo.pdb; align cadena_A, plantilla; align cadena_B, plantilla; save cadenas_alineadas.pdb"
```

### In the main.py file, in the inputs section we can choose:
o	File preparation for protein-ligand simulations.

> [!NOTE]  
  INPUTS:
  
    * main_chain: string            # amino acid sequence
    * ws_phase_1: float             # average interaction of an amino acid with solvent in phase 1

### Protein complex formation:
- Obtain .pdb files of each unit of the complex. 



#### Example
From [AlphaFold](https://alphafoldserver.com), from [UniProt](https://www.uniprot.org/uniprotkb/Q15562/entry), from [PDB](https://www.rcsb.org) or from [310.ai](https://310.ai/copilot/login)
-  .pdb from [YAP](): 
    ```txt
    ETDLEALFNAVMNPKTANVPQTVPMRLRKLPDSFFKPPE

-  .pdb from [TAZ](): 
    ```txt
    PLDTDLEALFNSVMNPKPSSWRKKILPESFFKEPD 

- .pdb from [TEAD]():
    ```txt
    SMRSIASSKLWMLEFSAFLERQQDPDTYNKHLFVHISQSSPSYSDPYLETVDIRQIYDKFPEKKGGLKELFERGPSNAFFLVKFWADLNTNIDDEGSAFYGVSSQYESPENMIITCSTKVCSFGKQVVEKVETEYARYENGHYLYRIHRSPLCEYMINFIHKLKHLPEKYMMNSVLENFTILQVVTNRDTQETLLCIAYVFEVSASEHGAQHHIYRLVKE
    ```
   * Obtained from 3JUA or 
 
    ```txt
    DLNWISMRSIASSKLWMLEFSAFLERQQDPDTYNKHLFVHISQSSPSYSDPYLETVDIRQIYDKFPEKKGGLKELFERGPSNAFFLVKFWADLNTNIDDEGSAFYGVSSQYESPENMIITCSTKVCSFGKQVVEKVETEYARYENGHYLYRIHRSPLCEYMINFIHKLKHLPEKYMMNSVLENFTILQVVTNRDTQETLLCIAYVFEVSASEHGAQHHIYRLVKE
    ```
    * Obtained from 5GN0

> [!NOTE] 
The YAP and TAZ human sequences from [Li, Ze, et al. "Structural insights into the YAP and TEAD complex." Genes & development 24.3 (2010): 235-240.](https://genesdev.cshlp.org/content/24/3/235.short) are, respectively:
```txt
AGHQIVHVRGDSETDETDLEALFNAVMNPKTANVPQTVPMRLRKLPDSFFKPPE
```


```txt
PGQQVIHVTQDLDTDLEALFNSVMNPKPSSWRKKILPESFFKEPD
```


> [!IMPORTANT]
> Here is the table with the amino acid codes in both one-letter and three-letter formats:

| Amino Acid    | 1-Letter Code | 3-Letter Code |
|----------------|---------------|---------------|
| Alanine       | A             | ALA           |
| Cysteine      | C             | CYS           |
| Aspartic acid | D             | ASP           |
| Glutamic acid | E             | GLU           |
| Phenylalanine | F             | PHE           |
| Glycine       | G             | GLY           |
| Histidine     | H             | HIS           |
| Isoleucine    | I             | ILE           |
| Leucine       | L             | LEU           |
| Lysine        | K             | LYS           |
| Methionine    | M             | MET           |
| Asparagine    | N             | ASN           |
| Proline       | P             | PRO           |
| Glutamine     | Q             | GLN           |
| Arginine      | R             | ARG           |
| Serine        | S             | SER           |
| Threonine     | T             | THR           |
| Valine        | V             | VAL           |
| Tryptophan    | W             | TRP           |
| Tyrosine      | Y             | TYR           |



      
From AlphaFold a .cif file will be obtained, this can be easy converted to a .pdb file with:

```
make cif_to_pdb cif_file=input_file.cif  output_file=output_file.pdb
```
or 
```
 python cif_to_pdb_converter.py -cif input_file.cif -o output_file.pdb
 ```

> [!IMPORTANT]
Allows us to modify different parameters related to the creation of the interface.

``` 
*
```             


Therefore, in the first case the importance of the interface will be ``` cls._instance._weight_interface * -1 ``` and in the second case ``` cls._instance._weight_interface * 1 ```.



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
Project coordinator: Ándres Gómez, Ángel Piñeiro and Rebeca García-Fandino

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
