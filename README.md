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


### Goal
<p align="justify"> 
The goal of this repository is to conduct molecular simulations to study a protein-protein complex implicated in cancer and to identify potential drugs to disrupt its formation. 
</p>

#### About this repository
<p align="justify"> 
Escribiendo...
</p>

<!-- Pre-requirements -->
## General Pre-requirements :computer:

In order to run the program it is necessary to have the following requirements installed:

The program is written in Python language so Python version 3 or higher is required. Also, for python programs to work properly the following libraries are needed:
- [Matplotlib](https://matplotlib.org)
- [MDAnalysis](https://www.mdanalysis.org)
- [Mdtraj](https://www.mdtraj.org/1.9.8.dev0/index.html)
- [Numpy](https://numpy.org)
- [Re](https://docs.python.org/3/library/re.html)
- [Os](https://docs.python.org/3/library/os.html)
- [Glob](https://docs.python.org/3/library/glob.html)
- [Fileinput](https://docs.python.org/es/3/library/fileinput.html)
- [Gromacs](https://gromacs.bioexcel.eu)

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
git clone https://github.com/qiskit-research/qiskit-research.git

```
- Forcefield: **charmm 36**
- Docking software: **  **
- Input files:
    * YAP–TEAD1 complex (PDB: [3KYS](https://www.rcsb.org/structure/3KYS))
    * mYAP-TEAD4 (PDB: [3JUA](https://www.rcsb.org/structure/3JUA))
    * Vgll1-TEAD4  (PDB [5Z2Q](https://www.rcsb.org/structure/5Z2Q))
    * TAZ-TEAD complex (PDB [5GN0](https://www.rcsb.org/structure/5GN0))
- Ligands topology can be obtained from: [CGenFF](https://cgenff.com)
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
    MEPSSWSGSESPAENMERMSDSADKPIDNDAEGVWSPDIEQSFQEALAIYPPCGRRKIIL
    SDEGKMYGRNELIARYIKLRTGKTRTRKQVSSHIQVLARRKSRDFHSKLKDQTAKDKALQ
    HMAAMSSAQIVSATAIHNKLGLPGIPRPTFPGAPGFWPGMIQTGQPGSSQDVKPFVQQAY
    PIQPAVTAPIPGFEPASAPAPSVPAWQGRSIGTTKLRLVEFSAFLEQQRDPDSADLNCNI
    QDDAGAFYGVTSQYESSENMTVTCSTKVCSFGKQVVEKVETEYARFENGRFVYRINRSPM
    CEYMINFIHKLKHLPEKYMMNSVLENFTILLVVTNRDTQETLLCMACVFEVSNSEHGAQH
    HIYRLVKD
 
 
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
