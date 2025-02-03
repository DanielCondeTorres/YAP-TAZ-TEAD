# YAP-TAZ-TEAD


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
    <li><a href="#output-files-">Output-files</a></li>
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

<p align="justify"> 
The goal of this repository is to conduct molecular simulations to study a protein-protein complex implicated in cancer and to identify potential drugs to disrupt its formation. 
</p>

<!-- Pre-requirements -->
## General Pre-requirements :computer:

In order to run the program it is necessary to have the following requirements installed:

The program is written in Python language so Python version 3 or higher is required. Also, for python programs to work properly the following libraries are needed:
- [Matplotlib](https://matplotlib.org)
- [Numpy](https://numpy.org)
- [Re](https://docs.python.org/3/library/re.html)
- [Os](https://docs.python.org/3/library/os.html)
- [Mayavi](http://docs.enthought.com/mayavi/mayavi/)
- [Imageio](https://pypi.org/project/imageio/) 
- [Sys](https://docs.python.org/3/library/sys.html)
- [Time](https://docs.python.org/3/library/time.html)
- [Collections](https://docs.python.org/3/library/collections.html)
- [Glob](https://docs.python.org/3/library/glob.html)
- [Fileinput](https://docs.python.org/es/3/library/fileinput.html)
- [Gromacas](https://gromacs.bioexcel.eu)

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
Then you can change the IBM **protein_folding** directory, for the one here!
<!-- Usage -->
## Usage ⚙️
In order to run this program, the following command has to be used in the **Working_Area**:
```
make run  # To perform the search for the most stable state of your amino acid sequence.
```
### In the main.py file, in the inputs section we can choose:
o	File preparation for protein-ligand simulations.

> [!NOTE]  
  INPUTS:
  
    * main_chain: string            # amino acid sequence
    * ws_phase_1: float             # average interaction of an amino acid with solvent in phase 1




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




## References to start with MD simulations

[Braun E, Gilmer J, Mayes HB, Mobley DL, Monroe JI, Prasad S, Zuckerman DM. Best Practices for Foundations in Molecular Simulations [Article v1.0]. Living J Comput Mol Sci. 2019;1(1):5957. doi: 10.33011/livecoms.1.1.5957. Epub 2018 Nov 29. PMID: 31788666; PMCID: PMC6884151](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6884151/)

[Gromacs Tutorial](http://www.mdtutorials.com/gmx/lysozyme/01_pdb2gmx.html)

## Importance of colors to make plots for everybody

[ColorBrewer](https://colorbrewer2.org/#type=sequential&scheme=GnBu&n=5)

## Acknowledgments


To **Cesga**, for allowing us the use of their facilities to hold this seminar.

This work has received financial support from the Spanish Agencia Estatal de Investigación (AEI) and the European Regional Development Fund - ERDF (RTI2018-098795-A-I00, PID2019-111327GB-I00, PDC2022-133402-I00 and RYC-2016-20335), by the (FPU22/00636)-

I would also like to thank Martín Calvelo, for having taught me with so much patience how to perform my first Molecular Dynamics simulations.

## Team

Daniel Conde-Torres, Alejandro Seco, Hugo A. L. Filipe, Ángel Piñeiro y Rebeca García Fandiño 


## Contact
danielconde.torres@usc.es

alejandro.seco.gonzalez@usc.es

## Social

### LinkedIn

[Daniel Conde-Torres](https://www.linkedin.com/in/daniel-conde-torres-4683b521a)
