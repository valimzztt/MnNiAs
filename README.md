# Studying the thermodynamics of (Mn,Ni)Sb system
Repository that collects Density Functional Theory (DFT) calculations, Cluster Expansion (CE) and Monte Carlo simulations for MnNiAs. 
Since we only vary the Mn to Ni ratio, while keeping As fixed, the material system can be viewed as pseudo-binary. The reason for choosing Mn-Ni-As as a model system is motivated by the recent discovery of Magnetically Segregated Nanolayering in Mn-Ni-As intermetallics. 

Reference paper:
1. https://pubs.acs.org/doi/pdf/10.1021/acs.chemmater.1c00760

Functionality: 
-------------
The repository currently includes the following functionality:
1. Generate either random or special quasi random structures (SQS) that will provide the fitting data to perform cluster expansion.
2. Perform DFT using VASP on the generated structures: computed energies are stored in an ASE database (of the form SQLite) as well as a ComputedEnergyEntry class provided by the Pymatgen package
3. Perform cluster expansion using CLEASE or SMOL


Required packages: 
-------------
-   **smol** is a minimal implementation of computational methods to calculate
statistical mechanical and thermodynamic properties of crystalline
material systems based on the *cluster expansion* method from alloy theory and
related methods (https://cedergrouphub.github.io/smol/)
-   **pymatgen**  
-   **clease**  

## Important considerations: 
### Choosing between random and SQS:
1. Pros:  SQS generation is the best periodic supercell approximations to the true disordered state for a given number of atoms per supercell. The method is based on a Monte Carlo simulated annealing loop with an objective function that seeks to
perfectly match the maximum number of correlation functions (as opposed to merely minimizing the
distance between the SQS correlation and the disordered state correlations for a pre-specified set of
correlations). SQS are optimal according to the criterion that a specified set of correlations between neighboring
site occupations in the SQS match the corresponding correlation of the true, fully disordered, state.The method optimizes the shape of the supercell jointly with the occupation of
the atomic sites, thus ensuring that the configurational space searched is exhaustive and not biased by a
pre-specified supercell shape. 

2. Cons: literature reported different predicted properties based on different initial SQS. 

### Truncation of the CE expansion: 
 An unsuccessfully truncated CE model due to a limited number of considered training structures may explain the discrepancy in why CE fails to predict the correct outcome in ref. 31. Such inconsistency demonstrates the importance of using an appropriately truncated CE model when studying the mixing of metals in MAX phases.


## Cluster expansion formalism: 
Reference paper: https://www.nature.com/articles/s41524-023-00971-3

An alternative to the computationally demanding crystal structure prediction (CSP) method is the use of methods such as CE which, in contrast to CSP, requires an a priori defined crystal structure. The expansion is carried out on one or multiple sublattices where parameterization is used to express the configurational dependence of physical properties such as energy27, band gap28,29, and magnetic interactions.

### Limitations of CE: 
Despite CE being computationally efficient it is limited by its dependence on an a priori defined input structure, which limits the considered chemical phase space. It would be interesting to combine CE with CSP so that the drawbacks of each method can be eliminated reciprocally. CE and CSP offers an efficient framework that yields reliable results when searching for low-energy basins as the drawbacks of each framework cancel out. The CSP is herein initially used to identify input structures for the CE models whereas the latter is used to explore the mixing and/or stability tendencies in (Mn,Ni)As

The limitations of the cluster expansion formalism when used to pave the path towards stable and possibly synthesizable materials is the restriction defined by the input structure. Relying on an input structure may hinder the exploration of the phase space. Using CE alone is thus prone to miss valuable information hidden within the complete chemical space. An alternative approach that may circumvent this problem, and where the dependence of any initial structure is of lesser importance, or even neglected, is thus preferred. [CITATION NEEDED]
