# Performing cluster expansion using the SMOL package

## Workflow will be as follows: 
1. Generate random structures using CLEASE (run generate.py)
2. Perform DFT using VASP on the LISA cluster
3. Retrieve relaxed structure and energy from computed energies (run retrieve-energies.py)
4. Transform ASE Atoms objects into Pymatgen structures (retrieve-energies.py)
5. Create cluster subspace (first fully understand how different parameters for the ClusterSubspace can impact calculated ECI values)    (run cluster-exp.py)
6. Run Monte Carlo simulations on CEDAR cluster

## A few notes on SQS generation in smol:

The default generator object will search for SQS structures at the composition given in the disordered structure, so users must make sure the chosen supercell size is compatible: for Mn at 0.6 and Ni at 0.4 (50%) and As at 1.0, one must choose for example supercell_size equal to 5 (20 atoms), 10 atoms (40 atoms) etc.

SQS generation can be carried out by matching correlation vectors, or cluster interaction vectors (with all ECI implicitly set to one, though this can also be overriden)

The reference paper for SQS in SMOL is: https://www.sciencedirect.com/science/article/abs/pii/S0364591613000540?via%3Dihub


## General workflow for cluster expansion using ICET
The typical workflow involves the following steps:

initialize a cluster space (via ClusterSpace) by providing a prototype structure (typically a primitive cell), the species that are allowed on each site as well as cutoff radii for clusters of different orders

initialize a structure container (via StructureContainer) using the cluster space created previously and add a set of input structures with reference data for the property or properties of interest

fit the parameters using an optimizer (e.g., Optimizer, EnsembleOptimizer, or CrossValidationEstimator from the trainstation package)

construct a cluster expansion (via ClusterExpansion) by combining the cluster space with a set of parameters obtained by optimization