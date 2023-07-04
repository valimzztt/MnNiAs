# Performing cluster expansion using the SMOL package

## Workflow will be as follows: 
1. Generate random structures using CLEASE (run generate.py) or run the script that generates the special quasi random structures
provided by SMOL (run generate-sqs-clusterpy or generate-sqs-corr.py)
2. Perform DFT using VASP on the LISA cluster
3. Retrieve relaxed structure and energy from computed energies (run retrieve-energies.py)
4. Transform ASE Atoms objects into Pymatgen structures (retrieve-energies.py)
5. Create cluster subspace (first fully understand how different parameters for the ClusterSubspace can impact calculated ECI values)    (run cluster-exp.py)
6. Run Monte Carlo simulations on CEDAR cluster

## A few notes on SQS generation in smol:

The default generator object will search for SQS structures at the composition given in the disordered structure, so users must make sure the chosen supercell size is compatible: for Mn at 0.6 and Ni at 0.4 (50%) and As at 1.0, one must choose for example supercell_size equal to 5 (20 atoms), 10 atoms (40 atoms) etc.

SQS generation can be carried out by matching correlation vectors, or cluster interaction vectors (with all ECI implicitly set to one, though this can also be overriden)

The reference paper for SQS in SMOL is: https://www.sciencedirect.com/science/article/abs/pii/S0364591613000540?via%3Dihub

