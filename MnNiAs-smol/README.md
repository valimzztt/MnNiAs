# Performing cluster expansion using the SMOL package

## Workflow will be as follows: 
1. Generate random structures using CLEASE (run generate.py) or run the script that generates the special quasi random structures
provided by SMOL (run generate-sqs-clusterpy or generate-sqs-corr.py)
2. Perform DFT using VASP on the LISA cluster
3. Retrieve relaxed structure and energy from computed energies (run retrieve-energies.py)
4. Transform ASE Atoms objects into Pymatgen structures (retrieve-energies.py)
5. Create cluster subspace    (run cluster-exp.py)
6. Run Monte Carlo simulations on CEDAR cluster

## A few notes on SQS generation in smol:

The default generator object will search for SQS structures at the composition given in the disordered structure, so users must make sure the chosen supercell size is compatible: for Mn at 0.6 and Ni at 0.4 (50%) and As at 1.0, one must choose for example supercell_size equal to 5 (20 atoms), 10 atoms (40 atoms) etc.

SQS generation can be carried out by matching correlation vectors, or cluster interaction vectors (with all ECI implicitly set to one, though this can also be overriden).

a. By matching correlation vectors, one can ensure that certain structural correlations observed in the target system are reproduced in the generated SQS.
b. Cluster interaction vectors are similar to correlation vectors but with all ECI (Effective Cluster Interactions) implicitly set to one. ECI represents the influence of different atomic clusters on the stability or energy of the system. By setting all ECI to one, the focus is solely on reproducing the distribution and arrangement of specific atomic clusters in the target system.

Both correlation vectors and cluster interaction vectors are used as input to SQS generation algorithms or techniques. These vectors guide the generation process by constraining the probabilities of certain atomic arrangements and correlations.

The reference paper for SQS in SMOL is: https://www.sciencedirect.com/science/article/abs/pii/S0364591613000540?via%3Dihub

## General workflow for cluster expansion using SMOL

1. Create a ClusterSubspace based on a disordered primitive pymatgen Structure, a given set of diameter cutoffs for clusters, and a specified type of basis set.

2. Use the ClusterSubspace to create a StructureWrangler to generate fitting data in the form of correlation vectors and a normalized property (usually energy). The training data, energy and additional properties are added to the StructureWrangler as pymatgen entries of type ComputedStructureEntry.

3. Fitting data in the form of a correlation StructureWrangler.feature_matrix and a normalized property StructureWrangler.get_property_vector() can be used as input to a linear regression estimator from any choice of third party package, such as scikit-learn, glmnet or sparse-lm.

4. Using the fitted coefficients and the ClusterSubspace instance, a ClusterExpansion is constructed. 

5. Using a ClusterExpansion instance, an Ensemble object can be created to sample the corresponding Hamiltonian for a given supercell size and shape that is specified as a supercell matrix of the unit cell corresponding to the disordered structure used in the first step.

6. Finally, an Ensemble can be sampled in a Monte Carlo simulation by using a Sampler.

## General workflow for cluster expansion using ICET

1. Initialize a cluster space (via ClusterSpace) by providing a prototype structure (typically a primitive cell), the species that are allowed on each site as well as cutoff radii for clusters of different orders

2. Initialize a structure container (via StructureContainer) using the cluster space created previously and add a set of input structures with reference data for the property or properties of interest

3. Fit the parameters using an optimizer (e.g., Optimizer, EnsembleOptimizer, or CrossValidationEstimator from the trainstation package)

4. construct a cluster expansion (via ClusterExpansion) by combining the cluster space with a set of parameters obtained by optimization
