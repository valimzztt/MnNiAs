# Performing cluster expansion using the SMOL package

The initial structures are generated using the CLEASE method 
ns.generate_gs_structure(), which generates two structures at the extrema and then generates completely random structures

We can change the code so that instead of generating completely random structures, we can generate GS structures

## Workflow will be as follows: 
1. Generate random structures using CLEASE (run generate.py)
2. Perform DFT using VASP on the LISA cluster
3. Retrieve relaxed structure and energy from computed energies (run retrieve-energies.py)
4. Transform ASE Atoms objects into Pymatgen structures (retrieve-energies.py)
5. Create cluster subspace (first fully understand how different parameters for the ClusterSubspace can impact calculated ECI values)    (run cluster-exp.py)
6. Run Monte Carlo simulations on CEDAR cluster