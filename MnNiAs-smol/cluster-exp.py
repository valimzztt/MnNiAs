

import random
import numpy as np
from monty.serialization import loadfn, dumpfn
from smol.cofe import ClusterSubspace, StructureWrangler, ClusterExpansion, RegressionData
import os 
from pymatgen.core.composition import Element, Composition
from pymatgen.core import Lattice, Structure
from pymatgen.entries.computed_entries import ComputedStructureEntry
from ase.db import connect
from pymatgen.io.ase import AseAtomsAdaptor

"""
    Scripts that demonstrate a basic CLUSTER EXPANSION using SMOL
"""
# directory setup
cwd = os.getcwd()
directory = os.path.join(cwd, "MnNiAs-smol")

# first create the primitive lattice: use fractional occupancies to have a disordered primitive cell
species = [{'Mn': 0.50, 'Ni': 0.5},{'As':1.0}]
my_lattice = Lattice.from_parameters(3.64580405, 3.64580405, 5.04506600, 90, 90, 120)
prim = Structure.from_spacegroup(194, my_lattice,  species, coords=[[0, 0, 0],[0.33333333, 0.66666667, 0.25]], site_properties={"oxidation": [0, 0]})
supercell = prim *(8,8,8)

# Now create a cluster subspace for that structure 
# including pair, triplet and quadruplet clusters up to given cluster diameter cutoffs.
from smol.cofe import ClusterSubspace
cutoffs = {2: 7, 3: 6, 4: 6}  # this is equal to CLEASE'S max_cluster_dia=(7,6,6), list of int or float containing the maximum diameter of clusters (in Å)
# single site and empty orbits are always included
subspace = ClusterSubspace.from_cutoffs(struct, cutoffs=cutoffs)

#the structure wrangler
from monty.serialization import loadfn
from smol.cofe import StructureWrangler

#load the fitting data that was created by looping over VASP output file
energies_file = os.path.join(directory,"comp-struct-energy.json" )
entries = loadfn(energies_file)
wrangler = StructureWrangler(subspace)
for entry in entries:
    wrangler.add_entry(entry)


from sklearn.linear_model import LinearRegression
reg = LinearRegression(fit_intercept=False)
reg.fit(wrangler.feature_matrix, wrangler.get_property_vector("energy"))

from smol.cofe import ClusterExpansion, RegressionData

reg_data = RegressionData.from_sklearn(
    estimator=reg,
    feature_matrix=wrangler.feature_matrix,
    property_vector=wrangler.get_property_vector("energy"),
)
expansion = ClusterExpansion(subspace, coefficients=reg.coef_, regression_data=reg_data)

print(expansion)