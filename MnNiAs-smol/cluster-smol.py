
# Scripts that demonstrate a basic CLUSTER EXPANSION using smol
import random
import numpy as np
import os 
from monty.serialization import loadfn, dumpfn
from pymatgen.core.structure import Structure
from pymatgen.core.structure import Lattice
from smol.cofe import ClusterSubspace, StructureWrangler, ClusterExpansion, RegressionData

# first create the lattice
from pymatgen.core.composition import Element, Composition
from pymatgen.core import Lattice, Structure
# [{"element": "Mn", "oxidation_state": 0.0, "occu": 0.5},{"element": "Ni", "oxidation_state": 0.0, "occu": 0.5},{"element": "Ni", "oxidation_state": 0.0, "occu": 0.5} ]
species = [{'Mn': 0.60, 'Ni': 0.4},{'As':1.0}]
#species = [[{"element": "Mn", "oxidation_state": 0.0, "occu": 0.6},{"element": "Ni", "oxidation_state": 0.0, "occu": 0.4}],
# {"element": "As", "oxidation_state": 0.0, "occu": 0.5} ]


my_lattice = Lattice.from_parameters(3.64580405, 3.64580405, 5.04506600, 90, 90, 120)
prim = Structure.from_spacegroup(194, my_lattice,  species, coords=[[0, 0, 0],[0.33333333, 0.66666667, 0.25]])
supercell = prim *(8,8,8) #where supercell_factor=8

# Now create a cluster subspace for that structure 
# including pair, triplet and quadruplet clusters up to given cluster diameter cutoffs.
from smol.cofe import ClusterSubspace


cutoffs = {2: 7, 3: 6, 4: 6}
print("Before creating the subspace")
subspace = ClusterSubspace.from_cutoffs(supercell,
                                        cutoffs=cutoffs, # will include orbits of 2 and 3 sites.
                                        basis='indicator', # sets the site basis type, default is indicator
                                        supercell_size='As')
# supercell_size specifies the method to determine the supercell size
# when trying to match a structure.
# (See pymatgen.structure_matcher.StructureMatcher for more info)
print(subspace) # single site and empty orbits are always included

#the structure wrangler
from monty.serialization import loadfn
from smol.cofe import StructureWrangler

directory = os.path.join(os.getcwd(), "smol-data")
dft_data = os.path.join(directory, "dft-data")
energy_file  = os.path.join(dft_data,"energies.json")
entries = loadfn(energy_file)
wrangler = StructureWrangler(subspace)
for entry in entries:
    wrangler.add_entry(entry)


# The verbose flag will print structures that fail to match.
print(f'\nTotal structures that match {wrangler.num_structures}/{len(entries)}')

from sklearn.linear_model import LinearRegression
# Set fit_intercept to False because we already do this using
# the empty cluster.
estimator = LinearRegression(fit_intercept=False)
estimator.fit(wrangler.feature_matrix,
              wrangler.get_property_vector('energy'))
coefs = estimator.coef_

from sklearn.metrics import mean_squared_error, max_error

train_predictions = np.dot(wrangler.feature_matrix, coefs)
rmse = mean_squared_error(wrangler.get_property_vector('energy'),
                          train_predictions, squared=False)
maxer = max_error(wrangler.get_property_vector('energy'),
                  train_predictions)

print(f'RMSE {1E3 * rmse} meV/prim')
print(f'MAX {1E3 * maxer} meV/prim')


reg_data = RegressionData.from_sklearn(
    estimator, wrangler.feature_matrix,
    wrangler.get_property_vector('energy'))


expansion = ClusterExpansion(subspace,
                             coefficients=coefs,
                             regression_data=reg_data)

structure = random.choice(wrangler.structures)
prediction = expansion.predict(structure)

print(f'The predicted energy for a structure with composition '
      f'{structure.composition} is {prediction} eV/prim.\n')
print(f'The fitted coefficients are:\n{expansion.coefs}\n')
print(f'The effective cluster interactions are:\n{expansion.eci}\n')

from smol.io import save_work

file_path = 'MnNiAs-smol/ce_MnNiAs.mson'
# we can save the subspace as well, but since both the wrangler
# and the expansion have it, there is no need to do so.
save_work(file_path, wrangler, expansion)

