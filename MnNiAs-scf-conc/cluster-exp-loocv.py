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
import warnings 


"""
    Scripts that demonstrate a basic CLUSTER EXPANSION using SMOL using as a primitive cell the one defined by SMOL, which 
    should be equivalent to the one of CLEASE
"""

# directory setup
cwd = os.getcwd()
dft_data = os.path.join(cwd , "dft-data")

# first create the primitive lattice: use fractional occupancies to have a disordered primitive cell
species = [{'Mn': 0.60, 'Ni': 0.40},{'As':1.0}]
my_lattice = Lattice.from_parameters(3.64580405, 3.64580405, 5.04506600, 90, 90, 120)
struct = Structure.from_spacegroup(194, my_lattice,  species, coords=[[0, 0, 0],[0.33333333, 0.66666667, 0.25]], 
                                   site_properties={"oxidation": [0, 0]})
prim = struct.get_primitive_structure()
supercell = struct *(8,8,8)
cwd = os.getcwd()
# this is the database that has the primitive cell

db_name = "clease_MnNiAs-conc.db"
db = connect(db_name)


"""
1. Create a ClusterSubspace based on a disordered primitive pymatgen Structure,
   a given set of diameter cutoffs for clusters, and a specified type of basis set.
   Here we are including pair, triplet and quadruplet clusters up to given cluster diameter cutoffs (in Å) of 7, 6, 6, respectively, 
   and single site and empty orbits are always included
   -----------------------------------------------------------------------------------------------------------------------------
   A ClusterSubspace is the main work horse used in constructing a cluster expansion. 
   A cluster subspace holds a finite set of orbits that contain symmetrically equivalent clusters. 
   The orbits also contain the set of orbit basis functions (also known as correlation functions) that represent 
   the terms in the cluster expansion. Taken together the set of all orbit functions for all orbits included span 
   a subspace of the total function space over the configurational space of a given crystal structure system.
"""


# Most straightforward way of generating a clusterSubspace.from_cutoffs
# This will auto-generate the orbits from diameter cutoffs
cutoffs = {2: 7, 3: 7, 4: 7} 
subspace = ClusterSubspace.from_cutoffs(prim, cutoffs=cutoffs, use_concentration=True, supercell_size="As")
warnings.warn("We have created the cluster subspace")

"""
2.  Use the ClusterSubspace to create a StructureWrangler
    to generate fitting data in the form of correlation vectors and a normalized property (usually energy)
    We will load as fitting data the energy of each DFT Vasp run that we previously stored in a json file 
-----------------------------------------------------------------------------------------------------------------------------
    The StructureWrangler handles input data structures and properties to fit to the cluster expansion. 
"""
from monty.serialization import loadfn
from smol.cofe import StructureWrangler

energies_file = os.path.join(dft_data ,"struct-energy-scf.json" )
entries = loadfn(energies_file)
wrangler = StructureWrangler(subspace)
for entry in entries:
    wrangler.add_entry(entry)
# The verbose flag will print structures that fail to match.
warnings.warn(f'\nTotal structures that match {wrangler.num_structures}/{len(entries)}')

"""
3. Fitting data in the form of a correlation StructureWrangler.feature_matrix and a normalized property
  StructureWrangler.get_property_vector() can be used as
 input to a linear regression estimator from any choice of third party package, such as scikit-learn, glmnet or sparse-lm."""
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import mean_squared_error, max_error
estimator = 0
init_cv = 1e6
e = RidgeCV(alphas=(1e-6, 100, 1000), fit_intercept=False, store_cv_values=True, gcv_mode = 'eigen', alpha_per_target = True)
X = wrangler.feature_matrix
y = wrangler.get_property_vector("energy")
e.fit(X=X, y=y)
if (-1*e.best_score_ < init_cv):
    estimator = e
    init_cv = -1*e.best_score_ 

coefs = estimator.coef_

train_predictions = np.dot(wrangler.feature_matrix, coefs)

rmse = mean_squared_error(
    wrangler.get_property_vector("energy"), train_predictions, squared=False
)
maxer = max_error(wrangler.get_property_vector("energy"), train_predictions)

# if rmse <= init_error:
#     true_est = estimator
#     coefficients = coefs
#     init_error = rmse
#     init_max_error = maxer

reg_data = RegressionData.from_sklearn(
    estimator, wrangler.feature_matrix, wrangler.get_property_vector("energy")
)
expansion = ClusterExpansion(
    subspace, coefficients=coefs, regression_data=reg_data
)
print(expansion)
structure = random.choice(wrangler.structures)
prediction = expansion.predict(structure, normalized=True)
print(
    f"The predicted energy for a structure with composition "
    f"{structure.composition} is {prediction} eV/prim.\n"
)
print(f"The fitted coefficients are:\n{expansion.coefs}\n")
print(f"The effective cluster interactions are:\n{expansion.eci}\n")

from smol.io import save_work

file_path = 'ce-data/ce_MnNiAs_scf_loocv.mson'
# We can save the subspace as well, but since both the wrangler
# and the expansion have it, there is no need to do so.
save_work(file_path, wrangler, expansion)
