import matplotlib.pyplot as plt
import crystal_toolkit
from pymatgen.core import Lattice, Structure
from smol.capp.generate.special.sqs import StochasticSQSGenerator
from monty.serialization import loadfn, dumpfn
from pymatgen.core.composition import Element, Composition
from pymatgen.core import Lattice, Structure
import json 
import os

"""
    Script that generates the special quasirandom structures (via a cluster vector based SQS generator)
    for the (Mn,Ni)Sb system with specific concentration and
    stores them as Structure objects inside the structures folder
"""
# directory setup
cwd = os.getcwd()
directory = os.path.join(cwd, "MnNiAs-smol")

species = [{'Mn': 0.60, 'Ni': 0.4},{'As':1.0}]
my_lattice = Lattice.from_parameters(3.64580405, 3.64580405, 5.04506600, 90, 90, 120)
struct = Structure.from_spacegroup(194, my_lattice,  species, coords=[[0, 0, 0],[0.33333333, 0.66666667, 0.25]])
prim = struct.get_primitive_structure()

# create a cluster interaction based SQS generator
generator_cint = StochasticSQSGenerator.from_structure(
    structure=prim,
    cutoffs={2: 7, 3: 6, 4: 6},  # cluster cutoffs as passed to cluster subspaces
    supercell_size=5,
    feature_type="cluster-interaction",
    match_weight=1.0,
)

# generate SQS using cluster interaction vector based score
generator_cint.generate(
    mcmc_steps=100000, # steps per temperature
    temperatures=None,  # use default, but any sequence of decreasing temperatures can be passed for further control of SA
    max_save_num= None, # the default in this case will be 1000 (1% of mcmc_steps), the actual value of SQSs will likely be much less than that
    progress=True       # show progress bar for each temperature
)

# get SQS structures (removing symmetrically equivalent duplicates)
sqs_cint_list = generator_cint.get_best_sqs(
    num_structures=generator_cint.num_structures,
    remove_duplicates=True,
)


# the lists are sorted as we can see in the plot above
sqs_cint = sqs_cint_list[0]
i = 0
structure_folder = os.path.join(directory,"structures_cint")
os.mkdir(structure_folder)

for sqs_cint in sqs_cint_list:
    i = i +1 
    structure = sqs_cint.structure
    print(structure)
    struct_name = "structure" + str(i) + ".json"
    structure_filename = os.path.join(structure_folder,struct_name)
    with open(structure_filename,'w') as f:
        json.dump(structure.as_dict(), f)