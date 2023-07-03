import matplotlib.pyplot as plt
import crystal_toolkit
from pymatgen.core import Lattice, Structure
from smol.capp.generate.special.sqs import StochasticSQSGenerator
from monty.serialization import loadfn, dumpfn
import os
import json
from pymatgen.core.composition import Element, Composition

# directory setup
cwd = os.getcwd()
directory = os.path.join(cwd, "MnNiAs-smol")


# comp1 = Composition({'Mn': 0.60, 'Ni': 0.4})
species = [{'Mn': 0.60, 'Ni': 0.4},{'As':1.0}]
my_lattice = Lattice.from_parameters(3.64580405, 3.64580405, 5.04506600, 90, 90, 120)
struct = Structure.from_spacegroup(194, my_lattice,  species, coords=[[0, 0, 0],[0.33333333, 0.66666667, 0.25]])
prim = struct.get_primitive_structure()

# create a correlation vector based SQS generator
generator_corr = StochasticSQSGenerator.from_structure(
    structure=prim,
    cutoffs={2: 7, 3: 6, 4: 6},  # cluster cutoffs as passed to cluster subspaces
    supercell_size=5,   # the search will be over supercells of 36 atoms
    feature_type="correlation",
    match_weight=1.0,  # weight given to the maximum diameter of perfectly matched vectors (see original publication for details)
)
# generate SQS using correlation vector based score
generator_corr.generate(
    mcmc_steps=100000,  # steps per temperature
    temperatures=None,  # use default, but any sequence of decreasing temperatures can be passed for further control of SA
    max_save_num= None, # the default in this case will be 1000 (1% of mcmc_steps), the actual value of SQSs will likely be much less than that
    progress=True       # show progress bar for each temperature
)


# get SQS structures (removing symmetrically equivalent duplicates)
sqs_corr_list = generator_corr.get_best_sqs(
    num_structures=generator_corr.num_structures,
    remove_duplicates=True,
)

sqs_corr = sqs_corr_list[0]
i = 0

structure_folder = os.path.join(directory,"structures_corr")
os.mkdir(structure_folder)

for sqs_corr in sqs_corr_list:
    i = i +1 
    structure = sqs_corr.structure
    print(structure)
    struct_name = "structure" + str(i) + ".json"
    structure_filename = os.path.join(structure_folder,struct_name)
    with open(structure_filename,'w') as f:
        json.dump(structure.as_dict(), f)