from clease.settings import Concentration
import os 
from ase.io import vasp, db
from pymatgen.io.ase import AseAtomsAdaptor
from clease.settings import Concentration
from scipy.constants import k, electron_volt
# Run a Monte Carlo adding constraints? 
import os
import json
from pymatgen.core.structure import Structure
from smol.io import load_work
from smol.moca import Ensemble
from smol.moca import Sampler
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.transformations.standard_transformations import (
    OrderDisorderedStructureTransformation,
)
import warnings


# These contain the ECI values obtained by fitting the SMOL data
file_path = 'ce-data/ce_MnNiAs_scf.mson'
warnings.warn("We are starting the job")
work = load_work(file_path)
warnings.warn("We have loaded the job")
for name, obj in work.items():
    print(f'{name}: {type(obj)}\n')
expansion  = work.get("ClusterExpansion")
from pymatgen.core.composition import Composition
from pymatgen.core.periodic_table import Element

sc_matrix = [[8, 0, 0], [0, 8, 0], [0, 0, 8]]


ensemble = Ensemble.from_cluster_expansion(expansion, supercell_matrix=sc_matrix)
# This will take care of setting the defaults
# for the supplied canonical ensemble
sampler = Sampler.from_ensemble(ensemble, temperature=2000)
warnings.warn(f"Sampling information: {sampler.samples.metadata}")

from pymatgen.transformations.standard_transformations import (
    OrderDisorderedStructureTransformation,
)

transformation = OrderDisorderedStructureTransformation(no_oxi_states=True)
structure = expansion.cluster_subspace.structure.copy()
structure.make_supercell(sc_matrix)
warnings.warn("About to apply transformation")
structure = transformation.apply_transformation_fast(structure)

print(structure)

import random

def randomizer():
    if random.random() < 0.6:
        return 1
    else:
        return 0

for i in range(0,len(structure.sites),1):
    value = randomizer()
    if(structure[i].species.get_atomic_fraction(Element("Mn"))==1.0 and value == 1):
        structure[i]="Mn"
    if(structure[i].species.get_atomic_fraction(Element("Mn"))==1.0 and value  == 0):
        structure[i]="Ni"

print(structure)
from ase.build import sort
from ase.io import vasp
atoms1 = AseAtomsAdaptor.get_atoms(structure)
sorted_atoms=sort(atoms1)
vasp.write_vasp('SORTED.vasp', sorted_atoms, direct=False, wrap=False)
vasp.write_vasp('INITIAL.vasp', atoms1 , direct=False, wrap=False)
init_occu = ensemble.processor.occupancy_from_structure(structure)

#atoms = AseAtomsAdaptor.get_atoms(structure)
T = 2000
# convert current energy structure to ASE atoms 
from pymatgen.io.xyz import XYZ
from pymatgen.io.cif import CifWriter

#vasp.write_vasp('POSCARsorted{}.vasp'.format(T), sorted_atoms, direct=False, wrap=False)
#sorted_atoms=sort(sorted_atoms)
writer = CifWriter(structure)
writer.write_file('POSCAR_smol_cif{}.cif'.format(T))

directory = 'MC_2000K_smol3'
cwd = os.getcwd()
parent_dir = cwd
path = os.path.join(parent_dir, directory)

if not os.path.exists(directory):
    os.makedirs(directory)
os.chdir(path)

mcmc_steps = 100
T = 2000
heat_caps = []
kb = k / electron_volt

sampler.run(mcmc_steps,init_occu,
            thin_by=1000, # thin_by will save every 100th sample only
            progress=False) # progress will show progress bar

f = open('MC_eq_{}'.format(T),'w')
# You can get the minimum energy structure and current structure
# by using the ensemble processor
samples = sampler.samples
print(samples)
print("These are the occupancies")
print(samples.get_occupancies())
print(samples.get_minimum_energy_occupancy())
curr_s = ensemble.processor.structure_from_occupancy(samples.get_occupancies()[-1])
min_s = ensemble.processor.structure_from_occupancy(samples.get_minimum_energy_occupancy())
# convert mininum energy structure to ASE atoms 
writer = CifWriter(min_s)
writer.write_file('RELAXED_MIN{}.cif'.format(T))
writer = CifWriter(curr_s)
writer.write_file('CURR_s{}.cif'.format(T))
# save minimum energy among all samples
min_e = samples.get_minimum_energy()
# convert current energy structure to ASE atoms 

# save minimum energy among all samples
min_e = samples.get_minimum_energy()
print(min_e)
# save average energy of samples
avg_e = samples.mean_energy()
mean_composition = samples.mean_composition()
energy_var = samples.energy_variance()
heat_capacity = energy_var/(kb*T**2)
thermo = {'temperature': T, 'minimum_energy':min_e, 'avg_energy': avg_e,
             'energy_var': energy_var, 'heat_capacity': heat_capacity}
print(thermo, file=f)
f.close()

# this is the minimum energy after equilibrium at 2000K
with open('MC_2000K_eq.json','w') as f:
    json.dump(min_s.as_dict(), f)

