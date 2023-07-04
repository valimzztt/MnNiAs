
# Run a Monte Carlo adding constraints? 
import numpy as np
import os
import json
from pymatgen.core.structure import Structure
from smol.io import load_work
from smol.cofe import ClusterSubspace, StructureWrangler, ClusterExpansion, RegressionData
from smol.io import save_work
from smol.moca import Ensemble
from smol.moca import Sampler
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.transformations.standard_transformations import (
    OrderDisorderedStructureTransformation,
)
import warnings

# These contain the ECI values obtained by fitting the CLEASE data
file_path = 'ce-data/ce_MnNiAs_clease.mson'
warnings.warn("We are starting the job")
work = load_work(file_path)
warnings.warn("We have loaded the job")
for name, obj in work.items():
    print(f'{name}: {type(obj)}\n')
expansion  = work.get("ClusterExpansion")
#wrangler = work.get("StructureWrangler")

sc_matrix = [[6, 0, 0], [0, 6, 0], [0, 0, 6]]
ensemble = Ensemble.from_cluster_expansion(expansion, supercell_matrix=sc_matrix)

# This will take care of setting the defaults
# for the supplied canonical ensemble
sampler = Sampler.from_ensemble(ensemble, temperature=3000)
warnings.warn(f"Sampling information: {sampler.samples.metadata}")


transformation = OrderDisorderedStructureTransformation(no_oxi_states=True)
# this gets you the primitive structure associated to cluster expansion 
structure = expansion.cluster_subspace.structure.copy()
warnings.warn("Applying transformation")
structure = transformation.apply_transformation(structure)
warnings.warn("We have applied transformation")
structure.make_supercell(sc_matrix)

with open('transformed_structure1.json','w') as f:
    json.dump(structure.as_dict(), f)

init_occu = ensemble.processor.occupancy_from_structure(structure)
temperatures = np.linspace(3000, 0, 10)

from ase.build import sort
from ase.io import vasp
# Monte Carlo results will be created in new directory

directory = 'MC_1500K_clease'
cwd = os.getcwd()
parent_dir = cwd
path = os.path.join(parent_dir, directory)
os.mkdir(path)  
os.chdir(path)

mcmc_steps = 1000000
T = 1500
sampler.run(mcmc_steps,init_occu,
            thin_by=100, # thin_by will save every 100th sample only
            progress=False) # progress will show progress bar

f = open('MC_eq_{}'.format(T),'w')
# You can get the minimum energy structure and current structure
# by using the ensemble processor
samples = sampler.samples
curr_s = ensemble.processor.structure_from_occupancy(samples.get_occupancies()[-1])
min_s = ensemble.processor.structure_from_occupancy(samples.get_minimum_energy_occupancy())
# convert mininum energy structure to ASE atoms 
atoms = AseAtomsAdaptor.get_atoms(min_s)
sorted_atoms=sort(atoms)
vasp.write_vasp('POSCARmin1500.vasp', sorted_atoms, direct=False, wrap=False)
# convert current energy structure to ASE atoms 
atoms = AseAtomsAdaptor.get_atoms(curr_s)
sorted_atoms=sort(atoms)
vasp.write_vasp('POSCARcurr1500.vasp', sorted_atoms, direct=False, wrap=False)
# save minimum energy among all samples
min_e = samples.get_minimum_energy()
# save average energy of samples
avg_e = samples.mean_energy()
mean_composition = samples.mean_composition()
energy_var = samples.energy_variance()
thermo = {'temperature': T, 'minimum_energy':min_e, 'avg_energy': avg_e,
             'energy_var': energy_var}
print(thermo, file=f)
f.close()

# this is the minimum energy after equilibrium at 2000K
with open('struct_1500K_eq.json','w') as f:
    json.dump(min_s.as_dict(), f)


for i in range(1500, 0, -50):
    T = i
    sampler.run(mcmc_steps,init_occu,
            thin_by=100, # thin_by will save every 100th sample only
            progress=False) # progress will show progress bar
    # Samples are saved in a sample container
    samples = sampler.samples
    # You can get the minimum energy structure and current structure
    # by using the ensemble processor
    curr_s = ensemble.processor.structure_from_occupancy(samples.get_occupancies()[-1])
    min_s = ensemble.processor.structure_from_occupancy(samples.get_minimum_energy_occupancy())
    # convert mininum energy structure to ASE atoms 
    atoms = AseAtomsAdaptor.get_atoms(min_s)
    sorted_atoms=sort(atoms)
    vasp.write_vasp('POSCAR{}.vasp'.format(i), sorted_atoms, direct=False, wrap=False)
    POSCARstring='POSCAR{}.vasp'.format(i)
    # save minimum energy among all samples
    min_e = samples.get_minimum_energy()
    # save average energy of samples
    avg_e = samples.mean_energy()
    mean_composition = samples.mean_composition()
    energy_var = samples.energy_variance()
    thermo = {'temperature': T, 'minimum_energy':min_e, 'avg_energy': avg_e,
             'energy_var': energy_var}
    f = open('MC_{}'.format(i),'w')
    print(thermo, file=f)
    f.close()
os.chdir(cwd)
   
# save_work("MnNiAs_montecarlo.json", wrangler, expansion, ensemble, sampler.samples)