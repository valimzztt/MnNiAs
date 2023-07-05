
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

sc_matrix = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
structure = expansion.cluster_subspace.structure.copy()
print(structure)
ensemble = Ensemble.from_cluster_expansion(expansion, supercell_matrix=sc_matrix)

# This will take care of setting the defaults
# for the supplied canonical ensemble
sampler = Sampler.from_ensemble(ensemble, temperature=2000)
warnings.warn(f"Sampling information: {sampler.samples.metadata}")



#transformation = OrderDisorderedStructureTransformation(no_oxi_states=True)
# structure = transformation.apply_transformation(structure)
structure.make_supercell(sc_matrix)
print(structure)

for i in range(0,len(structure.sites),1):
    if(structure[i].species.get_atomic_fraction(Element("Mn"))==0.6 and i%2 == 0):
        structure[i]="Ni"
    if(structure[i].species.get_atomic_fraction(Element("Mn"))==0.6 and i%2 != 0):
        structure[i]="Ni"


# this gets you the primitive structure associated to cluster expansion 

with open('transformed_structure.json','w') as f:
    json.dump(structure.as_dict(), f)

print(structure)
init_occu = ensemble.processor.occupancy_from_structure(structure)
from ase.build import sort
from ase.io import vasp
#atoms = AseAtomsAdaptor.get_atoms(structure)
T = 2000
# convert current energy structure to ASE atoms 
from pymatgen.io.xyz import XYZ
from pymatgen.io.cif import CifWriter

#vasp.write_vasp('POSCARsorted{}.vasp'.format(T), sorted_atoms, direct=False, wrap=False)
#sorted_atoms=sort(sorted_atoms)
writer = CifWriter(structure)
writer.write_file('POSCARNi_cif{}.cif'.format(T))