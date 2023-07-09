
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
file_path = 'ce-data/ce_MnNiAs_smol.mson'
warnings.warn("We are starting the job")
work = load_work(file_path)
warnings.warn("We have loaded the job")
for name, obj in work.items():
    print(f'{name}: {type(obj)}\n')
expansion  = work.get("ClusterExpansion")
from pymatgen.core.composition import Composition
from pymatgen.core.periodic_table import Element

sc_matrix = [[6, 0, 0], [0, 6, 0], [0, 0, 6]]
structure = expansion.cluster_subspace.structure.copy()*(6,6,6)

for i in range(0,len(structure.sites),1):
    print(structure[i].species)
    print(structure[i].species.get_atomic_fraction(Element("Mn")))
    if(structure[i].species.get_atomic_fraction(Element("Mn"))==0.6 and i%2 == 0):
        structure[i]="Mn"
    if(structure[i].species.get_atomic_fraction(Element("Mn"))==0.6 and i%2 != 0):
        structure[i]="Ni"

print(structure)
ensemble = Ensemble.from_cluster_expansion(expansion, supercell_matrix=sc_matrix)


from pymatgen.io.xyz import XYZ
from pymatgen.io.cif import CifWriter

#vasp.write_vasp('POSCARsorted{}.vasp'.format(T), sorted_atoms, direct=False, wrap=False)
#sorted_atoms=sort(sorted_atoms)
T = 2000
writer = CifWriter(structure)
writer.write_file('POSCARNi_cif{}.cif'.format(T))
