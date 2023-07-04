# This scripts runs in about 16 seconds on an i7-6700K CPU.

from ase.db import connect
from icet import ClusterSpace, StructureContainer, ClusterExpansion
from trainstation import CrossValidationEstimator
import ase 


# this database contains the primitive cell 
db = connect('MnNiAs-clease/clease_MnNiAs.db')
primitive_structure = db.get(id=1).toatoms()  # primitive structure

cs = ClusterSpace(structure=primitive_structure,
                  cutoffs=[7.0, 6.0, 6.0],
                  chemical_symbols=['Mn', 'Ni', 'As'])

sc = StructureContainer(cluster_space=cs)

# this database contains the DFT computed energies 
db = connect('MnNiAs-clease/clease_MnNiAs_energy.db')
for row in db.select():
    atoms = row.toatoms()
    print(atoms.calc)
    if(atoms.calc != None):
        print(atoms.calc.energy)
    
print(sc)