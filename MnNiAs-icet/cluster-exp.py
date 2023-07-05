from ase.db import connect
from pymatgen.io.ase import AseAtomsAdaptor
from monty.serialization import loadfn, dumpfn
import os
from ase.db import connect
from icet import ClusterSpace, StructureContainer, ClusterExpansion
from trainstation import CrossValidationEstimator
from ase.calculators.vasp import Vasp
from ase.db import connect
from clease.tools import update_db
import warnings
# This scripts runs in about 16 seconds on an i7-6700K CPU.



cwd = os.getcwd()
directory = os.path.join(cwd, "MnNiAs-icet")
# this database contains the DFT computed energies 
database = os.path.join(directory,'clease_MnNiAs-conc.db')
db = connect(database)
database_prim = os.path.join(directory,'clease_MnNiAs-conc.db')
db_prim = connect(database_prim)

"""Extracts energy info given Atoms object given row id id once DFT is done"""
def extract_energy(row_id):
    # get the id of the final structure whose calculation converged
    with db.managed_connection() as con:
            cur = con.cursor()
            cur.execute("SELECT energy FROM systems WHERE id = ? AND energy IS NOT NULL", (row_id,))
            # Fetch the result
            result = cur.fetchone()
            if result != None: 
                energy = result[0]
                return energy

# this database contains the primitive cell 
primitive_structure = db_prim.get(id=1).toatoms()  # primitive structure

cs = ClusterSpace(structure=primitive_structure,
                  cutoffs=[4.0, 4.0, 4.0],
                  chemical_symbols=['Mn', 'Ni', 'As'])

sc = StructureContainer(cluster_space=cs)

# step 3: Parse the input structures and set up a structure container
for row in db.select():
    atoms = row.toatoms()
    energy = extract_energy(row.id) 
    if(energy != None):
        sc.add_structure(structure=row.toatoms(),
                        user_tag="valmzztt",
                        properties={'energy': energy})
warnings.warn(str(sc))

# step 4: Train parameters
opt = CrossValidationEstimator(fit_data=sc.get_fit_data(key='energy'), fit_method='ardr')
opt.validate()
opt.train()
warnings.warn(str(opt))

# step 5: Construct cluster expansion and write it to file
ce = ClusterExpansion(cluster_space=cs, parameters=opt.parameters, metadata=opt.summary)
print(ce)
warnings.warn(str(ce))
ce.write('energy.ce')