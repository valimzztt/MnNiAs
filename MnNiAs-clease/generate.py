"""
    Script that generates random structures using the CLEASE package
    The initial structures are generated using the CLEASE method  ns.generate_gs_structure(), which generates two structures at the extrema and then generates completely random structures
    We can change the code so that instead of generating completely random structures, we can generate GS structures
"""
#setup directory
import os
import clease

#define concentration range of all elements
from clease.settings import Concentration
conc = Concentration(basis_elements=[['Mn', 'Ni'], ['As']])
conc.set_conc_ranges(ranges=[[(0,1),(0,1)], [(1,1)]])


#define crystal structure
from clease.settings import CECrystal
settings = CECrystal(concentration=conc,
    spacegroup=194,
    basis=[(0.00000, 0.00000, 0.00000), (0.33333333, 0.66666667, 0.25)],
    cell=[3.64580405, 3.64580405,   5.04506600, 90, 90, 120],
    supercell_factor=8,
    db_name="clease_MnNiAs.db",
    basis_func_type='binary_linear',
    max_cluster_dia=(7,6,6))


#generate first round of structures
from clease import NewStructures
ns = NewStructures(settings, generation_number=0, struct_per_gen=50)
ns.generate_initial_pool()

from ase.db import connect
db = connect('clease_MnNiAs.db')