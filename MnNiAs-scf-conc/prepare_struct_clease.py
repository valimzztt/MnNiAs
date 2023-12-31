from clease.settings import Concentration
import os 
from ase.io import vasp, db
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
    db_name="clease_MnNiAs-conc.db",
    basis_func_type='binary_linear',
    max_cluster_dia=(7,7,7))


import json
eci_file=open('MnNiAs_vasp.json', 'r')
eci = json.load(eci_file)

from clease.calculator import attach_calculator
# this means that we are changing the atoms manager to a larger DEFAULT VALUE! 
atoms = settings.atoms.copy()*(6,6,6)
i = 1
vasp.write_vasp('before_calc.vasp', atoms, direct=False, wrap=False)


atoms = attach_calculator(settings, atoms=atoms, eci=eci)
vasp.write_vasp('after_calc.vasp', atoms, direct=False, wrap=False)
