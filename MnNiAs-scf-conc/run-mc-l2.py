from clease.settings import Concentration
import os 

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
eci_file=open('MnNiAs_vaspl2.json', 'r')
eci = json.load(eci_file)

from clease.calculator import attach_calculator
atoms = settings.atoms.copy()*(6, 6, 6)
atoms = attach_calculator(settings, atoms=atoms, eci=eci)
print(atoms)

for i in range(0,len(atoms),1):
    if atoms[i].symbol=='Mn' and i%2 == 0:
        atoms[i].symbol='Ni'
print(atoms)

from clease.montecarlo import Montecarlo
from clease.montecarlo.observers import EnergyEvolution
# Monte Carlo results will be created in new directory
directory = 'MC_2000K_results_l2_clease2'
cwd = os.getcwd()
parent_dir = cwd
path = os.path.join(parent_dir, directory)
os.mkdir(path)  
os.chdir(path)
from clease.montecarlo import Montecarlo
T = 2000
nsteps=1000000
mc = Montecarlo(atoms, T)

from clease.montecarlo.constraints import FixedElement
cnst = FixedElement('As')
mc.generator.add_constraint(cnst)

from clease.montecarlo.observers import EnergyEvolution, Snapshot, CorrelationFunctionObserver, AcceptanceRate
obs = EnergyEvolution(mc)
mc.attach(obs, interval=1000)
corr = CorrelationFunctionObserver(atoms.calc)
mc.attach(corr)
accrate = AcceptanceRate()
mc.attach(accrate)
snap = Snapshot(atoms, fname='snapshot')
mc.attach(snap, interval=nsteps)

mc.run(steps=nsteps)
energies_eq = obs.energies
thermo_eq = mc.get_thermodynamic_quantities()

f = open('MC_eq_{}'.format(T),'w')
print(thermo_eq, file=f)
print(energies_eq, sep='\n', file=f)
f.close()

mc.run(steps=nsteps)
energies = obs.energies
thermo = mc.get_thermodynamic_quantities()

f = open('MC_{}'.format(T),'w')
print(thermo, file=f)
print(energies, sep='\n', file=f)
f.close()


for i in range(2000, 0, -5):
    T = i
    #atoms = vasp.read_vasp(POSCARstring)
    #atoms = attach_calculator(settings, atoms=atoms, eci=eci)
    mc = Montecarlo(atoms, T)
    mc.generator.add_constraint(cnst)
    obs = EnergyEvolution(mc)
    mc.attach(obs, interval=1000)
    snap = Snapshot(atoms, fname='snapshot')
    mc.attach(snap, interval=nsteps)
    corr = CorrelationFunctionObserver(atoms.calc)
    mc.attach(corr)
    accrate = AcceptanceRate()
    mc.attach(accrate)
    
    mc.run(steps=nsteps)
    mc.run(steps=nsteps)

    energies = obs.energies
    rate=accrate.get_averages()
    thermo = mc.get_thermodynamic_quantities()
    corr_dict=corr.get_averages()
    
    from ase.build import sort
    from ase.io import vasp, db
    sorted_atoms=sort(atoms)
    vasp.write_vasp('POSCAR{}.vasp'.format(i), sorted_atoms, direct=False, wrap=False)
    POSCARstring='POSCAR{}.vasp'.format(i)
    f = open('MC_{}'.format(i),'w')
    print(thermo, file=f)
    print(energies, sep='\n', file=f)
    f.close()
os.chdir(cwd)


