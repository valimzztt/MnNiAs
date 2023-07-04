
from clease.settings import Concentration
import os 

cwd = os.getcwd()
curr_directory = os.path.join(cwd, "MnNiAs-clease")
database = os.path.join(curr_directory , "clease_MnNiAs_energy.db")

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
    db_name=database,
    basis_func_type='binary_linear',
    max_cluster_dia=(7,6,6))

from clease import Evaluate
eva = Evaluate(settings=settings, scoring_scheme='loocv', nsplits=10)
# scan different values of alpha and return the value of alpha that yields
# the lowest CV score
eva.set_fitting_scheme(fitting_scheme='l2')
alpha = eva.plot_CV(alpha_min=1E-6, alpha_max=100.0, num_alpha=50)
eva.set_fitting_scheme(fitting_scheme='l2', alpha=alpha)
print("The chosen value of alpha is")
print(alpha)
eva.fit()

# plot ECI values
import clease.plot_post_process as pp
import json 
import matplotlib.pyplot as plt
import json
fig = pp.plot_fit(eva)
alpha_image = os.path.join(curr_directory,'alpha-fit-l2.png')
fig.savefig(alpha_image)

# plot ECI values
fig = pp.plot_eci(eva)
eci_image = os.path.join(curr_directory,'ECI-values-l2.png')
fig.savefig(eci_image)

# save a dictionary containing cluster names and their ECIs
cluster_expansion = curr_directory + '/MnNiAs_vasp'
eva.save_eci(fname=cluster_expansion)
filename = os.path.join(curr_directory, 'MnNiAs_vasp.json')
eci_file=open(filename, 'r')
eci = json.load(eci_file)
print(eci)

