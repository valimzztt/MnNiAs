import os
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.outputs import Outcar, Vasprun
from ase.io.vasp import read_vasp_out
import json
from pymatgen.entries.computed_entries import ComputedStructureEntry

"""
    Script that retrieves the relaxed structure and the final computed energy from the VASP runs, converts
    it to a Pymatgen structures and saves each DFT Vasp run as a ComputedStructureEntry in the form of a dictionary
    that is then converted into a JSON file
"""
cwd = os.getcwd()
src = os.path.join(cwd, "MnNiAs-clease")
dest = os.path.join(cwd, "MnNiAs-smol")
dest =  os.path.join(dest, "dft-data")
energies = []
output_json = os.path.join(dest, "comp-struct-energy.json")

for foldername in os.listdir(src):
    if(foldername.startswith("vasp")):
        try: 
            folder = os.path.join(src, foldername)
            poscar_file = os.path.join( folder , "POSCAR")
            contcar_file = os.path.join( folder , "CONTCAR")
            outcar_file = os.path.join( folder, "OUTCAR")
            vasprun_file = os.path.join(folder, "vasprun.xml")
            vasprun  = Vasprun(filename=vasprun_file,
                            exception_on_bad_xml=True)
            if(vasprun.converged):
                struct = vasprun.final_structure
                poscar = Poscar.from_file(poscar_file)
                initial_struct = poscar.structure
                contcar = Poscar.from_file(contcar_file)
                initial_struct = poscar.structure
                outcar = Outcar(outcar_file)
                # We can both get the relaxed structure from the VASPRUN as well as the intial structure: 
                # how much will the ECI values change?
                # comp_entry = ComputedStructureEntry(struct, outcar.final_energy)
                comp_entry = ComputedStructureEntry(contcar.structure, outcar.final_energy)
                comp_entry_dict = comp_entry.as_dict()
                energies.append(comp_entry_dict)
        except:
            print(foldername)
            continue

with open(output_json, 'w') as fout:
    json.dump(energies, fout)