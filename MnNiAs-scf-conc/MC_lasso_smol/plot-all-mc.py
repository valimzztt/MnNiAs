import os
import json 
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import cm
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 24})

# Get the current directory
current_directory =  os.getcwd()

init_temperature = 2000
# Loop through the current directory
rounds = 10

import seaborn as sns
mc_dir = os.path.join(current_directory, "MC_lasso_smol")
colors = sns.color_palette('Set1', n_colors=rounds)
# here we are creating sub plots

fig, ax = plt.subplots(figsize=(10, 8))  

i = 0
for file in os.listdir(mc_dir):
    filename = os.path.join(mc_dir, file)
    print(filename)
    data_dict = {}
    temperatures =[]
    info = []
    round = file[-6]
    if(filename.endswith("json")):
        i = i +1 
        with open(filename) as f:
            read = f.read()
            # Convert the string to a JSON object
            mc_data = json.loads(read)
            temperatures = mc_data.keys()
            temperatures = [float(temp) for temp in temperatures]
            info = list(mc_data.values())
            energies = [point[0] for point in info]
            heat_capacity =[point[1] for point in info]
            ax.plot(temperatures, energies, color=colors[i],  marker='.',  label="run " + str(round),linewidth=5)
  
plt.xlabel('Temperature (K)')
plt.ylabel('Energy (eV))')
plt.legend()
hc_plot = os.path.join(mc_dir, "energy-mc-result.png")
plt.savefig(hc_plot)
plt.close()

fig2, ax2 = plt.subplots(figsize=(12, 10))  
j = 0
for file in os.listdir(mc_dir):
    filename = os.path.join(mc_dir, file)
    print(filename)
    data_dict = {}
    temperatures =[]
    info = []
    round = file[-6]
    if(filename.endswith("json")):
        j = j + 1
        with open(filename) as f:
            read = f.read()
            # Convert the string to a JSON object
            mc_data = json.loads(read)
            temperatures = mc_data.keys()
            temperatures = [float(temp) for temp in temperatures]
            info = list(mc_data.values())
            energies = [point[0] for point in info]
            heat_capacity =[point[1] for point in info]
            ax2.plot(temperatures, heat_capacity, "bo",  color=colors[j],  marker='.',  label="run: " + str(round),linewidth=5)

#plt.plot(x_data, hc_rich, label ="Richardson Extrap. order %d step size %f" % (n, h))
plt.xlabel('Temperature (K)')
plt.ylabel('Heat capacity (J*K)')
plt.legend()
# fig.tight_layout()
hc_plot = os.path.join(mc_dir, "heatcapacity-mc-result.png")
plt.savefig(hc_plot)
