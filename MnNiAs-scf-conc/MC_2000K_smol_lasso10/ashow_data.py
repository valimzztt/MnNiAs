import os
import json
import matplotlib.pyplot as plt
import numpy as np

def calculate_heat_capacity(data):
    temperatures = list(data.keys())
    data = list(data.values())
    energies = [point[0] for point in data]
    num_points = len(temperatures)

    # Initialize arrays for storing central differences and heat capacity
    central_diffs = [0.0] * num_points
    heat_capacity = [0.0] * num_points

    # Calculate central differences for the interior points
    for i in range(1, num_points - 1):
        central_diffs[i] = (energies[i + 1] - energies[i - 1]) / (temperatures[i + 1] - temperatures[i - 1])

    # Calculate forward difference for the first point
    central_diffs[0] = (energies[1] - energies[0]) / (temperatures[1] - temperatures[0])

    # Calculate backward difference for the last point
    central_diffs[-1] = (energies[-1] - energies[-2]) / (temperatures[-1] - temperatures[-2])

    # Calculate heat capacity using central differences
    for i in range(num_points):
        heat_capacity[i] = central_diffs[i]
    return heat_capacity
import numpy as np

def richardson(f, x, n, h):
    # Initialize the Richardson extrapolation table
    d = np.zeros((n + 1, n + 1))
    # Compute the first column of the table using central difference method
    for i in range(n + 1):
        d[i, 0] = (f(x + h) - f(x - h)) / (2 * h)
        h /= 2

    # Perform Richardson extrapolation
    for j in range(1, n + 1):
        for i in range(j, n + 1):
            d[i, j] = d[i, j - 1] + (d[i, j - 1] - d[i - 1, j - 1]) / (4 ** j - 1)

    # Return the most accurate estimate
    return d[n, n]


# assign directory
cwd = os.getcwd()
import re
dir_hubb_mc = 'MC_2000K_smol_lasso10'
directory = os.path.join(cwd, dir_hubb_mc)
data_dict = {}
temperatures =[]
info = []
for filename in os.listdir(directory):
    f = os.path.join(directory, filename)
    if os.path.isfile(f):
        if(filename.startswith("MC") and not filename.endswith(".json")):
            with open(f) as f:
                read = f.read()
                print(read)
                
                pattern = r"'(\w+)': ([\-\d.]+)"
                matches = re.findall(pattern, read)

                mc_data = {}

                for match in matches:
                    key = match[0]
                    value = float(match[1])
                    mc_data [key] = value

                temperature = mc_data["temperature"]
                energy = mc_data["avg_energy"]
                heat_capacity = mc_data["heat_capacity"]
                if(temperature > 300): 
                    temperatures.append(temperature)
                    info.append((energy, heat_capacity))
                   
                
data = dict(zip(temperatures, info))

temps = list(data.keys())
temps.sort()
data= {i: data[i] for i in temps}
print(data)
temperatures = list(data.keys())
info = list(data.values())
energies = [point[0] for point in info]
max_energy = max(energies)
energies = [energy- max_energy for energy in energies]
heat_capacities = [point[1] for point in info]
hc_calc = calculate_heat_capacity(data)
# Degree of polynomial function used as fitting func
deg = 30
x_data = np.linspace(min(temperatures),max(temperatures),num=200)
mymodel = np.poly1d(np.polyfit(temperatures, energies, deg))

#Using Richardson Extrapolation to avoid oscillations: 
hc_rich = np.zeros(len(x_data))
n = 3  # Number of levels of extrapolation
h = 1 # Initial stepsize
for i in range(len(x_data)):
    x = x_data[i]
    approx_derivative = richardson(mymodel, x, n, h)
    hc_rich[i] = approx_derivative

#store energy for reference
json_filename = os.path.join(directory,"energies.json")
with open(json_filename, "w") as fp:
    json.dump(data,fp)
    
    
fig, ax = plt.subplots(2, 2)
plt.subplot(2, 1, 1) # row 1, col 2 index 1
plt.plot(temperatures, energies, "o",  marker='o',  label="measured")
plt.plot(x_data , mymodel(x_data ), label="Polyfit of degree: " + str(deg))
plt.grid(True)
plt.xlabel('Temperature (K)')
plt.ylabel('Energy (eV)')
plt.legend()
plt.subplot(2, 1, 2) # index 2
plt.plot(temperatures, heat_capacities, "bo",  marker='.',  label ="measured")
#plt.plot(temperatures, hc_calc, "ro",  marker='.', label ="Central difference approximation")

#plt.plot(x_data, hc_rich, label ="Richardson Extrap. order %d step size %f" % (n, h))
plt.xlabel('Temperature (K)')
plt.ylabel('Heat capacity (J*K)')
plt.grid(True)
plt.legend(loc=2, prop={'size': 6})
fig.tight_layout()
energy_plot = os.path.join(directory, "amc-result.png")
plt.savefig(energy_plot)
plt.close()
