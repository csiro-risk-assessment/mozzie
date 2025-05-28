######################################################
# Read all daily abundances for all species, sex and 
# phenotypes, between 2022-01-01 and 2023-12-31.
# Plot a map of the first day where the phenotypes 
# with the construct ("wc", "wr", "cc", "cr", "rr")
# of any species ("Ac", "Ag") of any sex ("male", "female")
# reach each cell
######################################################

import numpy as np
import csv
import os
import sys
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

sys.stdout.write("Initialising paths and output directory\n")
######################################################
# Paths and filenames
#
# Set the working directory,
# which contains files for defining active cells, etc, as well as the output
# Here, we assume your working directory is mozzie/example_full
working_dir = os.path.dirname(os.path.realpath(sys.argv[0]))
os.chdir(os.path.expanduser(working_dir))
# Directory, relative to working_directory, where the daily abundances files are
# These data will have been created by runner.py
reading_dir = "output"
reading_dir = os.path.join(working_dir, reading_dir)

# Different categories of populations (total 24 populations)
# species
species = ["Ac", "Ag"]
# genotype
genotypes = ["ww", "wc", "wr", "cc", "cr", "rr"]
# sex
sexes = ["male", "female"]
# years simulated
years = list(range(2022, 2024))
# days in each month
monthdays = [31,28,31,30,31,30,31,31,30,31,30,31]

# multidimensional array to store 100x100 matrix of daily abundances 
# for all populations on two years
mozzies = np.zeros((2,6,2,2,12,31,100,100), dtype=float) 

# read the csv abundances files
sys.stdout.write("Reading abundances files:\n")
for s in range(2):
    S = species[s]
    sys.stdout.write(S + "\n")
    for g in range(6):
        G = genotypes[g]
        sys.stdout.write("  " + G + "\n")
        for e in range(2):
           E = sexes[e]
           sys.stdout.write("    " + E + "\n")
           for y in range(2):
              Y = years[y]
              sys.stdout.write("      " + str(Y) + "\n")
              for m in range(12):
                 mo = m + 1
                 for d in range(monthdays[m]):
                    day = d + 1
                    with open(f"{reading_dir}/{S}_{G}_{E}_{Y}_{mo:02d}_{day:02d}.csv", newline='') as file:
                       reader = csv.reader(file)
                       data = list(reader)[2:]  # Skip the first two rows
                       data = np.array(data, dtype=float)  # Convert list of strings to float array
                       data = data.reshape(100, 100) 
                       mozzies[s,g,e,y,m,d,:,:] = data

# keep only genotypes with construct (c or r)
mozzies = np.delete(mozzies, 0, axis=1)
# sum of populations with c or r for all sex and species
gm = np.sum(mozzies, axis=(0, 1, 2))
# 3D map with daily abundances
gm = gm.reshape(2*12*31, 100, 100)
# transform in 2D map where each element corresponds the first day where the
# abundance is more than one
def safe_min(vec):
    return np.nan if vec[0].size == 0 else np.min(vec)
gm = np.apply_along_axis(lambda x: safe_min(np.where(x > 1)), axis=0, arr=gm)
gm = np.flipud(gm) - sum(monthdays[:5]) # release was on 1 June 2022


# Get terrain colormap
terrain_cmap = plt.get_cmap("terrain")
terrain = ListedColormap(terrain_cmap(np.linspace(0, 1, 255))[::-1])
# Plot the map (roughly in months)
plt.imshow(gm * 12 / sum(monthdays), cmap=terrain, extent=(-25, 25, -25, 25), vmin = 0, vmax = 18)
plt.colorbar(label = "Months")
plt.title("Time taken for gene drive to invade domain")
# Add contour lines (in km)
# define centres of grid cells
x = np.arange(-24.75, 25, 0.5)
y = np.arange(24.75, -25, -0.5)
X, Y = np.meshgrid(x, y)
plt.contour(X, Y, gm, levels=10, colors="black", linewidths = 0.5)
plt.xlabel("Distance east of release point (km)")
plt.ylabel("Distance north of release point (km)")
# plot waterbody
wb = np.genfromtxt('active.csv', delimiter = ',', skip_header = True)
wb = np.flipud(wb) # because of upside-down convention used in mozzie
custom = ListedColormap(['cyan', 'none'])
plt.imshow(wb, extent=(-25, 25, -25, 25), cmap = custom)
plt.text(23, -25, "Water body\nshown in cyan", ha = 'right', va = 'bottom')
# Save the figure in png format
plt.savefig(os.path.join(working_dir, "demo.pdf"), bbox_inches = 'tight')
plt.close()

