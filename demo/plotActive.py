import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

with open('active.csv', 'r') as f:
    header = f.readline().strip().split(",")
xmin, ymin, cell_size, nx, ny = [float(x.split("=")[1]) for x in header]
    
data = np.genfromtxt('active.csv', delimiter = ',', skip_header = True)
data = np.flipud(data) # because of upside-down convention used in mozzie

plt.figure()
custom = ListedColormap(['cyan', 'none'])
plt.imshow(data, extent = (0, nx * cell_size, 0, ny * cell_size), cmap = custom)
plt.text(45000, 0, "Water body\nshown in cyan", ha = 'right', va = 'bottom')
plt.xlabel("Easting + " + str(int(-xmin + 0.5)))
plt.ylabel("Northing + " + str(int(-ymin + 0.5)))
plt.title("Inactive cells shown in blue")
plt.savefig("active.pdf", bbox_inches = 'tight')
plt.show()
