import os
import sys
import numpy as np
import matplotlib.pyplot as plt

with open('activeDemo.csv', 'r') as f:
    header = f.readline().strip().split(",")
xmin, ymin, cell_size, nx, ny = [float(x.split("=")[1]) for x in header]
    
data = np.genfromtxt('activeDemo.csv', delimiter = ',', skip_header = True)
data = np.flipud(data) # because of upside-down convention used in mozzie

plt.figure()
plt.imshow(data, extent = (0, nx * cell_size, 0, ny * cell_size), cmap = 'terrain')
plt.xlabel("Easting + " + str(int(-xmin + 0.5)))
plt.ylabel("Northing + " + str(int(-ymin + 0.5)))
plt.title("Inactive cells shown in blue")
plt.savefig("activeDemo.pdf", bbox_inches = 'tight')
plt.show()
