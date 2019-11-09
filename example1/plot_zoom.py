#!/usr/bin/env python

# plots the active cells as dark-grey areas and the inactive cells as light-grey areas
import sys
import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) != 3:
    sys.stderr.write("You must supply 2 arguments to this script: the CSV file containing the active/inactive information; the output filename\n")
    sys.stderr.write("Eg " + sys.argv[0] + " data.csv picture.png\n")
    sys.exit(1)

with open(sys.argv[1], 'r') as f:
    data = [line for line in f.readlines() if line.strip()]

for line in data:
    if line.startswith("#xmin"):
        line = line.split(",")
        xmin = float(line[0].split("=")[1])
        ymin = float(line[1].split("=")[1])
        cell_size = float(line[2].split("=")[1])
        nx = int(line[3].split("=")[1])
        ny = int(line[4].split("=")[1])
        break

z = [[float(s) for s in line.strip().split(",")] for line in data if not line.startswith("#")]

x = np.arange(xmin, xmin + nx * cell_size, cell_size)
y = np.arange(ymin, ymin + ny * cell_size, cell_size)

plt.figure()
# for carrying_zoom: plt.contourf(x, y, z, levels = [3, 3.25, 3.5, 3.75, 4, 4.25, 4.5, 4.75, 5, 5.25, 5.5], vmin=3, vmax=5.5, cmap='coolwarm')
plt.contourf(x, y, z, cmap='coolwarm')
plt.axis('equal')
plt.colorbar()
plt.xlabel("x")
plt.ylabel("y")
plt.gca().set_xlim([-2100, -1700])
plt.gca().set_ylim([600, 1000])
plt.title(sys.argv[1])
plt.savefig(sys.argv[2])

sys.exit(0)
