######################################################
# Example simulation with
# - diffusion of a single species
# - the species does not breed, evolve or die
#
######################################################
# Import libraries
import os
import sys
import math

# so we can find our ../code no matter how we are called
findbin = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(findbin + "/../code")

from grid import Grid
from cellDynamics import CellDynamicsLogistic1_1
from spatialDynamics import SpatialDynamics
from populationsAndParameters import PopulationsAndParameters

######################################################
# Setup the 1km x 1km grid with 101 x 101 spatial locations
# The center of this grid is at (x, y) = (50km, 50km)
# The units (km) are arbitrary at this stage, but
# later, the diffusion coefficient must be defined in the
# same units
g1 = Grid(0, 0, 1.0, 101, 101)

######################################################
# define the lifecycle dynamics to be governed by the
# logistic equation.  Since the the species is not
# allowed to breed, evolve or die, this is largely
# irrelevant.  Logistic1_1 is just the simplest
# type of lifecycle that can be made into a "dummy"
# that does not breed, evolve or die.
cell = CellDynamicsLogistic1_1()

######################################################
# Using CellDynamicsLogistic1_1() means that at each
# spatial location there is:
# - 1 population
# - 1 parameter, which is the carrying capacity
# Define the populations and parameters over the grid:
all_pops = PopulationsAndParameters(g1, cell)
# Introduce 100 individuals at the center.
# The carrying capacity is set at 1000, but this
# is irrelevant because no lifecycle evolution is
# used in this example
all_pops.setPopulationAndParametersFromXY(50, 50, [1E6, 1000])

######################################################
# put these populations into the spatial structure
spatial = SpatialDynamics(g1, all_pops)

######################################################
# Evolve in time
dt = 0.5 # time step size, with units of days
t = 0
for i in range(201):
   # note the following is commented, so there is no lifecycle evolution
   # spatial.evolveCells(dt)
   # Use the diffusion coefficient 0.3 km^2/day
   # The units are the same as the Grid and dt
   spatial.diffuse(dt, 0.3)
   t += dt

######################################################
# Output useful information
spatial.outputCSV("diffusion_100days.csv", 0, "0", "")
spatial.outputCSVsubset(0, 50, 101, 51, "diffusion_100days_line.csv", 0, "0", "")

######################################################
# Contour and line plots
import numpy as np
import matplotlib.pyplot as plt

xvals, yvals = np.meshgrid(np.linspace(0, 100, 101), np.linspace(0, 100, 101))
data = []
with open("diffusion_100days.csv", "r") as f:
   f.readline()
   f.readline()
   for y in range(101):
      data.append(list(map(float, f.readline().strip().split(","))))
dmin = min([min(d) for d in data])
dmax = max([max(d) for d in data])
fig, ax = plt.subplots()
c = ax.pcolormesh(xvals, yvals, data, cmap = 'Greys', vmin = dmin, vmax = dmax)
fig.colorbar(c, ax = ax)
plt.grid()
plt.xlabel("x")
plt.ylabel("y")
plt.title("Population after diffusion from central position")
plt.savefig("diffusion.pdf", bbox_inches = 'tight')
plt.show()
plt.close()

with open("diffusion_100days_line.csv", "r") as f:
   f.readline()
   f.readline()
   data = list(map(float, f.readline().strip().split(",")))

plt.figure()
xvals = np.linspace(0, 100, 101)
plt.plot(xvals, data, linewidth = 3, label = "Mozzie")
plt.plot(xvals, 1E6 * np.exp(-np.power(xvals - 50, 2) / 4 / 0.3 / 100) / 4 / np.pi / 0.3 / 100, 'r--', label = "Analytic solution")
plt.grid()
plt.legend()
plt.xlabel("x")
plt.ylabel("Population")
plt.title("Population along y = 50, after diffusion from central position")
plt.savefig("diffusion_line.pdf", bbox_inches = 'tight')
plt.show()

sys.exit(0)
