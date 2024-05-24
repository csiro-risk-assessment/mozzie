######################################################
# Example simulation with
# - lifecycle, diffusion and advection of a single species
# The lifecycle is governed by the logistic equation.
# An island is separated from the mainland.
# A species initially exists only on the mainland,
# but advection carries it to the island.
#
######################################################
# Import libraries
import os
import sys
import array

# so we can find our ../code no matter how we are called
findbin = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(findbin + "/../code")

from grid import Grid
from cellDynamics import CellDynamicsLogistic1_1
from spatialDynamics import SpatialDynamics
from wind import Wind
from populationsAndParameters import PopulationsAndParameters

######################################################
# Setup the 10km x 10km grid with 5 x 3 spatial locations
# The units (km) are arbitrary at this stage, but
# later, the diffusion coefficient and wind velocity
# must be defined in the same units
g1 = Grid(0, 0, 10.0, 5, 3)

######################################################
# Set the active and inactive cells.
# Pictorially, these are:
# m m m * *
# m m m * i
# m m m * *
# where
# - "m" denotes the mainland
# - "*" an inactive cell (the ocean)
# - "i" denotes the island
# Note that the code always assumes that a "instant-kill zone"
# always exists around the grid, so the situation really
# looks like:
# + + + + + + +
# + m m m * * +
# + m m m * i +
# + m m m * * +
# + + + + + + +
# where "+" indicates a cell that instantly kills all
# species if they advect or diffuse there.
g1.setActiveAndInactive("island_active_inactive.csv")

######################################################
# define the lifecycle dynamics to be governed by the
# logistic equation.
cell = CellDynamicsLogistic1_1()

######################################################
# Using CellDynamicsLogistic1_1() means that at each
# spatial location there is:
# - 1 population
# - 1 parameter, which is the carrying capacity
# Define the populations and parameters over the grid:
all_pops = PopulationsAndParameters(g1, cell)
# The carrying capacity is set at 100
pop_and_param_array = all_pops.getQuantities()
for i in range(g1.getNumActiveCells()):
   pop_and_param_array[2 * i + 1] = 100
# Introduce 50 individuals at (x, y) = (10, 0)
all_pops.setPopulationAndParametersFromXY(10, 0, [50, 100])

######################################################
# put these populations into the spatial structure
spatial = SpatialDynamics(g1, all_pops)

######################################################
# Define the wind
# advection_wind.csv sets the advection velocity to
# (80, 0) km/day in all cells
# With the PDF = [[0.25, 1.0]] this means that 100% of
# individuals that advect will move 0.25*80 = 20km
# (which happens to be the distance to the island)
# towards the positive x direction every timestep
wind = Wind("island_wind.csv", "dummy.csv", [[0.25, 1.0]], g1)
wind.parseRawFile()

######################################################
# Evolve in time
dt = 1.0 # time step size, with units of days
t = 0
for i in range(730): # 2 years
   # perform the lifecycle evolution (solve the logistic equation)
   spatial.evolveCells(dt)
   # diffuse with diffusion coefficient 0.1 km^2/day
   spatial.diffuse(dt, 0.1)
   # The following means that 0.1% of the population will
   # advect using the wind: the remaining will not advect
   spatial.advect(array.array('f', [0.001]), wind)
   t += dt

######################################################
# Output useful information
spatial.outputCSV("island.csv", 0, "0", "")

######################################################
# Plot results
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

xvals, yvals = np.meshgrid(np.linspace(0, 40, 5), np.linspace(0, 20, 3))
data = []
with open("island.csv", "r") as f:
   f.readline()
   f.readline()
   for y in range(3):
      data.append(list(map(float, f.readline().strip().split(","))))
fig, ax = plt.subplots()
# blank out the surrounding "ocean"
ax.add_patch(patches.Rectangle((-10, -10), 60, 40, linewidth = 0, facecolor = 'black', alpha = 0.1))
c = ax.pcolormesh(xvals, yvals, data, cmap = 'viridis', vmin = 40, vmax = 80, shading = 'auto')
# blank out the intervening "ocean"
ax.add_patch(patches.Rectangle((25, -5), 10, 30, linewidth = 0, facecolor = 'white'))
ax.add_patch(patches.Rectangle((25, -5), 10, 30, linewidth = 0, facecolor = 'black', alpha = 0.1))
ax.add_patch(patches.Rectangle((35, -5), 10, 10, linewidth = 0, facecolor = 'white'))
ax.add_patch(patches.Rectangle((35, -5), 10, 10, linewidth = 0, facecolor = 'black', alpha = 0.1))
ax.add_patch(patches.Rectangle((35, 15), 10, 10, linewidth = 0, facecolor = 'white'))
ax.add_patch(patches.Rectangle((35, 15), 10, 10, linewidth = 0, facecolor = 'black', alpha = 0.1))
plt.text(25, -5, "mainland", horizontalalignment = 'right', verticalalignment = 'bottom')
plt.text(45, 6, "island", bbox=dict(facecolor='white', alpha=0.6, linewidth=0), horizontalalignment = 'right', verticalalignment = 'bottom')
plt.text(30, 0, "ocean", horizontalalignment = 'left', verticalalignment = 'bottom')
fig.colorbar(c, ax = ax)
plt.xlabel("x (km)")
plt.ylabel("y (km)")
plt.title("Population after growth, diffusion and advection")
plt.gca().set_aspect('equal')
plt.xticks([0, 10, 20, 30, 40])
plt.yticks([0, 10, 20])
plt.savefig("island.pdf", bbox_inches = 'tight')
plt.show()
plt.close()

sys.exit(0)
