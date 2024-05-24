######################################################
# Example simulation with
# - lifecycle, diffusion and advection of a single species
# The lifecycle is governed by the logistic equation, and
# wind velocity is spatio-temporal dependent, and
# carrying capacity is spatio-temporal dependent
#
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
from spatialDependence import SpatialDependence
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
# define the carrying capacity
# This is higher on the island than the mainland,
# and varies seasonally (first 100 days of year it
# is high, remainder of the year it is low).
cc_parser = SpatialDependence(0, 0, 10.0, 5, 3)
cc_parser.parse("spatial_cc_high.csv", "generic_float", [])
cc_parser.restrictToActive(g1.getGlobalIndex())
cc_high = cc_parser.getData0()
cc_parser = SpatialDependence(0, 0, 10.0, 5, 3)
cc_parser.parse("spatial_cc_low.csv", "generic_float", [])
cc_parser.restrictToActive(g1.getGlobalIndex())
cc_low = cc_parser.getData0()

######################################################
# Using CellDynamicsLogistic1_1() means that at each
# spatial location there is:
# - 1 population
# - 1 parameter, which is the carrying capacity
# Define the populations and parameters over the grid:
all_pops = PopulationsAndParameters(g1, cell)
# Introduce 50 individuals at (x, y) = (10, 0)
# Carrying capacity at this point is irrelevant here
# because it is overwritten later on
all_pops.setPopulationAndParametersFromXY(10, 0, [50, 100])

######################################################
# put these populations into the spatial structure
spatial = SpatialDynamics(g1, all_pops)

######################################################
# Define the wind
# advection_wind.csv sets the advection velocity to
# (80, 0) km/timestep on the mainland
# (0, 0) km/timestep on the island
# Note the timestep in the denominator (time step =
# 1 day, below, because there is 0.5 days of lifecycle
# evolution and 0.5 days of diffusion).  It is written
# in this way because the wind is applied every time step.
#
# With the PDF = [[0.1, 0.625], [0.2, 0.375]],
# then, of the individuals that advect:
# - 62.5% of them will move 0.1*80 = 8km
# - 37.5% of them will move 0.2*80 = 16km
# (because of the spatial discretisation to 10km,
# 8km -> 10km and 16km -> 20km)
# towards the positive x direction every timestep
wind = Wind("spatial_wind.csv", "dummy.csv", [[0.1, 0.625], [0.2, 0.375]], g1)
wind.parseRawFile()

######################################################
# Evolve in time
dt = 0.5 # time step size, with units of days
t = 0
for year in range(5):
   # The first 100 days are windy and have high carrying-capacity
   # Set the carrying capacity to the "high" values
   pop_and_param_array = all_pops.getQuantities()
   for i in range(g1.getNumActiveCells()):
      pop_and_param_array[2 * i + 1] = cc_high[i]
   for i in range(100):
      # The following may be interpreted as:
      # - allow lifecycle evolution to occur for 0.5 days
      # - then diffuse for a remaining 0.5 days
      # - then apply wind corresponding to the whole day
      # You will have to consider whether this interpretation
      # is appropriate for your particular species
      #
      # perform the lifecycle evolution (solve the logistic equation)
      spatial.evolveCells(dt)
      t += dt
      # diffuse with diffusion coefficient 0.1 km^2/day
      spatial.diffuse(dt, 0.1)
      # Finally, apply advection.
      # The following means that 0.1% of the population will
      # advect using the wind: the remaining will not advect
      spatial.advect(array.array('f', [0.001]), wind)
      t += dt
   # the next 265 days have no wind, smaller diffusion, and
   # low carrying capacity
   # Set the carrying capacity to the "low" values
   pop_and_param_array = all_pops.getQuantities()
   for i in range(g1.getNumActiveCells()):
      pop_and_param_array[2 * i + 1] = cc_low[i]
   for i in range(365 - 100): 
      spatial.evolveCells(dt)
      t += dt
      spatial.diffuse(dt, 0.01)
      t += dt

######################################################
# Output useful information
spatial.outputCSV("spatial.csv", 0, "0", "")

######################################################
# Plot results
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

xvals, yvals = np.meshgrid(np.linspace(0, 40, 5), np.linspace(0, 20, 3))
data = []
with open("spatial.csv", "r") as f:
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
fig.colorbar(c, ax = ax, shrink = 0.5)
plt.xlabel("x (km)")
plt.ylabel("y (km)")
plt.title("Population after spatio-temporal growth, diffusion and advection")
plt.gca().set_aspect('equal')
plt.xticks([0, 10, 20, 30, 40])
plt.yticks([0, 10, 20])
plt.savefig("spatial.pdf", bbox_inches = 'tight')
plt.show()
plt.close()

sys.exit(0)
