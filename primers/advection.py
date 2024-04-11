######################################################
# Example simulation with
# - advection of a single species
# - the species does not breed, evolve or die
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
# Setup the 10km x 10km grid with 1 x 8 spatial locations
# The units (km) are arbitrary at this stage, but
# later, the diffusion coefficient must be defined in the
# same units
g1 = Grid(0, 0, 10.0, 1, 8)

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
# Introduce 32 individuals at y = 0
# The carrying capacity is set at 100, but this
# is irrelevant because no lifecycle evolution is
# used in this example
all_pops.setPopulationAndParametersFromXY(0, 0, [32, 100])

######################################################
# put these populations into the spatial structure
spatial = SpatialDynamics(g1, all_pops)

######################################################
# Define the wind
# advection_wind.csv sets the advection velocity to
# (0, 80) km/day in all cells
# With the PDF = [[0.25, 1.0]] this means that 100% of
# individuals that advect will move 0.25*80 = 20km
# towards the positive y direction every timestep
wind = Wind("advection_wind.csv", "dummy.csv", [[0.25, 1.0]], g1)
wind.parseRawFile()

######################################################
# Evolve in time
dt = 1.0 # time step size, with units of days
t = 0
for i in range(5):
   # note the following is commented, so there is no lifecycle evolution
   # spatial.evolveCells(dt)
   # The following means that 50% of the population will
   # advect using the wind: the remaining will not advect
   spatial.advect(array.array('f', [0.5]), wind)
   t += dt

######################################################
# Output useful information
spatial.outputCSV("advection.csv", 0, "0", "")

######################################################
# Comparison with expected value
import matplotlib.pyplot as plt

yvals = list(range(0, 71, 10))
with open("advection.csv", "r") as f:
   f.readline()
   f.readline()
   data = [float(p.strip()) for p in f.readlines()]

plt.figure()
plt.bar(yvals, data, width = 10, label = 'Mozzie')
plt.scatter(yvals, [1, 0, 5, 0, 10, 0, 10, 0], label = 'Expected')
plt.legend()
plt.grid()
plt.xlabel("y")
plt.ylabel("population")
plt.suptitle("Population after advection from y = 0", fontsize = 16)
plt.title("25% of individuals advect each time step", fontsize = 12)
plt.savefig("advection.pdf", bbox_inches = 'tight')
plt.show()
plt.close()

sys.exit(0)
