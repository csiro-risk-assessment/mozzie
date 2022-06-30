# Allows assessment of runtimes for Mosquito23 growth only, with only death, emergence and aging, using the solve_ivp time integration
import os
import sys

sys.stdout.write("NOTE: You will have to run this a few times to properly assess runtimes\n")
sys.stdout.write("      This is because your runtime will be impacted by other processes on your system\n\n")

import timeit

# so we can find our ../code no matter how we are called
findbin = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(findbin + "/../code")

from wind import Wind
from grid import Grid
from cellDynamics import CellDynamicsMosquito23
from spatialDynamics import SpatialDynamics
from spatialDependence import SpatialDependence
from populationsAndParameters import PopulationsAndParameters

# setup the grid and the active/inactive information
sys.stdout.write("Initialising grid (without building adjacency list)...")
start = timeit.default_timer()
g1 = Grid(-4614.0, -3967.0, 5.0, 1517, 1667, False)
sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")

sys.stdout.write("Setting active and inactive...")
start = timeit.default_timer()
g1.setActiveAndInactive("../example1/active.csv")
sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")
sys.stdout.write("There are " + str(g1.getNumActiveCells()) + " active cells\n")


# read a file corresponding to carrying capacity
sys.stdout.write("Reading the carrying capacity and restricting to active cells...")
start = timeit.default_timer()
cc_parser = SpatialDependence(-4614.0, -3967.0, 5.0, 1517, 1667)
cc_parser.parse("../example1/carrying.csv", "generic_float", [])
cc_parser.restrictToActive(g1.getGlobalIndex())
cc = cc_parser.getData0()
sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")

# define a zeroed populations and parameters array
sys.stdout.write("Defining the populations and parameters array...")
start = timeit.default_timer()
cell = CellDynamicsMosquito23()
cell.setMuLarvae(0.1)
cell.setMuAdult(0.2)
cell.setFecundity(0.0)
cell.setAgingRate(0.3)
cell.setNumAges(2)
cell.setTimeIntegrationMethod("solve_ivp")
all_pops = PopulationsAndParameters(g1, cell)
sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")

# introduce carrying capacity (which, for CellDynamicsMosquito23, is 12th slot in the populations and params array)
# and a population of carrying-capacity for every population (just for ease of coding - this isn't realistic!)
sys.stdout.write("Populating the initial populations and parameters array...")
start = timeit.default_timer()
for i in range(13):
   all_pops.setOverActiveGrid(i, cc)
sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")

# put these populations into the spatial structure
sys.stdout.write("Defining the spatial structure, ready for ODE growth...")
start = timeit.default_timer()
spatial = SpatialDynamics(g1, all_pops)
sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")

sys.stdout.write("Doing ODE growth\n")
start = timeit.default_timer()
for i in range(1, 2):
   sys.stdout.write("Time step " + str(i) + "\n")
   spatial.evolveCells(1E-6)
sys.stdout.write("Time for 1 evolve step = " + str((timeit.default_timer() - start) / 1.0) + "s\n")

sys.exit(0)
