# Allows assessment of runtimes for logistic growth only
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
from cellDynamics import CellDynamicsLogistic1_1
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
g1.setActiveAndInactive("active.csv")
sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")
sys.stdout.write("There are " + str(g1.getNumActiveCells()) + " active cells\n")


# read a file corresponding to carrying capacity
sys.stdout.write("Reading the carrying capacity and restricting to active cells...")
start = timeit.default_timer()
cc_parser = SpatialDependence(-4614.0, -3967.0, 5.0, 1517, 1667)
cc_parser.parse("carrying.csv", "generic_float", [])
cc_parser.restrictToActive(g1.getGlobalIndex())
cc = cc_parser.getData0()
sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")

# define a zeroed populations and parameters array
sys.stdout.write("Defining the populations and parameters array...")
start = timeit.default_timer()
cell = CellDynamicsLogistic1_1()
all_pops = PopulationsAndParameters(g1, cell)
sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")

# introduce carrying capacity (which, for CellDynamicsLogistic1_1, is second slot in the populations and params array), which is max(1.0, cc)
# and a population of 1 mosquito everywhere that's active
sys.stdout.write("Populating the initial (zero) populations and parameters array...")
start = timeit.default_timer()
pop_and_param_array = all_pops.getQuantities()
for i in range(g1.getNumActiveCells()):
   pop_and_param_array[2 * i] = 1.0
   pop_and_param_array[2 * i + 1] = max(1.0, cc[i]) # cannot have zero carrying capacity
sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")

# put these populations into the spatial structure
sys.stdout.write("Defining the spatial structure, ready for logistic growth...")
start = timeit.default_timer()
spatial = SpatialDynamics(g1, all_pops)
sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")

sys.stdout.write("Doing logistic growth\n")
start = timeit.default_timer()
for i in range(1, 11):
   sys.stdout.write("Time step " + str(i) + "\n")
   # diffuse with timestep = 100 day
   spatial.evolveCells(100)
sys.stdout.write("Time for 1 diffusion step = " + str((timeit.default_timer() - start) / 10.0) + "s\n")

if True:
   # output something that can be visualised
   sys.stdout.write("Outputting population result...")
   start = timeit.default_timer()
   spatial.outputCSV("logistic.csv", 0, "0", "")
   sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")
   



sys.exit(0)
