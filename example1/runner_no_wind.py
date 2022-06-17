######################################################
# Example simulation with
# - 5km x 5km grid defined that covers africa entirely (nx=1517, ny=1667)
# - active cells defined simply by 'ocean=inactive' 'land=active', resulting in 1.4 million cells
# - diffusing and advecting timestep of 0.5days, then lifecycle timestep of 0.5days, repeated cyclically.  This means that logistic growth only occurs for 0.5 days out of each day, diffusion only occurs for 0.5 days out of each day, etc
# - diffusion (diffusion coefficient = 0.1 km^2/day)
# - In this example, there is no advection by wind
# - logistic growth with growth rate = 0.1/day, and carrying capacity defined by a file


######################################################
# Import libraries
import os
import sys
import array
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


######################################################
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


######################################################
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
sys.stdout.write("Populating the initial (zero) populations and parameters array...")
start = timeit.default_timer()
pop_and_param_array = all_pops.getQuantities()
for i in range(g1.getNumActiveCells()):
   pop_and_param_array[2 * i + 1] = max(1.0, cc[i]) # cannot have zero carrying capacity
# introduce 10000 mosquito at (x, y) = (-2000, 700), and set the carrying cap=150000 at that cell
all_pops.setPopulationAndParametersFromXY(-2000, 700, [10000.0, 150000.0])
sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")


######################################################
# put these populations into the spatial structure
sys.stdout.write("Defining the spatial structure, ready for diffusion, advection and logistic growth...")
start = timeit.default_timer()
spatial = SpatialDynamics(g1, all_pops)
sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")


######################################################
# define wind
sys.stdout.write("Initialising wind...")
start = timeit.default_timer()
import math
# this is the PDF i use
xa = [0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.17, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5]
ya = [math.exp(-x * 6) for x in xa]
su = sum(ya)
ya = [y / sum(ya) for y in ya]
pdf = [[xa[i], ya[i]] for i in range(len(xa))]
wind = Wind(os.path.join(findbin, "raw_wind.csv"), os.path.join(findbin, "wind_8_subdivisions.csv"), pdf, g1)
sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")
sys.stdout.write("Parsing raw_wind.csv and doing particle tracking (not necessary if particle-tracking pre-done)...")
start = timeit.default_timer()
wind.parseRawFile()
sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")


######################################################
sys.stdout.write("Doing logistic growth, diffusion and advection\n")
start = timeit.default_timer()
timestep_size = 0.5
for i in range(501):
   sys.stdout.write("Time step " + str(i) + "\n")
   # logistically grow
   spatial.evolveCells(0.5 * 10) # the multiplication here is just because Logistic1_1 has a hard-coded 0.01 growth rate which is too tiny
   # diffuse with diffusion coefficient 0.1 km^2/day
   spatial.diffuse(timestep_size, 0.1)
   # NO ADVECTION
   # spatial.advect(array.array('f', [0.01]), wind)
   if (i%100 == 0):
      sys.stdout.write("Outputting population result...\n")
      spatial.outputCSV("runner_no_wind_" + str(i) + "_days.csv", 0, "0", "")

sys.stdout.write("CPU time taken = " + str(timeit.default_timer() - start) + "s\n")

######################################################
# output populations
   



sys.exit(0)
