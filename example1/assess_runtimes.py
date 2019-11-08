# Processes wind vectors given by Nick into a single file, wind.csv, and changing to m/day
import os
import sys
import array

sys.stdout.write("NOTE: You will have to run this a few times to properly assess runtimes\n")
sys.stdout.write("      This is because your runtime will be impacted by other processes on your system\n\n")

import timeit

# so we can find our ../code no matter how we are called
findbin = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(findbin + "/../code")

from wind import Wind
from grid import Grid
from cellDynamics import CellDynamicsStatic15_9_3_2
from spatialDynamics import SpatialDynamics
from populationsAndParameters import PopulationsAndParameters

# setup the grid and the active/inactive information
sys.stdout.write("Initialising grid (without building adjacency list)...")
start = timeit.default_timer()
g1 = Grid(-4614384.0, -3967418.0, 5000.0, 1517, 1667, False)
sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")

sys.stdout.write("Setting active and inactive...")
start = timeit.default_timer()
g1.setActiveAndInactive("active.csv")
sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")
sys.stdout.write("There are " + str(g1.getNumActiveCells()) + " active cells\n")


# define mosquito populations and parameters
sys.stdout.write("Defining the populations...")
start = timeit.default_timer()
cell = CellDynamicsStatic15_9_3_2()
all_pops = PopulationsAndParameters(g1, cell)
sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")
# some place has a non-zero population
sys.stdout.write("Introducing mosquitoes...")
start = timeit.default_timer()
pop = list(range(all_pops.getCell().getNumberOfPopulations() + all_pops.getCell().getNumberOfParameters()))
all_pops.setPopulationAndParameters(g1.getNumActiveCells() // 2, pop)
sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")

# initialise the spatial structure with timestep = 1 day and diffusion coefficient = 1E4 m^2/day and advection fraction = 0.5
# No wind
sys.stdout.write("Defining the spatial structure, ready for diffusion...")
start = timeit.default_timer()
spatial = SpatialDynamics(1.0, 1E4, 0.5, g1, array.array('f', []), [], all_pops)
sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")

sys.stdout.write("Diffusing\n")
start = timeit.default_timer()
for i in range(1, 11):
   sys.stdout.write("Time step " + str(i) + "\n")
   spatial.diffuse()
sys.stdout.write("Time for 1 diffusion step = " + str((timeit.default_timer() - start) / 10.0) + "s\n")

if True:
   sys.stdout.write("Initialising wind...")
   start = timeit.default_timer()
   w = Wind(os.path.join(findbin, "raw_wind.csv"), os.path.join(findbin, "wind_8_subdivisions.csv"), [[0.05, 0.125], [0.1, 0.125], [0.15, 0.125], [0.2, 0.125], [0.25, 0.125], [0.3, 0.125], [0.4, 0.125], [0.5, 0.125]], g1)
   sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")
   particle_tracking_needed = True
   if particle_tracking_needed:
      sys.stdout.write("Parsing raw_wind.csv and doing particle tracking (not necessary if particle-tracking pre-done)...")
      start = timeit.default_timer()
      w.parseRawFile()
      sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")
   output_required = False
   if output_required:
      sys.stdout.write("Outputting particle-tracking result (not necessary if particle-tracking pre-done)...")
      start = timeit.default_timer()
      w.outputProcessedCSV()
      sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")
   input_required = False
   if input_required:
      sys.stdout.write("Inputing the particle-tracking result...")
      start = timeit.default_timer()
      w.parseProcessedFile()
      sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")



sys.exit(0)
