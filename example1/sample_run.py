# Processes wind vectors given by Nick into a single file, wind.csv, and changing to m/day
import os
import sys
import time

# so we can find our ../code no matter how we are called
findbin = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(findbin + "/../code")

from cell import Cell
from wind import Wind
from grid import Grid
from populations import Populations
from spatial import Spatial

# setup the grid and the active/inactive information
sys.stdout.write("Initialising grid...")
start = time.clock()
g1 = Grid(-4614384.0,-3967418.0,5000.0,1517,1667)
sys.stdout.write(" " + str(time.clock() - start) + "s\n")

sys.stdout.write("Setting active and inactive...")
start = time.clock()
g1.setActiveAndInactive("active.csv")
sys.stdout.write(" " + str(time.clock() - start) + "s\n")
sys.stdout.write("There are " + str(g1.getNumActiveCells()) + " active cells\n")


if True:
   sys.stdout.write("Initialising wind...")
   start = time.clock()
   w = Wind(os.path.join(findbin, "raw_wind.csv"), os.path.join(findbin, "wind_8_subdivisions.csv"), [[0.05, 0.125], [0.1, 0.125], [0.15, 0.125], [0.2, 0.125], [0.25, 0.125], [0.3, 0.125], [0.4, 0.125], [0.5, 0.125]], g1)
   sys.stdout.write(" " + str(time.clock() - start) + "s\n")
   sys.stdout.write("Parsing raw_wind.csv and doing particle tracking (not necessary if particle-tracking pre-done)...")
   start = time.clock()
   w.parseRawFile()
   sys.stdout.write(" " + str(time.clock() - start) + "s\n")
   sys.stdout.write("Outputting particle-tracking result (not necessary if particle-tracking pre-done)...")
   start = time.clock()
   w.outputProcessedCSV()
   sys.stdout.write(" " + str(time.clock() - start) + "s\n")
   sys.stdout.write("Inputing the particle-tracking result...")
   start = time.clock()
   w.parseProcessedFile()
   sys.stdout.write(" " + str(time.clock() - start) + "s\n")

# define mosquito populations
sys.stdout.write("Defining the populations...")
start = time.clock()
all_pops = Populations(g1)
sys.stdout.write(" " + str(time.clock() - start) + "s\n")
# some place has a non-zero population
sys.stdout.write("Introducing mosquitoes...")
start = time.clock()
pop = list(range(Cell().getNumberOfPopulations()))
all_pops.setPopulation(g1.getNumActiveCells() // 2, pop)
sys.stdout.write(" " + str(time.clock() - start) + "s\n")

# initialise the spatial structure with timestep = 1 day and diffusion coefficient = 1E4 m^2/day
sys.stdout.write("Defining the spatial structure, ready for diffusion...")
start = time.clock()
spatial = Spatial(1.0, 1E4, g1, all_pops)
sys.stdout.write(" " + str(time.clock() - start) + "s\n")

sys.stdout.write("Diffusing\n")
start = time.clock()
for i in range(1, 11):
   sys.stdout.write("Time step " + str(i) + "\n")
   spatial.diffuse()
sys.stdout.write("Time for 1 diffusion step = " + str((time.clock() - start) / 10.0) + "s\n")


sys.exit(0)
