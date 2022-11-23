######################################################
# Example simulation with
# - 5km x 5km grid defined that covers africa entirely (nx=1517, ny=1667)
# - active cells defined simply by 'ocean=inactive' 'land=active', resulting in 1.4 million cells (read from ../example1)
# - diffusing and advecting timestep of 0.5days, then lifecycle timestep of 0.5days, repeated cyclically.  This means that mosquito population evolution via the ODE only occurs for 0.5 days out of each day, diffusion only occurs for 0.5 days out of each day, etc
# - diffusion (diffusion coefficient = 0.1 km^2/day)
# - advection by wind that is defined by 1 wind file only, and has exponentially-distributed 'drop off' timesx (1% of the population of each cell is advected) (read from ../example1)
# - mosquito lifecycle defined by the ODEs defined in Eqns(4.1) and (4.2) of Beeton, Hosack, Wilkins, Forbes, Ickowicz and Hayes, Journal of Theoretical Biology 2019, parameters defined in Fig5(c) of that article, with carrying capacity defined as above, and cycling between that value and 0.2*above



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
from cellDynamics import CellDynamicsBeeton2_2
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
g1.setActiveAndInactive(os.path.join(findbin, "../example1/active.csv"))
sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")
sys.stdout.write("There are " + str(g1.getNumActiveCells()) + " active cells\n")


######################################################
# read a file corresponding to the "high" carrying capacity, and create the "low" carrying capacity from it
sys.stdout.write("Reading the carrying capacity and restricting to active cells...")
start = timeit.default_timer()
cc_parser = SpatialDependence(-4614.0, -3967.0, 5.0, 1517, 1667)
cc_parser.parse(os.path.join(findbin, "../example1/carrying.bin"), "generic_float_binary", [])
cc_parser.restrictToActive(g1.getGlobalIndex())
cc_high = cc_parser.getData0()
# some of carrying.csv is zero, but we divide by Kx and Ky in the ODEs, so set a minimum carrying capacity
cc_high = array.array('f', [max(1.0, c) for c in cc_high])
cc_low = array.array('f', [c * 0.2 for c in cc_high])
sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")

# define a zeroed populations and parameters array
sys.stdout.write("Defining the populations and parameters array...")
start = timeit.default_timer()
cell = CellDynamicsBeeton2_2()
all_pops = PopulationsAndParameters(g1, cell)
sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")

# introduce species x and y, and carrying capacity, equal for both populations
sys.stdout.write("Populating the initial (zero) populations and parameters array...")
start = timeit.default_timer()
initial_y = array.array('f', [c * 0.1 for c in cc_high])
all_pops.setPopulationAndParametersFromXY(-2000.0, 700.0, [10000, 0, 0, 0]) # introduce 10000 type X at (-2000,700) (the population of Y and carrying capacities get overwritten by the following lines)
all_pops.setOverActiveGrid(1, initial_y)
all_pops.setOverActiveGrid(2, cc_high)
all_pops.setOverActiveGrid(3, cc_high)
sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")


######################################################
# put these populations into the spatial structure
sys.stdout.write("Defining the spatial structure, ready for diffusion, advection and ODE evolution...")
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
wind = Wind(os.path.join(findbin, "../example1/raw_wind.bin"), os.path.join(findbin, "../example2/wind_8_subdivisions.csv"), pdf, g1)
wind.setBinaryFileFormat(1)
sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")
sys.stdout.write("Parsing raw_wind.csv and doing particle tracking (not necessary if particle-tracking pre-done)...")
start = timeit.default_timer()
wind.parseRawFile()
sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")


######################################################
# Simulate, with Kx varying between cc_high and cc_low, and Ky = cc_low always
sys.stdout.write("Doing ODE evolution, diffusion and advection, and time-varying Kx\n")
start = timeit.default_timer()
timestep_size = 0.5
the_time = 0.0
end_time = 200.0
for bigP in [20, 1]:
   all_pops.setPopulationAndParametersFromXY(-2000.0, 700.0, [10000, 0, 0, 0]) # introduce 10000 type X at (-2000,700) (the population of Y and carrying capacities get overwritten by the following lines)
   num_cycles = int(end_time / bigP + 0.5)
   all_pops.setOverActiveGrid(1, initial_y)
   for cycle_num in range(num_cycles + 1):
      sys.stdout.write("  Cycle " + str(cycle_num) + " of " + str(num_cycles + 1) + "\n")
      all_pops.setOverActiveGrid(2, cc_high)
      all_pops.setOverActiveGrid(3, cc_low)
      for sub_cycle in range(bigP):
         sys.stdout.write("    K_high subcycle " + str(sub_cycle) + "\n")
         spatial.evolveCells(timestep_size)
         spatial.diffuse(timestep_size, 0.1)
         spatial.advect(array.array('f', [0.01]), wind)
         the_time += timestep_size
      all_pops.setOverActiveGrid(2, cc_low)
      all_pops.setOverActiveGrid(3, cc_low)
      for sub_cycle in range(bigP):
         sys.stdout.write("    K_low subcycle " + str(sub_cycle) + "\n")
         spatial.evolveCells(timestep_size)
         spatial.diffuse(timestep_size, 0.1)
         spatial.advect(array.array('f', [0.01]), wind)
         the_time += timestep_size

   sys.stdout.write("Outputting population result...\n")
   spatial.outputCSV("runner_Kxvary_x_P_" + str(bigP) + "_" + str(cycle_num * bigP) + "_days.csv", 0, "0", "")
   spatial.outputCSV("runner_Kxvary_y_P_" + str(bigP) + "_" + str(cycle_num * bigP) + "_days.csv", 1, "0", "")


######################################################
# Simulate again, with both Kx and Ky varying in time
sys.stdout.write("Doing ODE evolution, diffusion and advection\n")
start = timeit.default_timer()
timestep_size = 0.5
the_time = 0.0
end_time = 200.0
for bigP in [20, 1]:
   all_pops.setPopulationAndParametersFromXY(-2000.0, 700.0, [10000, 0, 0, 0]) # introduce 10000 type X at (-2000,700) (the population of Y and carrying capacities get overwritten by the following lines)
   num_cycles = int(end_time / bigP + 0.5)
   all_pops.setOverActiveGrid(1, initial_y)
   for cycle_num in range(num_cycles + 1):
      sys.stdout.write("  Cycle " + str(cycle_num) + " of " + str(num_cycles + 1) + "\n")
      all_pops.setOverActiveGrid(2, cc_high)
      all_pops.setOverActiveGrid(3, cc_high)
      for sub_cycle in range(bigP):
         sys.stdout.write("    K_high subcycle " + str(sub_cycle) + "\n")
         spatial.evolveCells(timestep_size)
         spatial.diffuse(timestep_size, 0.1)
         spatial.advect(array.array('f', [0.01]), wind)
         the_time += timestep_size
      all_pops.setOverActiveGrid(2, cc_low)
      all_pops.setOverActiveGrid(3, cc_low)
      for sub_cycle in range(bigP):
         sys.stdout.write("    K_low subcycle " + str(sub_cycle) + "\n")
         spatial.evolveCells(timestep_size)
         spatial.diffuse(timestep_size, 0.1)
         spatial.advect(array.array('f', [0.01]), wind)
         the_time += timestep_size

   sys.stdout.write("Outputting population result...\n")
   spatial.outputCSV("runner_x_P_" + str(bigP) + "_" + str(cycle_num * bigP) + "_days.csv", 0, "0", "")
   spatial.outputCSV("runner_y_P_" + str(bigP) + "_" + str(cycle_num * bigP) + "_days.csv", 1, "0", "")

sys.stdout.write("CPU time taken = " + str(timeit.default_timer() - start) + "s\n")

######################################################
# Simulate again, but with just diffusion and advection
sys.stdout.write("Doing only diffusion and advection\n")
start = timeit.default_timer()
timestep_size = 0.5
the_time = 0.0
end_time = 200.0
for bigP in [1]: # this is irrelevant - i just keep it here for ease of comparison with the previous block
   all_pops.setPopulationAndParametersFromXY(-2000.0, 700.0, [10000, 0, 0, 0]) # introduce 10000 type X at (-2000,700) (the population of Y and carrying capacities get overwritten by the following lines)
   num_cycles = int(end_time / bigP + 0.5)
   all_pops.setOverActiveGrid(1, initial_y)
   for cycle_num in range(num_cycles + 1):
      sys.stdout.write("  Cycle " + str(cycle_num) + " of " + str(num_cycles + 1) + "\n")
      #all_pops.setOverActiveGrid(2, cc_high)
      #all_pops.setOverActiveGrid(3, cc_high)
      for sub_cycle in range(bigP):
         sys.stdout.write("    K_high subcycle " + str(sub_cycle) + "\n")
         #### spatial.evolveCells(timestep_size)
         spatial.diffuse(timestep_size, 0.1)
         spatial.advect(array.array('f', [0.01]), wind)
         the_time += timestep_size
      #all_pops.setOverActiveGrid(2, cc_low)
      #all_pops.setOverActiveGrid(3, cc_low)
      for sub_cycle in range(bigP):
         sys.stdout.write("    K_low subcycle " + str(sub_cycle) + "\n")
         #### spatial.evolveCells(timestep_size)
         spatial.diffuse(timestep_size, 0.1)
         spatial.advect(array.array('f', [0.01]), wind)
         the_time += timestep_size

   sys.stdout.write("Outputting population result...\n")
   spatial.outputCSV("runner_x_no_cell_evolve_" + str(cycle_num * bigP) + "_days.csv", 0, "0", "")
   spatial.outputCSV("runner_y_no_cell_evolve_" + str(cycle_num * bigP) + "_days.csv", 1, "0", "")

sys.stdout.write("CPU time taken = " + str(timeit.default_timer() - start) + "s\n")


######################################################
# Simulate again, but with "x" always having cc_high and "y" always having cc_low
sys.stdout.write("Doing ODE evolution, diffusion and advection\n")
start = timeit.default_timer()
timestep_size = 0.5
the_time = 0.0
end_time = 200.0
for bigP in [1]: # this is irrelevant - i just keep it here for ease of comparison with the previous block
   all_pops.setPopulationAndParametersFromXY(-2000.0, 700.0, [10000, 0, 0, 0]) # introduce 10000 type X at (-2000,700) (the population of Y and carrying capacities get overwritten by the following lines)
   num_cycles = int(end_time / bigP + 0.5)
   all_pops.setOverActiveGrid(1, initial_y)
   for cycle_num in range(num_cycles + 1):
      sys.stdout.write("  Cycle " + str(cycle_num) + " of " + str(num_cycles + 1) + "\n")
      all_pops.setOverActiveGrid(2, cc_high)
      all_pops.setOverActiveGrid(3, cc_low)
      for sub_cycle in range(bigP):
         sys.stdout.write("    K_high subcycle " + str(sub_cycle) + "\n")
         spatial.evolveCells(timestep_size)
         spatial.diffuse(timestep_size, 0.1)
         spatial.advect(array.array('f', [0.01]), wind)
         the_time += timestep_size
      all_pops.setOverActiveGrid(2, cc_high)
      all_pops.setOverActiveGrid(3, cc_low)
      for sub_cycle in range(bigP):
         sys.stdout.write("    K_low subcycle " + str(sub_cycle) + "\n")
         spatial.evolveCells(timestep_size)
         spatial.diffuse(timestep_size, 0.1)
         spatial.advect(array.array('f', [0.01]), wind)
         the_time += timestep_size

   sys.stdout.write("Outputting population result...\n")
   spatial.outputCSV("runner_x_highcc_" + str(cycle_num * bigP) + "_days.csv", 0, "0", "")
   spatial.outputCSV("runner_y_highcc_" + str(cycle_num * bigP) + "_days.csv", 1, "0", "")

sys.stdout.write("CPU time taken = " + str(timeit.default_timer() - start) + "s\n")

sys.exit(0)
