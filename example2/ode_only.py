######################################################
# Example simulation of Beeton, Hosack, Wilkins, Forbes, Ickowicz and Hayes, Journal of Theoretical Biology 2019
# performed on a single grid-cell only.
# That is, just the ODEs are solved


######################################################
# Import libraries
import os
import sys
import timeit
import array

# so we can find our ../code no matter how we are called
findbin = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(findbin + "/../code")

from grid import Grid
from spatialDynamics import SpatialDynamics
from cellDynamics import CellDynamicsBeeton2_2
from populationsAndParameters import PopulationsAndParameters


######################################################
# setup the grid, which is just 1 active cell
sys.stdout.write("Initialising grid (without building adjacency list)...\n")
g1 = Grid(0.0, 0.0, 1.0, 1, 1, False)


######################################################
# define a zeroed populations and parameters array
sys.stdout.write("Defining the populations and parameters array...\n")
cell = CellDynamicsBeeton2_2()
all_pops = PopulationsAndParameters(g1, cell)

# introduce carrying capacity (which, for CellDynamicsLogistic1_1, is second slot in the populations and params array), which is max(1.0, cc)
sys.stdout.write("Populating the initial populations and parameters array with zeroes...\n")
start = timeit.default_timer()
pop_and_param_array = all_pops.getQuantities()



######################################################
# put these populations into the spatial structure
sys.stdout.write("Defining the spatial structure, ready for diffusion, advection and logistic growth...\n")
spatial = SpatialDynamics(g1, all_pops)


######################################################
sys.stdout.write("Doing cell dynamics\n")
start = timeit.default_timer()
timestep_size = 1.0
total_num_steps = 200
for bigP in [20, 1]:
   sys.stdout.write("P = " + str(bigP) + "\n")
   # set initial conditions
   pop_and_param_array[0] = 0.02 # x
   pop_and_param_array[1] = 0.1   # y
   f = open("ode_only_P_" + str(bigP) + ".csv", "w")
   f.write("#time,x,y,kx,ky\n")
   outstr = ""
   the_time = 0.0
   num_cycles = total_num_steps // bigP
   for cycle_num in range(num_cycles):
      sys.stdout.write("  Cycle " + str(cycle_num) + "\n")
      pop_and_param_array[2] = 1.0   # k_x
      pop_and_param_array[3] = 1.0   # k_y
      for sub_cycle in range(bigP):
         spatial.evolveCells(timestep_size)
         the_time += timestep_size
         outstr += str(the_time) + "," + ",".join([str(p) for p in pop_and_param_array]) + "\n"
      pop_and_param_array[2] = 0.2   # k_x
      pop_and_param_array[3] = 0.2   # k_y
      for sub_cycle in range(bigP):
         spatial.evolveCells(timestep_size)
         the_time += timestep_size
         outstr += str(the_time) + "," + ",".join([str(p) for p in pop_and_param_array]) + "\n"
   f.write(outstr)
   f.close()

sys.stdout.write("CPU time taken = " + str(timeit.default_timer() - start) + "s\n")

sys.exit(0)
