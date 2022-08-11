######################################################
# Example of calculating qm values in BevertonHolt dynamics
# It is assumed the carrying capacity is provided as CSV files
# - 5km x 5km grid defined that covers africa entirely (nx=1517, ny=1667)
# - active cells defined simply by 'ocean=inactive' 'land=active', resulting in 1.4 million cells
######################################################

######################################################
# Lifecycle parameters
species_list = ["Aa", "Ac", "Ag"] # if you change to !=3 subspecies, a few other things below will need to be changed, as noted below.  If you change these names, then search throughout this file for "species_list" to see what else you might need to change
num_species = 3
delay_days = 10 # number of days in the delay DE
death_rate_ww = [0.1, 0.1, 0.1] # death rate of wild-types of each species, measured in 1/day.  Below we assume the other genotypes have the same death rate: if a bad assumption then just modify death_rate variable below.  Must be changed if num_species changes from 3
competition = [[1, 0.717, 0.562], [0.717, 1, 0.317], [0.562, 0.317, 1]] # competition between subspecies.  Must be changed if num_species changes from 3
emergence_rate = [9.0, 9.0, 9.0] # emergence rate for each species.  Must be changed if num_species changes from 3
activity = [[0.9982032342, 0.0008983829, 0.0008983829], 
[0.0008974781, 0.9971978740, 0.0019046479], 
[0.0008974781, 0.0019046479, 0.9971978740]] # activity[female_of_species1][male_of_species2] is activity level in the proportionate mixing.  Must be changed if num_species changes from 3
hM = 0. #0.9 # reduction R = (1 - s * h) * (1 - s * h) for genotype = wc, for instance: see the recut function below.  Using that function, it is assumed that hM = hF
sM = 0. # 0.9 # reduction R = (1 - s * h) * (1 - s * h) for genotype = wc, for instance: see the recut function below.  Using that function, it is assumed that sM = sF
hybridisation = [[[1, 0, 0], [0, 1, 0], [0, 0, 1]], # Aa
[[1, 0, 0], [0, 1, 0], [0, 0, 1]],  # Ac
[[1, 0, 0], [0, 1, 0], [0, 0, 1]]] # Ag
# hybridisation[mM][mF][m] = prob that offspring of species m results from male of species mM and female of species mF.  The current value means offspring is always same as mF.  Must be changed if num_species changes from 3
offspring_modifier = [[[1, 1, 1], [1, 1, 1], [1, 1, 1]], [[1, 1, 1], [1, 1, 1], [1, 1, 1]]] # offspring_modifier[s][mM][mF] = suppression (if <1, or increased vigour if >1) of sppring of sex s that arises from male of species mM and female of species mF
sex_ratio = 0.950 # probability that offspring of wc or cc fathers are male
female_bias = 1. - 0.449 # probability that offspring of (mother wc or cc + father ww) is female
m_w = 0 # This is necessary for this burnin simulation: if m_w != 0 then resistant types are produced
m_c = 1E-6
small_value = 1E-5 # when populations, carrying capacity, etc get smaller than this, set them to zero
monthdays = [31,28,31,30,31,30,31,31,30,31,30,31] # TODO: adjust for LEAP years

death_rate = [[death_rate_ww] * 6] * 2
def reduct(genotype):
    if genotype == 0 or genotype == 2 or genotype == 5:
        return 1.0
    elif genotype == 1 or genotype == 4:
        return 1.0 - sM * hM
    elif genotype == 3:
        return 1.0 - sM
    raise ValueError("invalid genotype")
reduction = [[reduct(gM) * reduct(gF) for gM in range(6)] for gF in range(6)]

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
from cellDynamics import CellDynamicsMosquitoBH26Delay
from spatialDynamics import SpatialDynamics
from spatialDependence import SpatialDependence
from populationsAndParameters import PopulationsAndParameters


######################################################
# setup the grid and the active/inactive information
sys.stdout.write("Initialising grid (without building adjacency list)...")
start = timeit.default_timer()
g1 = Grid(-4614384.0, -3967418.0, 5000.0, 1517, 1667, False)
sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")

sys.stdout.write("Setting active and inactive...")
start = timeit.default_timer()
g1.setActiveAndInactive(os.path.join("..", "../simulation_xmas2021", "active2.csv"))
sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")
sys.stdout.write("There are " + str(g1.getNumActiveCells()) + " active cells\n")


######################################################
# define a zeroed populations and parameters array
# Note the use of the above parameters to define the lifecycle
sys.stdout.write("Defining the populations and parameters array...")
start = timeit.default_timer()
cell = CellDynamicsMosquitoBH26Delay(num_species = num_species, delay = delay_days, current_index = 0, death_rate = death_rate, competition = competition, 
emergence_rate = emergence_rate, activity = activity, reduction = reduction, hybridisation = hybridisation, offspring_modifier = offspring_modifier, 
sex_ratio = sex_ratio, female_bias = female_bias, m_w = m_w, m_c = m_c)
all_pops = PopulationsAndParameters(g1, cell)
sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")
sys.stdout.flush()

# NICK: start your loop over days here: modify the input filename (carrying.csv) and output filename (qm*.csv)
######################################################
# Read 3 files corresponding to carrying capacity for each of the 3 species
# (In this simple example, the three files are actually identical files)
sys.stdout.write("Reading the 3 carrying capacity files and restricting to active cells...")
start = timeit.default_timer()
cc = [0] * num_species
cc_parser = [0] * num_species
for year in range(2005, 2016):
   for monthno in range(12):
      month = str(monthno + 1).zfill(2)
      for dayno in range(monthdays[monthno]):
         day = str(dayno + 1).zfill(2)
         
         for m in range(num_species):
            species = species_list[m]
            fn = str(year) + month + day + "." + species + ".csv"
            cc_parser[m] = SpatialDependence(-4614384.0, -3967418.0, 5000.0, 1517, 1667)
            cc_parser[m].parse(os.path.join("..", "../K", fn), "generic_float", [])
            cc_parser[m].restrictToActive(g1.getGlobalIndex())
            cc[m] = cc_parser[m].getData0()
         sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")

         ######################################################
         # introduce carrying capacities into the populations
         sys.stdout.write("Populating the equilibrium populations and parameters array with wild-types...")
         start = timeit.default_timer()
         pop_and_param_array = all_pops.getQuantities()
         num_pops_and_params = cell.getNumberOfPopulations() + cell.getNumberOfParameters()
         for i in range(g1.getNumActiveCells()):
            for m in range(num_species):
               for g in range(1): # only ww genotype.  Usually this would be range(6)
                  for s in range(2):
                     ind = m + g * num_species + s * num_species * 6 + i * num_pops_and_params
                     pop_and_param_array[ind] = max(1.0, 0.5 * cc[m][i]) # 0.5 is male-female split
         sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")

         ######################################################
         # put these populations into the spatial structure
         sys.stdout.write("Defining the spatial structure, ready for calculating qm...")
         start = timeit.default_timer()
         spatial = SpatialDynamics(g1, all_pops)
         sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")

         ######################################################
         sys.stdout.write("Calculating qm...")
         start = timeit.default_timer()
         spatial.calcQm()
         sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")

         ######################################################
         sys.stdout.write("Outputting result...")
         start = timeit.default_timer()
         for m in range(num_species):
            species = species_list[m]
            fn = str(year) + month + day + "." + species + ".csv"
            spatial.outputCSV(os.path.join("..", "../q_m", fn), cell.getNumberOfPopulations(), "0", "")
         sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")
         sys.stdout.flush()

sys.exit(0)
