######################################################
# Import libraries
import os
import sys
import array
import timeit
import math
import numpy as np

######################################################
# Paths and filenames
#
# Append the path of the mozzie code directory to the system path
# Here, we assume you are running this script from mozzie/example_full
findbin = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(os.path.join(findbin, "..", "code"))
#
# Set the working directory,
# which contains files for defining active cells, etc, as well as the output
# Here, we assume your working directory is mozzie/example_full
working_dir = os.path.dirname(os.path.realpath(sys.argv[0]))
os.chdir(os.path.expanduser(working_dir))
# Directory, relative to working_directory, where the daily carrying-capacity files are
# These data will have been created by genCC.py
dir_with_cc = "cc"
# Directory, relative to working_directory, where the wind CSV files are
dir_with_wind = "wind"
# The file (in the working_directory) that contains the active cells
file_defining_active_cells = "active.csv"
# The output directory
output_dir = "output"

from wind import Wind
from grid import Grid
from cellDynamics import CellDynamicsMosquito26
from spatialDynamics import SpatialDynamics
from spatialDependence import SpatialDependence
from populationsAndParameters import PopulationsAndParameters

######################################################
# Setup the grid and the active/inactive information
#
sys.stdout.write("Initialising grid (without building adjacency list)...")
start = timeit.default_timer()
g1 = Grid(-1799134, -1299134, 5000.0, 100, 100, False)
sys.stdout.write("  Time taken = " + str(timeit.default_timer() - start) + "s\n")

sys.stdout.write("Setting active and inactive...")
start = timeit.default_timer()
g1.setActiveAndInactive(os.path.join(working_dir, file_defining_active_cells))
sys.stdout.write("  Time taken = " + str(timeit.default_timer() - start) + "s\n")
sys.stdout.write("There are " + str(g1.getNumActiveCells()) + " active cells\n")


###################################################
# Model parameters
#
# Introduction location roughly in middle of grid
locx = -1799134 + 2500 + 49*5000
locy = -1299134 + 2500 + 49*5000
# Following are the parameters and input-files of this model
species_list = ["Ac", "Ag"] # if you change to a different number of species, a few other things below will need to be changed, as noted below.  If you change these names, then search throughout this file for "species_list" to see what else you might need to change
sex_list = ["male", "female"] # must be 2 sexes
genotype_list = ["ww", "wc", "wr", "cc", "cr", "rr"] # must be 6 genotypes because of the use of CellDynamicsMosquito26
num_species = len(species_list) 
# death rate in 1/day for each sex, genotype and species.  Must be changed if num_species changes from 2.
death_rate = [[[0.1 for m in range(2)] for g in range(6)] for s in range(2)]
competition = [[1.00000000,0.01], [0.01,1.00000000]] # competition between subspecies.  Must be changed if num_species changes from 2
emergence_rate = [9.0, 9.0] # emergence rate in 1/day.  Must be changed if num_species changes from 2
# activity[female_of_species1][male_of_species2] is activity level in the proportionate mixing.  Must be changed if num_species changes from 2
activity = [[0.998093641, 0.001906359], 
            [0.001906359, 0.998093641]]
# hybridisation[mM][mF][m] = probability that offspring of species m results from male of species mM and female of species mF.  The current value means offspring is always same as mF.  Must be changed if num_species changes from 2
hybridisation = [[[[1, 0], [0, 1]],  # Ac
                  [[1, 0], [0, 1]]]] # Ag
# offspring_modifier[s][mM][mF] = suppression (if <1, or increased vigour if >1) offspring of sex s that arises from male of species mM and female of species mF
offspring_modifier = [[[1, 1], [1, 1]], [[1, 1], [1, 1]]]
# diffusion coefficient in km^2/day.  This is massively turbocharged compared with something that is realistic, to best illustrate mozzie's capability.  Or, equivalently, the grid could be measured in a smaller unit than km
diffusion_coefficient = 100 * 190.0**2 / (math.pi * 7)
# wind, and advection fraction
windhrs = 2
if windhrs <= 0:
    advection_fraction = 0
elif windhrs == 2:
    advection_fraction = 2.6E-3
elif windhrs == 9:
    advection_fraction = 6.4E-4
else:
    raise ValueError("windhrs must be <=0, 2 or 9")

# Example simulation goes from 2022-3
years_simulated = list(range(2022, 2024))
year_begin = 2022
monthdays = [31,28,31,30,31,30,31,31,30,31,30,31] # No leap years

# add output folder
try:
    os.mkdir(os.path.join(working_dir, output_dir))
except:
    sys.stdout.write("Output folder already exists\n")
# If num_species != 2, you shouldn't have to change anything below here

# 26May2025: Andy to here: need to use all the above parameters!
######################################################
# Define the ODE
sys.stdout.write("Defining the ODE\n")
start = timeit.default_timer()
cell = CellDynamicsMosquito26()
cell.setTimeIntegrationMethod("runge_kutta4")
cell.setMinCarryingCapacity(1.0)
cell.setZeroCutoff(1.0)
        # h_e = h_n = 0.5
        # s_e = 0.1
        # s_n = 0.05
cell.setFitnessComponents26(0.0, 0.5, 0.0, 0.05) # NO EFFECTOR (h_e, h_n, s_e, s_n)
cell.setHybridisationRate(0, 1, 1, 1.0) 
cell.setHybridisationRate(0, 1, 0, 0.0) 
cell.setHybridisationRate(1, 0, 0, 1.0)
cell.setHybridisationRate(1, 0, 1, 0.0) 
cell.setAlphaComponent(0, 1, 0.1)
cell.setAlphaComponent(1, 0, 0.1)
sys.stdout.write("  Time taken = " + str(timeit.default_timer() - start) + "s\n")

######################################################
# define a zeroed populations and parameters array
start = timeit.default_timer()
sys.stdout.write("Defining the populations and parameters array...")
all_pops = PopulationsAndParameters(g1, cell)
sys.stdout.write("  Time taken = " + str(timeit.default_timer() - start) + "s\n")

######################################################
# initialise populations
# read a file corresponding to the cc for the first day, to set initial populations
start = timeit.default_timer()
sys.stdout.write("Reading a representative cc to set the initial populations...")
ini_cc_parsers = [0 for i in range(num_species)]
ini_ccs = [0 for i in range(num_species)]
first_day = str(years_simulated[0]) + "0101"

# ######################################################
# # load spatial variables

######################################################
# load first spatiotemporal variables

year = year_begin 

temp = SpatialDependence(-1799134, -1299134, 5000.0, 100, 100) 
temp.parse(os.path.join(working_dir, dir_with_cc, "cc.Ac." + str(year_begin) + "." + str(1) + ".csv"), "generic_float", [])
temp.restrictToActive(g1.getGlobalIndex())
ini_cc = temp.getData0()
num_actives = len(ini_cc) 
ini_ccs[0] = [ini_cc[i] + 1.e-3 for i in range(num_actives)]

temp = SpatialDependence(-1799134, -1299134, 5000.0, 100, 100) 
temp.parse(os.path.join(working_dir, dir_with_cc, "cc.Ag." + str(year_begin) + "." + str(1) + ".csv"), "generic_float", [])
temp.restrictToActive(g1.getGlobalIndex())
ini_cc = temp.getData0()
ini_ccs[1] = [ini_cc[i] + 1.e-3 for i in range(num_actives)]
sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")
sys.stdout.write("Initialising the populations, with male and female wild-types at the carrying capacity...\n")
start = timeit.default_timer()
num_actives = len(ini_ccs[0])

ini_pops = [array.array('f', [100.0 * 0.5 * ini_ccs[m][i] for i in range(num_actives)]) for m in range(num_species)]

for species in range(num_species):
    for genotype in [0]: # only wild-types
        for sex in range(2):
            for age in range(2):
                index = species + genotype * num_species + sex * num_species * 6 + age * num_species * 12
                all_pops.setOverActiveGrid(index, ini_pops[species])

sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")

######################################################
# put these populations into the spatial structure
sys.stdout.write("Defining the spatial structure, ready for diffusion, advection and ODE evolution...")
start = timeit.default_timer()
spatial = SpatialDynamics(g1, all_pops)
sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")
sys.stdout.flush()


######################################################
# Simulate
sys.stdout.write("Doing ODE evolution, diffusion and advection\n")
sys.stdout.flush()
start = timeit.default_timer()
num_diffusing = g1.getNumActiveCells() * cell.getNumberOfDiffusingPopulations()

timestep_size = 1.0

for year in years_simulated:
    doyno = 0
    for monthno in range(12):
        month = str(monthno + 1).zfill(2)
        
        ######################################################
        # define wind for the month
        sys.stdout.write("Reading and processing wind for month " + month + "...")
        sys.stdout.flush() 
        start = timeit.default_timer()
        if windhrs > 0:
            pdf = [[float(windhrs) / 24.0, 1.0]]
            wind = {}
            for dayno in range(monthdays[monthno]):
                if year >= year_begin:
                    day = str(dayno + 1).zfill(2)
                    fn = str(year) + "." + month + "." + day + ".windvec_m_per_day.csv"
                    processed_fn = "processed_windhrs_" + str(windhrs) + "_" + str(year) + "." + month + "." + day + ".bin"
                    wind[(month, day)] = Wind(os.path.join(working_dir, dir_with_wind, fn), os.path.join(working_dir, dir_with_wind, processed_fn), pdf, g1)
                    wind[(month, day)].setBinaryFileFormat(0) # just using csvs, haven't processed binary versions for example (but could do for speedup, then use parseProcessedFile()
                    wind[(month,day)].parseRawFile()
                    sys.stdout.write("    Day " + day + "\n")
                    sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")
                    sys.stdout.flush() 
        sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")
        sys.stdout.flush() 

        ######################################################
        # define cc for the month
        sys.stdout.write("Reading and processing cc for month " + month + "...")
        sys.stdout.flush() 
        start = timeit.default_timer()
        cc = {}
        for dayno in range(monthdays[monthno]):
                doyno = doyno + 1
                day = str(dayno + 1).zfill(2)
                for speciesno in range(num_species):
                    species = species_list[speciesno]
                    if speciesno == 0:
                        temp = SpatialDependence(-1799134, -1299134, 5000.0, 100, 100) 
                        temp.parse(os.path.join(working_dir, dir_with_cc, "cc.Ac." + str(year) + "." + str(doyno) + ".csv"), "generic_float", [])
                        temp.restrictToActive(g1.getGlobalIndex())
                        ini_cc = temp.getData0()
                        cc[(month, day, species)] = [100.0 * ini_cc[i] + 1.e-3 for i in range(num_actives)]
                    elif speciesno == 1:
                        temp = SpatialDependence(-1799134, -1299134, 5000.0, 100, 100) 
                        temp.parse(os.path.join(working_dir, dir_with_cc, "cc.Ag." + str(year) + "." + str(doyno) + ".csv"), "generic_float", [])
                        temp.restrictToActive(g1.getGlobalIndex())
                        ini_cc = temp.getData0()
                        cc[(month, day, species)] = [100.0 * ini_cc[i] + 1.e-3 for i in range(num_actives)]
                        
                    cc[(month, day, species)] = array.array('f', cc[(month, day, species)])

        sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")
        sys.stdout.flush() 

        ######################################################
        # evolve over the month
        for dayno in range(monthdays[monthno]):
                day = str(dayno + 1).zfill(2)
                sys.stdout.write("    Day " + day + "\n")

                all_pops.setOverActiveGrid(48, cc[(month, day, "Ac")]) # set carrying capacity for males AND females for this month
                all_pops.setOverActiveGrid(49, cc[(month, day, "Ag")])

                if year == 2022 and (monthno == 5) and (dayno == 0): # introduce on 1 June 2022
                    # introduce 10,000 genetically modified mosquitoes
                    pap = all_pops.getPopulationAndParametersFromXY(locx, locy)
                    # species 0 = Ac, genotype 1 = wc, sex 0 = male, age 1 = adult
                    species = 0
                    genotype = 1
                    sex = 0
                    age = 1
                    index = species + genotype * num_species + sex * num_species * 6 + age * num_species * 12
                    pap[index] = pap[index] + 100.0 * 10000.0
                    all_pops.setPopulationAndParametersFromXY(locx,locy, pap.tolist())

                spatial.evolveCells(timestep_size)
                sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")
                spatial.diffuse(timestep_size, diffusion_coefficient)
                sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")
                    
                spatial.advect(array.array('f', [advection_fraction]), wind[(month, day)])
                    
                sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")
                
                sys.stdout.write("Outputting population result...\n")
                
                age = 1 # adults only
                for genotype in range(6):
                    for species in range(num_species):
                        for sex in range(2):
                            ind = species + genotype * num_species + sex * num_species * 6 + age * num_species * 12
                            spatial.outputCSV(os.path.join(working_dir, output_dir, species_list[species] + "_" + genotype_list[genotype] + "_" + sex_list[sex] + "_" + str(year) + "_" + month + "_" + day + ".csv"), ind, "0", "")
                
                sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")
            
                sys.stdout.flush() # every day
        
    
        sys.stdout.write("CPU time taken = " + str(timeit.default_timer() - start) + "s\n")

sys.exit(0)
