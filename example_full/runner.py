######################################################
# Import libraries
import os
import sys
import array
import timeit
import math
import numpy as np

######################################################
# User input directories
# Set this to the path of your working directory
working_dir = "/scratch3/bee069/bee069/example_full" #"SET PATH HERE"
os.chdir(os.path.expanduser(working_dir))

# Set this to the path of the mozzie code directory
sys.path.append("/home/bee069/mozzie/mozzie/code")#"SET PATH HERE")

######################################################
# Default example directories
# Please enter the path of the directory where to read the .bin files of q
dir_with_cc = "cc" # corresponding to the "cc" of each species at each year at each month at each day.
# Please enter the path of the directory where to read the .bin files of processed wind fields
dir_with_wind = "wind/"
# Please enter the path of the directory where to read the .csv file of the active grid
file_defining_active_cells = "activeDemo.csv"
# Please enter the name of the output directory
output_dir = "output/"

from wind import Wind
from grid import Grid
from cellDynamics import CellDynamicsMosquito26
from spatialDynamics import SpatialDynamics
from spatialDependence import SpatialDependence
from populationsAndParameters import PopulationsAndParameters

######################################################
# setup the grid and the active/inactive information
sys.stdout.write("Initialising grid (without building adjacency list)...")
start = timeit.default_timer()
g1 = Grid(-1799134, -1299134, 5000.0, 100, 100, False)
sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")

sys.stdout.write("Setting active and inactive...")
start = timeit.default_timer()
g1.setActiveAndInactive(os.path.join(working_dir, file_defining_active_cells))
sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")
sys.stdout.write("There are " + str(g1.getNumActiveCells()) + " active cells\n")


###################################################
# Model parameters

# Introduction location in middle of grid somewhere
locx = -1799134 + 2500 + 49*5000
locy = -1299134 + 2500 + 49*5000


# Following are the parameters and input-files of this model
species_list = ["Ac", "Ag"] # if you change to !=3 subspecies, a few other things below will need to be changed, as noted below.  If you change these names, then search throughout this file for "species_list" to see what else you might need to change
sex_list = ["male", "female"] # must be 2 sexes
genotype_list = ["ww", "wc", "wr", "cc", "cr", "rr"] # must be 6 genotypes
num_species = len(species_list) 
death_rate_ww = [0.1, 0.1] # death rate of wild-types of each species, measured in 1/day.  Below we assume the other genotypes have the same death rate: if a bad assumption then just modify death_rate variable below.  Must be changed if num_species changes from 3
competition = [[1.00000000,0.01], [0.01,1.00000000]] # competition between subspecies.  Must be changed if num_species changes from 3
emergence_rate = [9.0, 9.0] 
activity = [[0.998093641, 0.001906359], 
[0.001906359, 0.998093641]]
# activity[female_of_species1][male_of_species2] is activity level in the proportionate mixing.  Must be changed if num_species changes from 3

hybridisation = [[[[1, 0], [0, 1]],  # Ac
[[1, 0], [0, 1]]]] # Ag
# hybridisation[mM][mF][m] = prob that offspring of species m results from male of species mM and female of species mF.  The current value means offspring is always same as mF.  Must be changed if num_species changes from 3
offspring_modifier = [[[1, 1], [1, 1]], [[1, 1], [1, 1]]]
# offspring_modifier[s][mM][mF] = suppression (if <1, or increased vigour if >1) offspring of sex s that arises from male of species mM and female of species mF

diffusion_coefficient = 100 * 190.0**2 / (math.pi * 7) # TURBOCHARGED 100 x 1641.57 m^2/day for male Ac wildtype

windhrs = 2;
if windhrs <= 0:
    advection_fraction = 0
elif windhrs == 2:
    advection_fraction = 2.6E-3
elif windhrs == 9:
    advection_fraction = 6.4E-4
else:
    raise ValueError("windhrs must be <=0, 2 or 9")

######################################################
# 50/50 sex advection ratio
advection_males = 0.5 

# Example simulation goes from 2022-3
years_simulated = list(range(2022, 2024))
year_begin = 2022
monthdays = [31,28,31,30,31,30,31,31,30,31,30,31] # No leap years

death_rate = [[[0.1 for m in range(2)] for g in range(6)] for s in range(2)]


# add output folder
try:
    os.mkdir(os.path.join(working_dir, output_dir))
except:
    sys.stdout.write("Output folder already exists\n")

# End parameter definitions
# If num_species != 3, you shouldn't have to change anything below here


######################################################
# Define the ODE
sys.stdout.write("Defining the populations and parameters array...")
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
self.setAlphaComponent(0, 1, 0.1)
self.setAlphaComponent(1, 0, 0.1)

######################################################
# define a zeroed populations and parameters array
start = timeit.default_timer()
sys.stdout.write("defining the populations and parameters array...")
all_pops = PopulationsAndParameters(g1, cell)
# get coords of introduction cell
sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")

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

ini_pops = [array.array('f', [0.5 * ini_ccs[m][i] for i in range(num_actives)]) for m in range(num_species)]

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
                        cc[(month, day, species)] = [ini_cc[i] + 1.e-3 for i in range(num_actives)]
                    elif speciesno == 1:
                        temp = SpatialDependence(-1799134, -1299134, 5000.0, 100, 100) 
                        temp.parse(os.path.join(working_dir, dir_with_cc, "cc.Ag." + str(year) + "." + str(doyno) + ".csv"), "generic_float", [])
                        temp.restrictToActive(g1.getGlobalIndex())
                        ini_cc = temp.getData0()
                        cc[(month, day, species)] = [ini_cc[i] + 1.e-3 for i in range(num_actives)]
                        
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
                    pap[index] = pap[index] + 10000.0
                    all_pops.setPopulationAndParametersFromXY(locx,locy, pap.tolist())

                spatial.evolveCells(timestep_size)
                sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")
                spatial.diffuse(timestep_size, diffusion_coefficient)
                sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")
                    
                if windhrs == 9:
                    spatial.advect(array.array('f', [26.65/41526.]), wind[(month, day)])
                elif windhrs == 2:
                    spatial.advect(array.array('f', [106.39/41526.]), wind[(month, day)])
                    
                sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")
                
                sys.stdout.write("Outputting population result...\n")
                
                age = 1 # adults only
                for genotype in range(6):
                    for species in range(num_species):
                        for sex in range(2):
                            ind = species + genotype * num_species + sex * num_species * 6 + age * num_species * 12
                            spatial.outputCSV(output_dir + species_list[species] + "_" + genotype_list[genotype] + "_" + sex_list[sex] + "_" + str(year) + "_" + month + "_" + day + ".csv", ind, "0", "")
                
                sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")
            
                sys.stdout.flush() # every day
        
    
        sys.stdout.write("CPU time taken = " + str(timeit.default_timer() - start) + "s\n")

sys.exit(0)
