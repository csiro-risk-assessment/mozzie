######################################################
# This demo uses CellDynamics26 with 2 species
# The simulation runs over the years 2022 and 2023
# with daily changes in spatially-varying carrying capacity
# and wind vectors.
# Genetically modified mosquitoes are released near the
# center of the domain during June 2022
# This demo models the spread of the genetic modification

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
g1 = Grid(-25000, -25000, 500.0, 100, 100, False)
sys.stdout.write(f"  Time taken = {(timeit.default_timer() - start):.2g}s\n")

sys.stdout.write("Setting active and inactive...")
start = timeit.default_timer()
g1.setActiveAndInactive(os.path.join(working_dir, file_defining_active_cells))
sys.stdout.write(f"  Time taken = {(timeit.default_timer() - start):.2g}s\n")
sys.stdout.write("There are " + str(g1.getNumActiveCells()) + " active cells\n")


###################################################
# Model parameters
#
species_list = ["Ac", "Ag"] # if you change to a different number of species, a few other things below will need to be changed, as noted below.  If you change these names, then search throughout this file for "species_list" to see what else you might need to change
sex_list = ["male", "female"] # must be 2 sexes
genotype_list = ["ww", "wc", "wr", "cc", "cr", "rr"] # must be 6 genotypes because of the use of CellDynamicsMosquito26
num_species = len(species_list) # useful
number_species_name = [(i, species_list[i]) for i in range(num_species)] # useful
# multiply the carrying capacities found in the *cc* files by this number
cc_scale = 0.5
# competition between subspecies.  Must be changed if num_species changes from 2
competition = [[1.0, 0.1], [0.1, 1.0]]
# hybridisation[mM][mF][m] = probability that offspring of species m results from male of species mM and female of species mF.  The current value means offspring is always same as mF.  Must be changed if num_species changes from 2
hybridisation = [[[1.0, 0.0],   # father Ac, mother Ac
                  [0.1, 0.9]],  # father Ac, mother Ag
                 [[0.9, 0.1],   # father Ag, mother Ac
                  [0.0, 1.0]]]  # father Ag, mother Ag
# diffusion coefficient in m^2/day (roughly equivalent to 30 m / day average mozzie movement)
diffusion_coefficient = 900 
# number of hours mosquitoes advect by wind, and advection fraction
windhrs = 2
advection_fraction = 0.005 #2.6E-3
# Introduction location roughly in middle of grid
locx = 0
locy = 0
# If num_species != 2, you shouldn't have to change anything below here

# Example simulation runs through 2022 and 2023
years_simulated = list(range(2022, 2024))
year_begin = years_simulated[0]
monthdays = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31] # No leap years

# add output folder
try:
    os.mkdir(os.path.join(working_dir, output_dir))
except:
    sys.stdout.write("Output folder already exists\n")

######################################################
# Define the ODE
sys.stdout.write("Defining the ODE...")
start = timeit.default_timer()
cell = CellDynamicsMosquito26()
# Change the default parameters to ones appropriate to this demo
# NO EFFECTOR (h_e, h_n, s_e, s_n)
cell.setFitnessComponents26(0.0, 0.5, 0.0, 0.05)
# offspring species is always equal to mother's species
for father in range(num_species):
    for mother in range(num_species):
        for offspring in range(num_species):
            cell.setHybridisationRate(father, mother, offspring, hybridisation[father][mother][offspring])
# competition
for sp0 in range(num_species):
    for sp1 in range(num_species):
        cell.setAlphaComponent(sp0, sp1, competition[sp0][sp1])
cell.setMinCarryingCapacity(1E-2)
cell.setZeroCutoff(1E-2)
sys.stdout.write(f"  Time taken = {(timeit.default_timer() - start):.2g}s\n")

######################################################
# Define a zeroed populations and parameters array
start = timeit.default_timer()
sys.stdout.write("Defining a zeroed populations and parameters array...")
all_pops = PopulationsAndParameters(g1, cell)
sys.stdout.write(f"  Time taken = {(timeit.default_timer() - start):.2g}s\n")

######################################################
# Initialise populations
# read a file corresponding to the cc for the first day, to set initial populations
start = timeit.default_timer()
sys.stdout.write("Reading a representative carrying capacity to set the initial populations...")
temp = SpatialDependence(-25000, -25000, 500.0, 100, 100) 
ini_ccs = [0 for i in range(num_species)]
for species_num, species_name in number_species_name:
    fn = os.path.join(working_dir, dir_with_cc, "cc." + species_name + "." + str(year_begin) + ".1.csv")
    if not os.path.exists(fn):
        sys.stdout.flush()
        sys.stderr.write("\n\nError: File " + fn + " does not exist.  You probably need to run genCC.py\nExiting\n")
        sys.exit(1)
    temp.parse(fn, "generic_float", [])
    temp.restrictToActive(g1.getGlobalIndex())
    ini_cc = temp.getData0()
    num_actives = len(ini_cc) 
    ini_ccs[species_num] = [ini_cc[i] + 1.e-3 for i in range(num_actives)]
sys.stdout.write(f"  Time taken = {(timeit.default_timer() - start):.2g}s\n")
# Initialise the populations
sys.stdout.write("Initialising the populations, with male and female wild-types at the carrying capacity...")
start = timeit.default_timer()
num_actives = len(ini_ccs[0])
ini_pops = [array.array('f', [cc_scale * 0.5 * ini_ccs[m][i] for i in range(num_actives)]) for m in range(num_species)]
for species in range(num_species):
    for genotype in [0]: # only wild-types
        for sex in range(2):
            for age in range(2):
                index = species + genotype * num_species + sex * num_species * 6 + age * num_species * 12
                all_pops.setOverActiveGrid(index, ini_pops[species])
sys.stdout.write(f"  Time taken = {(timeit.default_timer() - start):.2g}s\n")

######################################################
# put these populations into the spatial structure
sys.stdout.write("Defining the spatial structure, ready for diffusion, advection and ODE evolution...")
start = timeit.default_timer()
spatial = SpatialDynamics(g1, all_pops)
sys.stdout.write(f"  Time taken = {(timeit.default_timer() - start):.2g}s\n")
sys.stdout.flush()


######################################################
# Simulate
sys.stdout.write("\n")
sys.stdout.write("############################################\n")
sys.stdout.write("Doing ODE evolution, diffusion and advection\n")
sys.stdout.flush()
timestep_size = 1.0

for year in years_simulated:
    day_in_year = 0
    sys.stdout.write(str(year) + "\n")
    for monthno in range(12):
        month = str(monthno + 1).zfill(2)
        sys.stdout.write("  Month " + month + "\n")
        
        ######################################################
        # define wind for the month
        sys.stdout.write("    Reading and processing wind ")
        start = timeit.default_timer()
        pdf = [[float(windhrs) / 24.0, 1.0]]
        wind = {}
        for dayno in range(monthdays[monthno]):
            day = str(dayno + 1).zfill(2)
            fn = os.path.join(working_dir, dir_with_wind, str(year) + "." + month + "." + day + ".windvec_m_per_day.csv")
            processed_fn = os.path.join(working_dir, dir_with_wind, "processed_windhrs_" + str(windhrs) + "_" + str(year) + "." + month + "." + day + ".bin")
            if not os.path.exists(fn):
                sys.stdout.flush()
                sys.stderr.write("\n\nError: File " + fn + " does not exist.  You probably need to run genWind.py\nExiting\n")
                sys.exit(1)
            wind[(month, day)] = Wind(fn, processed_fn, pdf, g1)
            wind[(month, day)].setBinaryFileFormat(0) # just using csvs, haven't processed binary versions for example (but could do for speedup, then use parseProcessedFile()
            wind[(month,day)].parseRawFile()
            sys.stdout.write(".")
            sys.stdout.flush() 
        sys.stdout.write(f"  Time taken = {(timeit.default_timer() - start):.2g}s\n")
        sys.stdout.flush()

        ######################################################
        # define carrying capacity for the month
        sys.stdout.write("    Reading and processing cc ")
        sys.stdout.flush() 
        start = timeit.default_timer()
        temp = SpatialDependence(-25000, -25000, 500.0, 100, 100) 
        cc = {}
        for dayno in range(monthdays[monthno]):
            day = str(dayno + 1).zfill(2)
            day_in_year += 1
            for species_num, species_name in number_species_name:
                fn = os.path.join(working_dir, dir_with_cc, "cc." + species_name + "." + str(year) + "." + str(day_in_year) + ".csv")
                if not os.path.exists(fn):
                    sys.stdout.flush()
                    sys.stderr.write("\n\nError: File " + fn + " does not exist.  You probably need to run genCC.py\nExiting\n")
                    sys.exit(1)
                temp.parse(fn, "generic_float", [])
                temp.restrictToActive(g1.getGlobalIndex())
                ini_cc = temp.getData0()
                cc[(month, day, species_name)] = [cc_scale * ini_cc[i] + 1.e-3 for i in range(num_actives)]
                cc[(month, day, species_name)] = array.array('f', cc[(month, day, species_name)])
            sys.stdout.write(".")
            sys.stdout.flush() 
        sys.stdout.write(f"  Time taken = {(timeit.default_timer() - start):.2g}s\n")
        sys.stdout.flush()

        ######################################################
        # evolve over the month
        sys.stdout.write("    Evolving ")
        sys.stdout.flush() 
        start = timeit.default_timer()
        for dayno in range(monthdays[monthno]):
            day = str(dayno + 1).zfill(2)

            # set carrying capacity for males and females for this month and day
            all_pops.setOverActiveGrid(48, cc[(month, day, "Ac")])
            all_pops.setOverActiveGrid(49, cc[(month, day, "Ag")])

            # introduce genetically modified mosquitoes on 1 June 2022
            if (year == 2022) and (monthno == 5) and (dayno == 0):
                sys.stdout.write("\n    Introducing genetically modified mosquitoes!\n")
                # introduce 10000 * cc_scale genetically modified mosquitoes
                pap = all_pops.getPopulationAndParametersFromXY(locx, locy)
                # species 0 = Ac, genotype 1 = wc, sex 0 = male, age 1 = adult
                species = 0
                genotype = 1
                sex = 0
                age = 1
                index = species + genotype * num_species + sex * num_species * 6 + age * num_species * 12
                pap[index] = pap[index] + cc_scale * 10000.0
                all_pops.setPopulationAndParametersFromXY(locx, locy, pap.tolist())
                sys.stdout.write("    Now continuing evolution ")
                sys.stdout.flush() 
                
            # lifecycle evolution
            sys.stdout.write(".")
            sys.stdout.flush() 
            spatial.evolveCells(timestep_size)
            # diffusion
            sys.stdout.write("\b@")
            sys.stdout.flush() 
            spatial.diffuse(timestep_size, diffusion_coefficient)
            # advection
            sys.stdout.write("\b>")
            sys.stdout.flush() 
            spatial.advect(array.array('f', [advection_fraction]), wind[(month, day)])

            # output
            sys.stdout.write("\bo")
            sys.stdout.flush() 
            age = 1 # adults only
            for genotype in range(6):
                for species in range(num_species):
                    for sex in range(2):
                        ind = species + genotype * num_species + sex * num_species * 6 + age * num_species * 12
                        spatial.outputCSV(os.path.join(working_dir, output_dir, species_list[species] + "_" + genotype_list[genotype] + "_" + sex_list[sex] + "_" + str(year) + "_" + month + "_" + day + ".csv"), ind, "0", "")
                
        sys.stdout.write(f"  Time taken = {(timeit.default_timer() - start):.2g}s\n")
        sys.stdout.flush()


sys.exit(0)
