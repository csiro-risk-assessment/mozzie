######################################################
# Calculates qm values in BevertonHolt dynamics for the given year
# It is assumed the carrying capacity is provided as CSV files
# - 5km x 5km grid defined that covers africa entirely (nx=1517, ny=1667)
# - active cells defined simply by 'ocean=inactive' 'land=active', resulting in 1.4 million cells
######################################################
import sys
try:
    yr = sys.argv[1]
except:
    raise IndexError('Need to specify the year as an input argument to this program')
sys.stdout.write("Calculating qm values for year " + yr + "\n")

######################################################
# Lifecycle parameters
num_species = 3
delay_days = 0 # In principal can have > 0 delay days, but the code below is slightly slower and more complicated, since it is necessary to put the equilibrium populations into the correct slots
death_rate_ww = [0.1, 0.1, 0.1] # death rate of wild-types of each species, measured in 1/day.  Below we assume the other genotypes have the same death rate: if a bad assumption then just modify death_rate variable below.  Must be changed if num_species changes from 3
competition = [[1, 0.4, 0.4], [0.4, 1, 0.4], [0.4, 0.4, 1]] # competition between subspecies.  Must be changed if num_species changes from 3
emergence_rate = [9.0, 9.0, 9.0] # emergence rate for each species.  Must be changed if num_species changes from 3
activity = [[1, 0, 0], [0, 1, 0], [0, 0, 1]] # activity[female_of_species1][male_of_species2] is activity level in the proportionate mixing.  Must be changed if num_species changes from 3
hM = 0.0 # reduction R = (1 - s * h) * (1 - s * h) for genotype = wc, for instance: see the reduct function below.  Using that function, it is assumed that hM = hF
sM = 0.0 # reduction R = (1 - s * h) * (1 - s * h) for genotype = wc, for instance: see the reduct function below.  Using that function, it is assumed that sM = sF
hybridisation = [[[1, 0, 0], [0, 1, 0], [0, 0, 1]], [[1, 0, 0], [0, 1, 0], [0, 0, 1]], [[1, 0, 0], [0, 1, 0], [0, 0, 1]]] # hybridisation[mM][mF][m] = prob that offspring of species m results from male of species mM and female of species mF.  The current value means offspring is always same as mF.  Must be changed if num_species changes from 3
offspring_modifier = [[[1, 1, 1], [1, 1, 1], [1, 1, 1]], [[1, 1, 1], [1, 1, 1], [1, 1, 1]]] # offspring_modifier[s][mM][mF] = suppression (if <1, or increased vigour if >1) of sppring of sex s that arises from male of species mM and female of species mF
sex_ratio = 0.7 # probability that offspring of wc or cc fathers are male
female_bias = 0.6 # probability that offspring of (mother wc or cc + father ww) is female
m_w = 0
m_c = 0
small_value = 1E-5 # when populations, carrying capacity, etc get smaller than this, set them to zero
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
g1 = Grid(-4614.0, -3967.0, 5.0, 1517, 1667, False)
sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")

sys.stdout.write("Setting active and inactive...")
start = timeit.default_timer()
g1.setActiveAndInactive(os.path.join("..", "example1", "active.csv"))
sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")
sys.stdout.write("There are " + str(g1.getNumActiveCells()) + " active cells\n")


######################################################
# define a zeroed populations and parameters array
sys.stdout.write("Defining the populations and parameters array...")
start = timeit.default_timer()
cell = CellDynamicsMosquitoBH26Delay(num_species = num_species, delay = delay_days, current_index = 0, death_rate = death_rate, competition = competition, emergence_rate = emergence_rate, activity = activity, reduction = reduction, hybridisation = hybridisation, offspring_modifier = offspring_modifier, sex_ratio = sex_ratio, female_bias = female_bias, m_w = m_w, m_c = m_c)
all_pops = PopulationsAndParameters(g1, cell)
sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")

for day in range(1, 366):
    sys.stdout.write("year = " + yr + ", day = " + str(day) + "\n")
    ######################################################
    # Read 3 files corresponding to carrying capacity for each of the 3 species
    # (Currently, the three files are actually identical files, but they should eventually be indexed with yr and day)
    #sys.stdout.write("Reading the 3 carrying capacity files and restricting to active cells...")
    start = timeit.default_timer()
    cc = [0] * num_species
    cc_parser = [0] * num_species
    for m in range(num_species):
        cc_parser[m] = SpatialDependence(-4614.0, -3967.0, 5.0, 1517, 1667)
        # NOTE: following line needs to read the carrying-capacity files for the given year and day.  At the moment it just reads the same carrying-capacity file over and over again!
        cc_parser[m].parse(os.path.join("..", "example1", "carrying.csv"), "generic_float", [])
        cc_parser[m].restrictToActive(g1.getGlobalIndex())
        cc[m] = cc_parser[m].getData0()
    #sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")

    ######################################################
    # introduce carrying capacities into the populations
    #sys.stdout.write("Populating the equilibrium populations and parameters array with wild-types...")
    start = timeit.default_timer()
    pop_and_param_array = all_pops.getQuantities()
    num_pops_and_params = cell.getNumberOfPopulations() + cell.getNumberOfParameters()
    for i in range(g1.getNumActiveCells()):
        for m in range(num_species):
            for g in range(1): # only ww genotype.  Usually this would be range(6)
                for s in range(2):
                    ind = m + g * num_species + s * num_species * 6 + i * num_pops_and_params
                    pop_and_param_array[ind] = max(1.0, 0.5 * cc[m][i]) # 0.5 is male-female split
    #sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")

    ######################################################
    # put these populations into the spatial structure
    #sys.stdout.write("Defining the spatial structure, ready for calculating qm...")
    start = timeit.default_timer()
    spatial = SpatialDynamics(g1, all_pops)
    #sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")

    ######################################################
    #sys.stdout.write("Calculating qm...")
    start = timeit.default_timer()
    spatial.calcQm()
    #sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")

    ######################################################
    #sys.stdout.write("Outputting result...")
    start = timeit.default_timer()
    for m in range(num_species):
        # NOTE: output filenames should be corrected
        spatial.outputCSV("qm_" + yr + "_" + str(day) + ".csv", cell.getNumberOfPopulations(), "0", "")
    #sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")

    ######################################################
    #sys.stdout.write("Converting to binary...")
    # NOTE: comment the following if you don't want binary files
    start = timeit.default_timer()
    for m in range(num_species):
        # NOTE: the filenames need to be corrected in the following (and check the 1517 1667 is correct)
        os.system(findbin + "/../code/auxillary/ab_convert ascii2binary 1517 1667 generic_float qm_" + yr + "_" + str(day) + ".csv qm_" + yr + "_" + str(day) + ".bin")
    #sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")
    

sys.exit(0)
