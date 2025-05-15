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
dir_with_qm = "qm" # corresponding to the "q_m" of each species at each year at each month at each day.
# Please enter the path of the directory where to read the .bin files of processed wind fields
dir_with_wind = "wind/"
# Please enter the path of the directory where to read the .csv file of the active grid
file_defining_active_cells = "activeDemo.csv"
# Please enter the name of the output directory
output_dir = "output/"

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
g1 = Grid(-1799134, -1299134, 5000.0, 100, 100, False)
sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")

sys.stdout.write("Setting active and inactive...")
start = timeit.default_timer()
g1.setActiveAndInactive(os.path.join(working_dir, file_defining_active_cells))
sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")
sys.stdout.write("There are " + str(g1.getNumActiveCells()) + " active cells\n")


###################################################
# Model parameters
m_c = 0.0001136718  # taken from fault tree distribution for construct failure
# Following are the parameters and input-files of this model
species_list = ["Aa", "Ac", "Ag"] # if you change to !=3 subspecies, a few other things below will need to be changed, as noted below.  If you change these names, then search throughout this file for "species_list" to see what else you might need to change
sex_list = ["male", "female"] # must be 2 sexes
genotype_list = ["ww", "wc", "wr", "cc", "cr", "rr"] # must be 6 genotypes
num_species = len(species_list) 
death_rate_ww = [0.1, 0.1, 0.1] # death rate of wild-types of each species, measured in 1/day.  Below we assume the other genotypes have the same death rate: if a bad assumption then just modify death_rate variable below.  Must be changed if num_species changes from 3
competition = [[1.00000000,0.22,0.37], [0.22,1.00000000,0.23], [0.37,0.23,1.00000000]] # competition between subspecies.  Must be changed if num_species changes from 3
emergence_rate = [9.0, 9.0, 9.0] # no longer need conversion as changed to Geoff method
activity = [[0.9982032342, 0.0008983829, 0.0008983829], 
[0.0008974781, 0.9971978740, 0.0019046479], 
[0.0008974781, 0.0019046479, 0.9971978740]]
# activity[female_of_species1][male_of_species2] is activity level in the proportionate mixing.  Must be changed if num_species changes from 3
hM = 1.0 # 0. #0.9 # reduction R = (1 - s * h) * (1 - s * h) for genotype = wc, for instance: see the recut function below.  Using that function, it is assumed that hM = hF
# 0. # 0.9 # reduction R = (1 - s * h) * (1 - s * h) for genotype = wc, for instance: see the recut function below.  Using that function, it is assumed that sM = sF
fm = "F"
if fm.find("F")>=0:
    sM = 0.0 # FULL FECUNDITY FOR GE
else:
    sM = 1.0 - 0.752 
sM = 1.0
hybridisation = [[[1, 0, 0], [0, 1, 0], [0, 0, 1]], # Aa
[[1, 0, 0], [0, 1, 0], [0, 0, 1]],  # Ac
[[1, 0, 0], [0, 1, 0], [0, 0, 1]]] # Ag
# hybridisation[mM][mF][m] = prob that offspring of species m results from male of species mM and female of species mF.  The current value means offspring is always same as mF.  Must be changed if num_species changes from 3
offspring_modifier = [[[1, 1, 1], [1, 1, 1], [1, 1, 1]], [[1, 1, 1], [1, 1, 1], [1, 1, 1]]]
# offspring_modifier[s][mM][mF] = suppression (if <1, or increased vigour if >1) offspring of sex s that arises from male of species mM and female of species mF
sex_ratio = 0.950 # probability that offspring of wc or cc fathers are male
female_bias = 1. - 0.449 # probability that offspring of (mother wc or cc + father ww) is female
m_w = 0

diffusion_coefficient = 190.0**2 / (math.pi * 7) # 1641.57 m^2/day for male Ac wildtype
diffusion_coeffs = [0.0 for i in range(36)] #[0.0] * 36
windhrs = 9;
if windhrs <= 0:
    advection_fraction = 0
elif windhrs == 2:
    advection_fraction = 2.6E-3
elif windhrs == 9:
    advection_fraction = 6.4E-4
else:
    raise ValueError("windhrs must be <=0, 2 or 9")

######################################################
# Different sex advection scenarios
advectsex = "M";
if advectsex.find("L")>= 0: # low 95% Jeffreys Binomial HPDI from Huestis
    advection_males = 0.0001039736
elif advectsex.find("M")>=0: # advection is 1/22 males estimate from Huestis
    advection_males = 1.0/22.0
elif advectsex.find("H")>= 0: # high 95% Jeffreys Binomial HPDI from Huestis
    advection_males = 0.1645481
elif advectsex.find("O")>= 0: # original 50/50 sex ratio
    advection_males = 0.5 
else:
    raise ValueError("advectsex must be L, M, H or O")

# Example simulation goes from 2022-3
years_simulated = list(range(2022, 2024))
year_begin = 2022
monthdays = [31,28,31,30,31,30,31,31,30,31,30,31] # No leap years

death_rate = [[[0.1 for m in range(3)] for g in range(6)] for s in range(2)]

# Daily survival values are placed into death_rate array
# and then converted to death rate by death rate = -log(daily survival)
                                                                    
for genotype in range(6):
    if genotype == 0 or genotype == 2 or genotype == 5: # WT (ww wr rr) or full mortality (no GE effect)
        sex = 0
        index = genotype * num_species + sex * num_species * 6 # + species 
        death_rate[0][genotype][0] = 0.7958697 # 0.8049269  #male Aa daily survival
        death_rate[0][genotype][1] = 0.7778986   #0.7786237 #male Ac
        death_rate[0][genotype][2] = 0.7670186    #0.7720204 #male Ag
        diffusion_coeffs[index + 0] = 2051.0410 # 1830.3888 # male Aa diffusion coeff
        diffusion_coeffs[index + 1] = diffusion_coefficient # male Ac diffusion coeff
        diffusion_coeffs[index + 2] = 1394.6711 # 1513.4562 # male Ag diffusion coeff
        
        sex = 1
        index = genotype * num_species + sex * num_species * 6 # + species 
        death_rate[1][genotype][0] = 0.850   # 0.8535435    #female Aa
        death_rate[1][genotype][1] = 0.840   # 0.8408630 #female Ac
        death_rate[1][genotype][2] = 0.84   # 0.8393229 #female Ag
        diffusion_coeffs[index + 0] = 9981.4766 # 4037.5528 # female Aa diffusion coeff
        diffusion_coeffs[index + 1] = 7988.4408 # 3621.0799 # female Ac diffusion coeff
        diffusion_coeffs[index + 2] = 6786.7788 # 3338.4999 # female Ag diffusion coeff
    else: # GE (wc cc cr)
        sex = 0
        index = genotype * num_species + sex * num_species * 6 # + species 
        death_rate[0][genotype][0] = 0.7191003 # 0.7289444 #male Aa
        death_rate[0][genotype][1] = 0.6900155 # 0.6923257 #male Ac
        death_rate[0][genotype][2] = 0.6779111    # 0.6810786 #male Ag
        diffusion_coeffs[index + 0] = 829.4782 # 1166.6575 # male Aa diffusion coeff
        diffusion_coeffs[index + 1] = 663.8849 # 1046.3046 # male Ac diffusion coeff
        diffusion_coeffs[index + 2] = 564.0360 # 964.6459 # male Ag diffusion coeff
        
        sex = 1
        index = genotype * num_species + sex * num_species * 6 # + species 
        death_rate[1][genotype][0] = 0.8014861   # 0.8054594 #female Aa
        death_rate[1][genotype][1] = 0.7782918   # 0.7784780 #female Ac
        death_rate[1][genotype][2] = 0.7696894   # 0.7714454 #female Ag
        diffusion_coeffs[index + 0] = 4036.1969  # 2573.5473 # female Aa diffusion coeff
        diffusion_coeffs[index + 1] = 3230.3752  # 2308.0723 # female Ac diffusion coeff
        diffusion_coeffs[index + 2] = 2744.4966  # 2127.9474 # female Ag diffusion coeff
        
    for s in range(2):
        for m in range(3):
            death_rate[s][genotype][m] = -math.log(death_rate[s][genotype][m]) # convert to death rate
# ww wc wr cc cr rr
def reduct(genotype):
    if genotype == 0 or genotype == 2 or genotype == 5: # ww wr rr
        return 1.0
    elif genotype == 1 or genotype == 4: # wc cr
        return 1.0 - sM * hM
    elif genotype == 3: # cc
        return 1.0 - sM
    raise ValueError("invalid genotype")
reduction = [[min(reduct(gM), reduct(gF)) for gM in range(6)] for gF in range(6)]

# add output folder
try:
    os.mkdir(os.path.join(working_dir, output_dir, folder))
except:
    sys.stdout.write("Output folder already exists\n")

# add output folder
try:
    os.mkdir(os.path.join(localdir, folder))
except:
    sys.stdout.write("Output folder already exists\n")

# End parameter definitions
# If num_species != 3, you shouldn't have to change anything below here


######################################################
# Define the ODE
start = timeit.default_timer()
sys.stdout.write("Defining the ODE...")
timestep_size = 1.0
delay_days = int(10.0/timestep_size + 0.5) # number of days in the delay DE
steps = int(round(1./timestep_size))
cell = CellDynamicsMosquitoBH26Delay(num_species = num_species, delay = delay_days, current_index = 0, death_rate = death_rate, competition = competition, emergence_rate = emergence_rate, activity = activity, reduction = reduction, hybridisation = hybridisation, offspring_modifier = offspring_modifier, sex_ratio = sex_ratio, female_bias = female_bias, m_w = m_w, m_c = m_c)
cell.setSmallValue(1.0)
cell.setUseQm(1)


######################################################
# Set two different advection classes (0 for male, 1 for female)
for species in range(num_species):
   for genotype in range(6):
      for sex in range(2):
         ind = species + genotype * num_species + sex * num_species * 6
         cell.setAdvectionClass(ind, sex)

sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")

######################################################
# define a zeroed populations and parameters array
start = timeit.default_timer()
sys.stdout.write("defining the populations and parameters array...")
all_pops = PopulationsAndParameters(g1, cell)
# get coords of introduction cell
sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")

######################################################
# initialise populations
# read a file corresponding to the qm for the first day, to set initial populations
start = timeit.default_timer()
sys.stdout.write("Reading a representative qm to set the initial populations...")
ini_qm_parsers = [0 for i in range(num_species)]
ini_qms = [0 for i in range(num_species)]
first_day = str(years_simulated[0]) + "0101"

# ######################################################
# # load spatial variables

######################################################
# load first spatiotemporal variables

year = year_begin 

temp = SpatialDependence(-1799134, -1299134, 5000.0, 100, 100) 
temp.parse(os.path.join(working_dir, dir_with_qm, "qm.Aa." + str(year_begin) + "." + str(1) + ".csv"), "generic_float", [])
temp.restrictToActive(g1.getGlobalIndex())
ini_qm = temp.getData0()
num_actives = len(ini_qm) 
ini_qms[0] = [ini_qm[i] + 1.e-3 for i in range(num_actives)]

temp = SpatialDependence(-1799134, -1299134, 5000.0, 100, 100) 
temp.parse(os.path.join(working_dir, dir_with_qm, "qm.Ac." + str(year_begin) + "." + str(1) + ".csv"), "generic_float", [])
temp.restrictToActive(g1.getGlobalIndex())
ini_qm = temp.getData0()
ini_qms[1] = [ini_qm[i] + 1.e-3 for i in range(num_actives)]

temp = SpatialDependence(-1799134, -1299134, 5000.0, 100, 100) 
temp.parse(os.path.join(working_dir, dir_with_qm, "qm.Ag." + str(year_begin) + "." + str(1) + ".csv"), "generic_float", [])
temp.restrictToActive(g1.getGlobalIndex())
ini_qm = temp.getData0()
ini_qms[2] = [ini_qm[i] + 1.e-3 for i in range(num_actives)]
sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")
sys.stdout.write("Initialising the populations, with male and female wild-types at the carrying capacity...\n")
start = timeit.default_timer()
num_actives = len(ini_qms[0])

ini_pops = [array.array('f', [ini_qms[m][i] for i in range(num_actives)]) for m in range(num_species)]

# fix with matrix algebra (just using female wildtype)
A = np.array(competition)
#A = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1]]) # no competition
for i in range(num_actives):
    Q = np.array([[1./ini_qms[0][i], 0., 0.],
    [0., 1./ini_qms[1][i], 0.],
    [0., 0., 1./ini_qms[2][i]]])
    B = np.matmul(Q, A)
    b = np.array([emergence_rate[m]/2.0/(1.0 - math.exp(-death_rate[1][0][m])) - 1.0 for m in range(num_species)])
    z = np.linalg.solve(B, b)
    for m in range(num_species):
        ini_pops[m][i] = max(1e-3, z[m]) # make sure pops not below zero

for species in range(num_species):
    for genotype in [0]: # only wild-types
        for sex in range(2):
            for d in range(delay_days + 1):
                index = species + genotype * num_species + sex * num_species * 6 + d * num_species * 6 * 2
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

quants = spatial.getAllQuantities()[:]
num_quantities = len(quants)

for year in years_simulated:
    year_birth_quantities = np.zeros(num_diffusing) #[0.0 for i in range(num_diffusing)]
    year_yyp_quantities = np.zeros(num_diffusing) #[0.0 for i in range(num_diffusing)]
    year_quantities = np.zeros(num_quantities) #[0.0 for i in range(num_quantities)]
    
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
        # define q_m for the month
        sys.stdout.write("Reading and processing q_m for month " + month + "...")
        sys.stdout.flush() 
        start = timeit.default_timer()
        qm_parsers = {}
        qm = {}
        for dayno in range(30):
            doyno = doyno + 1
            if  year > year_begin or monthno >= 1 or dayno >= 30:
                day = str(dayno + 1).zfill(2)
                for speciesno in range(num_species):
                    species = species_list[speciesno]
                    if speciesno == 0:
                        temp = SpatialDependence(-1799134, -1299134, 5000.0, 100, 100) 
                        temp.parse(os.path.join(working_dir, dir_with_qm, "qm.Aa." + str(year) + "." + str(doyno) + ".csv"), "generic_float", [])
                        temp.restrictToActive(g1.getGlobalIndex())
                        ini_qm = temp.getData0()
                        qm[(month, day, species)] = [ini_qm[i] + 1.e-3 for i in range(num_actives)]
                    elif speciesno == 1:
                        temp = SpatialDependence(-1799134, -1299134, 5000.0, 100, 100) 
                        temp.parse(os.path.join(working_dir, dir_with_qm, "qm.Ac." + str(year) + "." + str(doyno) + ".csv"), "generic_float", [])
                        temp.restrictToActive(g1.getGlobalIndex())
                        ini_qm = temp.getData0()
                        qm[(month, day, species)] = [ini_qm[i] + 1.e-3 for i in range(num_actives)]
                    else:
                        temp = SpatialDependence(-1799134, -1299134, 5000.0, 100, 100) 
                        temp.parse(os.path.join(working_dir, dir_with_qm, "qm.Ag." + str(year) + "." + str(doyno) + ".csv"), "generic_float", [])
                        temp.restrictToActive(g1.getGlobalIndex())
                        ini_qm = temp.getData0()
                        qm[(month, day, species)] = [ini_qm[i] + 1.e-3 for i in range(num_actives)]
                        
                    qm[(month, day, species)] = array.array('f', qm[(month, day, species)])

        sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")
        sys.stdout.flush() 

        ######################################################
        # evolve over the month
        for dayno in range(monthdays[monthno]):
            #if year > year_begin or monthno >= 1 or dayno >= 10: # 10 delay days (start simulation on the 11th January)
                day = str(dayno + 1).zfill(2)
                sys.stdout.write("    Day " + day + "\n")
                for m in range(num_species):
                    all_pops.setOverActiveGrid(cell.getNumberOfPopulations() + m, qm[(month, day, species_list[m])])
                # Ac ww male
                #index = 1 + 0 * num_species + 0 * num_species * 6 + cell.getCurrentIndex() * num_species * 6 * 2
                for step in range(steps):
                    spatial.evolveCells(timestep_size)
                    cell.incrementCurrentIndex() # needed for Delay DE!
                    sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")
                    if windhrs >= 0:
                        spatial.diffuseVarying(timestep_size, array.array('f', diffusion_coeffs))
                    sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")
                    
                #spatial.advect(advection_fraction, wind[(month, day)])
                # makes 21 females advect for every 1 male while keeping overall probability the same
                if windhrs > 0 and year >= year_begin:
                    spatial.advect(array.array('f', [2.0 * advection_males * advection_fraction, 2.0 * (1.0 - advection_males) * advection_fraction]), wind[(month, day)])
                sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")
                
                curr_birth_quantities = spatial.getBirthQuantities()[:]
                curr_yyp_quantities = spatial.getYypQuantities()[:]
                year_birth_quantities = year_birth_quantities + np.array(curr_birth_quantities) #[year_birth_quantities[i] + curr_birth_quantities[i] for i in range(num_diffusing)]
                year_yyp_quantities = year_yyp_quantities + np.array(curr_yyp_quantities) #[year_birth_quantities[i] + curr_birth_quantities[i] for i in range(num_diffusing)]
                
                quants = spatial.getAllQuantities()[:]
                year_quantities = year_quantities + np.array(quants) #[year_quantities[i] + quants[i] for i in range(num_quantities)]
                                
                sys.stdout.write(" " + str(timeit.default_timer() - start) + "s\n")
            
                sys.stdout.flush() # every month
        
    
        sys.stdout.write("CPU time taken = " + str(timeit.default_timer() - start) + "s\n")
    
    # every year
    spatial.setBirthQuantities(year_birth_quantities.tolist())
    spatial.setYypQuantities(year_yyp_quantities.tolist())
    
               
    quants = spatial.getAllQuantities()[:]
    spatial.setAllQuantities(year_quantities.tolist())
    genotype = 0
    for species in range(num_species):
        for sex in range(2):
            ind = species + genotype * num_species + sex * num_species * 6
            m_str = species_list[species] + "_" + genotype_list[genotype] + "_" + sex_list[sex]
            spatial.outputCSV(os.path.join(working_dir, output_dir, folder, "popsum_" + str(windhrs) + "_" + m_str + "_" + str(year) + "_" + month + "_" + day + ".csv"), ind, "0", "")
    spatial.setAllQuantities(quants.tolist())
sys.exit(0)
