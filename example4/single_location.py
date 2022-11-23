# Following are the parameters and input-files of this model
species_list = ["ag", "ac", "ar"] # if you change to !=3 subspecies, a few other things below will need to be changed, as noted below.  If you change these names, then search throughout this file for "species_list" to see what else you might need to change
sex_list = ["male", "female"] # must be 2 sexes
genotype_list = ["ww", "wc", "wr", "cc", "cr", "rr"] # must be 6 genotypes
eqm_wild_pops = [140000, 140000, 140000] # for the 3 species.  Must be changed if num_species changes from 3
num_species = len(species_list)
delay_days = 10 # number of days in the delay DE
death_rate_ww = [0.1, 0.1, 0.1] # death rate of wild-types of each species, measured in 1/day.  Below we assume the other genotypes have the same death rate: if a bad assumption then just modify death_rate variable below.  Must be changed if num_species changes from 3
#competition = [[1, 0.4, 0.4], [0.4, 1, 0.4], [0.4, 0.4, 1]] # competition between subspecies.  Must be changed if num_species changes from 3
#competition = [[1, 0, 0], [0, 1, 0], [0, 0, 1]] # competition between subspecies.  Must be changed if num_species changes from 3
competition = [[1, 0.717, 0.562], [0.717, 1, 0.317], [0.562, 0.317, 1]]
emergence_rate = [9.0, 9.0, 9.0] # emergence rate for each species.  Must be changed if num_species changes from 3
activity = [[0.9982032342, 0.0008983829, 0.0008983829], 
[0.0008974781, 0.9971978740, 0.0019046479], 
[0.0008974781, 0.0019046479, 0.9971978740]] # activity[female_of_species1][male_of_species2] is activity level in the proportionate mixing.  Must be changed if num_species changes from 3
hM = 0.0 # reduction R = (1 - s * h) * (1 - s * h) for genotype = wc, for instance: see the reduct function below.  Using that function, it is assumed that hM = hF
sM = 0.0 # reduction R = (1 - s * h) * (1 - s * h) for genotype = wc, for instance: see the reduct function below.  Using that function, it is assumed that sM = sF
hybridisation = [[[1, 0, 0], [0, 1, 0], [0, 0, 1]], [[1, 0, 0], [0, 1, 0], [0, 0, 1]], [[1, 0, 0], [0, 1, 0], [0, 0, 1]]] # hybridisation[mM][mF][m] = prob that offspring of species m results from male of species mM and female of species mF.  The current value means offspring is always same as mF.  Must be changed if num_species changes from 3
offspring_modifier = [[[1, 1, 1], [1, 1, 1], [1, 1, 1]], [[1, 1, 1], [1, 1, 1], [1, 1, 1]]] # offspring_modifier[s][mM][mF] = suppression (if <1, or increased vigour if >1) of sppring of sex s that arises from male of species mM and female of species mF
sex_ratio = 0.950 # probability that offspring of wc or cc fathers are male
female_bias = 1. - 0.449 # probability that offspring of (mother wc or cc + father ww) is female
m_w = 1E-6
m_c = 1E-6
small_value = 1E-5 # when populations, carrying capacity, etc get smaller than this, set them to zero
intro = 1E4 # number of AC wc mosquitoes introduced


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
# End parameter definitions
# If num_species != 3, you shouldn't have to change anything below here


######################################################
# Import libraries
import os
import sys
import array
import matplotlib.pyplot as plt
findbin = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(findbin + "/../../mozzie/code")
from cellDynamics import CellDynamicsMosquitoBH26Delay


######################################################
# Define the ODE
# Note that m_w = 0 at this stage, since we don't want any resistant types being produced
cell = CellDynamicsMosquitoBH26Delay(num_species = num_species, delay = delay_days, current_index = 0, death_rate = death_rate, competition = competition, emergence_rate = emergence_rate, activity = activity, reduction = reduction, hybridisation = hybridisation, offspring_modifier = offspring_modifier, sex_ratio = sex_ratio, female_bias = female_bias, m_w = 0.0, m_c = m_c)
cell.setSmallValue(small_value)

 
######################################################
# Define initial condition, which is wild-only
initial_condition = [0 for i in range(cell.getNumberOfPopulations() + cell.getNumberOfParameters())]
for species in range(num_species):
    for genotype in [0]: # only wild-types
        for sex in range(2):
            for d in range(delay_days + 1):
                 index = species + genotype * num_species + sex * num_species * 6 + d * num_species * 6 * 2
                 initial_condition[index] = eqm_wild_pops[species]
sys.stdout.write(str(initial_condition) + "\n")
qm = cell.calcQm(array.array('f', initial_condition)) # calculate qm for this equilibrium
for m in range(num_species):
    initial_condition[cell.getNumberOfPopulations() + m] = qm[m] # set the qm values, which live at the end of the list
    sys.stdout.write("qm[" + str(m) + "] = " + str(qm[m]) + "\n")


######################################################
# Evolve to find the true steady-state with the given competition
pap = array.array('f', initial_condition)
dt = 1
for day in range(1000):
    sys.stdout.write("day = " + str(day) + ". (ww) = ")
    for species in [1]: # AC
        for genotype in [0]: # ww and wc
            for sex in [0]: # males
                for d in [cell.getCurrentIndex()]:
                    ind = species + genotype * num_species + sex * num_species * 6 + d * num_species * 6 * 2
                    sys.stdout.write(str(round(pap[ind], 0)) + " ")
    sys.stdout.write("\n")
    cell.evolve(dt, pap)
    cell.incrementCurrentIndex() # needed for Delay DE!

######################################################
# Define the ODE with m_w != 0 so resistant types can be created
cell = CellDynamicsMosquitoBH26Delay(num_species = num_species, delay = delay_days, current_index = 0, death_rate = death_rate, competition = competition, emergence_rate = emergence_rate, activity = activity, reduction = reduction, hybridisation = hybridisation, offspring_modifier = offspring_modifier, sex_ratio = sex_ratio, female_bias = female_bias, m_w = m_w, m_c = m_c)
cell.setSmallValue(small_value)


######################################################
# Introduce modified mosquitoes as current adults
# This means they will breed instantly, but their offspring won't been seen in adult populations for "delay" days
genotype = 1 # genotype = wc
sex = 0 # sex = male
species = 1 # species = AC
num = intro
pap[species + genotype * num_species + sex * num_species * 6 + cell.getCurrentIndex() * num_species * 6 * 2] = num


######################################################
# Evolve
results = []
for day in range(365):
    # output ww and wc, species=AC adult populations at the start of this day
    sys.stdout.write("day = " + str(day) + ". (ww, wc) = ")
    for species in [1]: # AC
        for genotype in [0, 1]: # ww and wc
            for sex in [0]: # males
                for d in [cell.getCurrentIndex()]:
                    ind = species + genotype * num_species + sex * num_species * 6 + d * num_species * 6 * 2
                    sys.stdout.write(str(round(pap[ind], 0)) + " ")
    sys.stdout.write("\n")
    results.append([day, pap[1 + 0 * num_species + cell.getCurrentIndex() * num_species * 6 * 2], pap[1 + 1 * num_species + cell.getCurrentIndex() * num_species * 6 * 2]])
    cell.evolve(dt, pap)
    cell.incrementCurrentIndex() # needed for Delay DE!

plt.figure()
plt.plot([x[0] for x in results], [x[1] for x in results], label = 'ww')
plt.plot([x[0] for x in results], [x[2] for x in results], label = 'wc')
plt.legend()
plt.xlabel("Time (day)")
plt.ylabel("Population")
plt.title("Male population, with hM = " + str(hM) + " sm = " + str(sM))
plt.grid()
plt.show()
plt.close()

