import os
import sys
import matplotlib.pyplot as plt
import array
import numpy as np
findbin = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(findbin + "/../code")

from cellDynamics import CellDynamicsMosquitoBH26Delay

print("Solves mosquito evolution according to BH26Delay for a single species with wild-types only, using various different timesteps")

num_species = 1
death_rate = 1.0
competition = 1.0
emergence_rate = 4.0
activity = 1.0
reduction = 1.0
hybridisation = 1.0
offspring_modifier = 1.0
sex_ratio = 0.5
female_bias = 0.5
m_w = 0.0 # so only wild-type mosquitos
m_c = 0.0
q_m = 10.0

num_sexes = 2
HORXprimeM = hybridisation * (1.0 / num_sexes) * reduction
steady_state = (1.0 / (num_sexes * competition) / (HORXprimeM * emergence_rate)) * q_m * (HORXprimeM * offspring_modifier * emergence_rate / death_rate - 1.0)
print("Male population at steady state =", steady_state)

end_time = 100.0

plt.figure()
plt.plot([0, end_time], [steady_state, steady_state], 'k:', label = "steady-state")

for dt in [1, 0.5, 0.001]:
    delay = int(10.0 / dt + 0.5)
    wild = CellDynamicsMosquitoBH26Delay(num_species = num_species,
                                         delay = delay,
                                         current_index = 0,
                                         death_rate = [[[death_rate] for i in range(6)], [[death_rate] for i in range(6)]],
                                         competition = [[competition]],
                                         emergence_rate = [emergence_rate],
                                         activity = [[activity]],
                                         reduction = [[reduction for i in range(6)] for j in range(6)],
                                         hybridisation = [[[hybridisation]]],
                                         offspring_modifier = [[[offspring_modifier]], [[offspring_modifier]]],
                                         sex_ratio = sex_ratio,
                                         female_bias = female_bias,
                                         m_w = m_w,
                                         m_c = m_c)

    # initialise with a small population of wild-types
    initial_condition = [0.1 * steady_state for i in range(wild.getNumberOfPopulations())] + [q_m]
    for delay in range(wild.getDelay() + 1):
        for species in range(num_species):
            for genotype in range(6):
               for sex in range(2):
                  ind = species + genotype * num_species + sex * num_species * 6 + delay * num_species * 6 * 2
                  if genotype != 0:
                     initial_condition[ind] = 0.0 # zero non-wildtype
    pap = array.array('f', initial_condition)

    sim_times = []
    sim_pops = []
    for t in np.arange(0, end_time, dt):
        wild.evolve(dt, pap)
        wild.incrementCurrentIndex()
        sim_times.append(t)
        delay, species, genotype, sex = (wild.getCurrentIndex(), 0, 0, 0)
        sim_pops.append(pap[species + genotype * num_species + sex * num_species * 6 + delay * num_species * 6 * 2])

    plt.plot(sim_times, sim_pops, label = "dt = " + str(dt))

plt.legend()
plt.grid()
plt.xlim([0, plt.xlim()[1]])
plt.xlabel("time (days)")
plt.ylabel("Male population")
plt.title("BH evolution using Andy method with different timestep sizes")
plt.savefig("wild_1_species.png", bbox_inches = "tight")
plt.show()
plt.close()

sys.exit(0)
