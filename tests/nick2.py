import os
import sys
import array

# so we can find our ../code no matter how we are called
findbin = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(findbin + "/../code")

from cellDynamics import CellDynamicsMosquito26

cell = CellDynamicsMosquito26()
dt=1.

cell.setTimeIntegrationMethod("runge_kutta4")
cell.setMinCarryingCapacity(1.0)
cell.setZeroCutoff(1.0)
cell.setMinimumDt(1.0)

## mother determines offspring subspecies
cell.setHybridisationRate(0, 1, 1, 1.0) 
cell.setHybridisationRate(0, 1, 0, 0.0) 
cell.setHybridisationRate(1, 0, 0, 1.0)
cell.setHybridisationRate(1, 0, 1, 0.0) 

#cell.setAccuracy(0.95)

#initial_condition = [0.1] * cell.getNumberOfPopulations() + [9.0 / 7.0]
initial_condition = [9.861140251159668, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.861140251159668, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 135.16949462890625, 4.195792198181152, 0.0036913766525685787, 0.0, 0.0, 0.0, 3.2018029116898106e-08, 0.0, 0.0, 0.0, 0.0, 0.0, 135.16949462890625, 4.195792198181152, 0.0036913766525685787, 0.0, 0.0, 0.0, 3.2018029116898106e-08, 0.0, 0.0, 0.0, 0.0, 0.0, 19.37061309814453, 8.370390892028809]
pap = array.array('f', initial_condition)
cell.evolve(dt, pap)

## If we want to run for longer
#import numpy as np
#dt = 1
#res = np.zeros((10000, 6))
#for i in range(10000):
#    cell.evolve(dt, pap)
#    res[i,:] = pap[0:6]
#import matplotlib.pyplot as plt
#fig, ax = plt.subplots()
#for i in range(6):
#    ax.plot(np.arange(0.01, 100.0001, 0.01), res[:,i], color = ((i%3)/2.2, ((i//3)%2)*0.9, (i//6)*0.9), label = str(i))
#ax.legend()
#plt.savefig('foo.png')

print("and the final population is....", pap)

# expected from Fortran/R code:
# 0.10800000 0.12600000 0.10800000 0.10231579 0.10326316 0.09094737
