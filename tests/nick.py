import os
import sys
import array

# so we can find our ../code no matter how we are called
findbin = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(findbin + "/../code")

from cellDynamics import CellDynamicsMosquito23

cell = CellDynamicsMosquito23()
dt=1
cell.setMuLarvae(0.1)
cell.setMuAdult(0.1)
cell.setFecundity(0.9)
cell.setAgingRate(0.1)
cell.setNumAges(1)
cell.setAccuracy(0.95)

initial_condition = [0.1] * cell.getNumberOfPopulations() + [9.0 / 7.0]
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
