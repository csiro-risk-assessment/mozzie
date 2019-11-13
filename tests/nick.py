import os
import sys
import array

# so we can find our ../code no matter how we are called
findbin = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(findbin + "/../code")

from cellDynamics import CellDynamicsMosquito23

cell = CellDynamicsMosquito23()
dt = 50.0
cell.setMuLarvae(0.1)
cell.setMuAdult(0.1)
cell.setFecundity(0.9)
cell.setAgingRate(0.1)
cell.setNumAges(1)
cell.setAccuracy(0.95)

initial_condition = [0.1] * cell.getNumberOfPopulations() + [9.0 / 7.0]
pap = array.array('f', initial_condition)
cell.evolve(dt, pap)
print("and the final population is....", pap)
