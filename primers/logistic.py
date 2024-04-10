######################################################
# Example simulation with
# - a single location (no spatial spread)
# - 'lifecycle' defined by the logistic equation with growth rate = 0.1/day and carrying capacity = 1E3
# The logistic equation for population P is
# dP/dt = r P (1 - P / K)
# where t is time, r is the growth rate and K the carrying capacity.
# Numerical values of the parameters are r = 0.01 and K = 1E3.
#
######################################################
# Import libraries
import os
import sys
import math

# so we can find our ../code no matter how we are called
findbin = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(findbin + "/../code")

from grid import Grid
from cellDynamics import CellDynamicsLogistic1_1
from spatialDynamics import SpatialDynamics
from populationsAndParameters import PopulationsAndParameters


######################################################
# setup the grid (just one spatial location)
g1 = Grid(0, 0, 1.0, 1, 1)

######################################################
# define the lifecycle dynamics to be governed by the
# logistic equation.
cell = CellDynamicsLogistic1_1()

######################################################
# Using CellDynamicsLogistic1_1() means that at each
# spatial location there is:
# - 1 population
# - 1 parameter, which is the carrying capacity
# Define the populations and parameters over the grid:
all_pops = PopulationsAndParameters(g1, cell)
# So, at each spatial location (there is only one of
# location in this example) there will be two "populations
# and parameters": the first being the population, the
# second being the carrying capacity.
# Insert an initial condition of 10 individuals
# and a carrying capacity of 1E3:
all_pops.setPopulationAndParametersFromXY(0, 0, [10, 1E3])

######################################################
# put these populations into the spatial structure
spatial = SpatialDynamics(g1, all_pops)

######################################################
# Evolve in time (ie, solve the logistic equation)
mozzie_result = [(0, all_pops.getQuantities()[0])]  # index [0] is the population
dt = 5 # time step size
t = 0
for i in range(200):
   spatial.evolveCells(dt)
   t += dt
   mozzie_result.append((t, all_pops.getQuantities()[0]))

######################################################
# For comparison, use an explicit solver
# The solution of (P - Pold) / dt = r P (1 - P / K)
# is P = (r dt - 1) + sqrt( (1 - r dt)^2 + 4 Pold r dt / K ) ) / (2 r dt / K)
implicit_result = [(0, 10.0)]
dt = 5 # time step size
t = 0
rdt = 0.01 * dt
for i in range(200):
   pold = implicit_result[-1][-1]
   p = ((rdt - 1) + math.sqrt( (1 - rdt) * (1 - rdt) + 4 * pold * rdt / 1E3 )) / (2 * rdt / 1E3)
   t += dt
   implicit_result.append((t, p))

######################################################
# Plot the comparison
import matplotlib.pyplot as plt
plt.figure()
plt.plot([x[0] for x in mozzie_result], [x[1] for x in mozzie_result], label = 'Mozzie (explicit time-stepping)')
plt.plot([x[0] for x in implicit_result], [x[1] for x in implicit_result], label = 'Implicit time-stepping')
plt.grid()
plt.legend()
plt.xlabel("Time")
plt.ylabel("Population")
plt.title("Logistic growth example (r = 0.01, K = 1000)")
plt.savefig("logistic.pdf", bbox_inches = 'tight')
plt.show()

sys.exit(0)
