import os
import sys
import unittest

# so we can find our ../code no matter how we are called
findbin = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(findbin + "/../code")

from grid import Grid
from cellDynamics import CellDynamicsLogistic1_1
from spatialDynamics import SpatialDynamics
from populationsAndParameters import PopulationsAndParameters

delete_output = True

def arrayfuzzyequal(a, b, eps):
   return (len(a) == len(b)) and all([(a[i] > b[i] - eps and a[i] < b[i] + eps) for i in range(min(len(a), len(b)))])

class TestEvolveCells(unittest.TestCase):

   def setUp(self):
      self.grid = Grid(0, 0, 1, 6, 1)

      self.cell = CellDynamicsLogistic1_1()

      # all the population and parameters
      self.all_quantities = PopulationsAndParameters(self.grid, self.cell)
      # include some interesting population numbers and carrying capacities (latter can't be zero because of logistic equation)
      self.all_quantities.setPopulationAndParameters(0, [0, 1])
      self.all_quantities.setPopulationAndParameters(1, [1, 2])
      self.all_quantities.setPopulationAndParameters(2, [0, 2])
      self.all_quantities.setPopulationAndParameters(3, [1, -1])
      self.all_quantities.setPopulationAndParameters(4, [2, 2])
      self.all_quantities.setPopulationAndParameters(5, [4, 1])

      self.spatial = SpatialDynamics(self.grid, self.all_quantities)


   def testEvolveCells1(self):
      # evolve with timestep = 12
      self.spatial.evolveCells(12)
      fn = os.path.join(findbin, "2D_evolution_out_1.csv")
      self.spatial.outputCSV(fn, 0, "0", "")
      with open(fn) as f:
         data = f.readlines()
      self.assertTrue(arrayfuzzyequal([float(d) for d in data[2].strip().split(",")], [0.0, 1.06, 0.0, 1.24, 2.0, 2.56], 1E-5))
      if os.path.isfile(fn) and delete_output: os.remove(fn)

      fn = os.path.join(findbin, "2D_evolution_out_2.csv")
      self.spatial.outputCSV(fn, 1, "0", "")
      with open(fn) as f:
         data = f.readlines()
      self.assertTrue(arrayfuzzyequal([float(d) for d in data[2].strip().split(",")], [1.0, 2.0, 2.0, -1.0, 2.0, 1.0], 1E-5))
      if os.path.isfile(fn) and delete_output: os.remove(fn)

      # evolve with timestep = 3
      self.spatial.evolveCells(3)
      fn = os.path.join(findbin, "2D_evolution_out_3.csv")
      self.spatial.outputCSV(fn, 0, "0", "")
      with open(fn) as f:
         data = f.readlines()
      self.assertTrue(arrayfuzzyequal([float(d) for d in data[2].strip().split(",")], [0.0, 1.0749459, 0.0, 1.323328, 2.0, 2.44019198], 1E-5))
      if os.path.isfile(fn) and delete_output: os.remove(fn)

      fn = os.path.join(findbin, "2D_evolution_out_4.csv")
      self.spatial.outputCSV(fn, 1, "0", "")
      with open(fn) as f:
         data = f.readlines()
      self.assertTrue(arrayfuzzyequal([float(d) for d in data[2].strip().split(",")], [1.0, 2.0, 2.0, -1.0, 2.0, 1.0], 1E-5))
      if os.path.isfile(fn) and delete_output: os.remove(fn)


if __name__ == '__main__':
   unittest.main()

