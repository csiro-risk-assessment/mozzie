import os
import sys
import array
import unittest

# so we can find our ../code no matter how we are called
findbin = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(findbin + "/../code")

from grid import Grid
from wind import Wind
from cellDynamics import CellDynamicsStatic15_9_3_2
from spatialDynamics import SpatialDynamics
from populationsAndParameters import PopulationsAndParameters

def arrayfuzzyequal(a, b, eps):
   return (len(a) == len(b)) and all([(a[i] > b[i] - eps and a[i] < b[i] + eps) for i in range(min(len(a), len(b)))])

class TestAdvection_3(unittest.TestCase):

   def setUp(self):
      # This is the grid and wind used in TestWind, so we can be sure it's correct
      self.g1 = Grid(1.0, 2.0, 3.0, 4, 3)
      self.w1 = Wind(os.path.join(findbin, "wind1.csv"), os.path.join(findbin, "wind1_processed.csv"), [[1.0, 0.3], [2.0, 0.7]], self.g1)
      self.w1.parseRawFile()
      self.g2 = Grid(1.0, 2.0, 3.0, 4, 3)
      self.g2.setActiveAndInactive(os.path.join(findbin, "inactive_active.csv"))
      self.w2 = Wind(os.path.join(findbin, "wind1.csv"), os.path.join(findbin, "wind2_processed.csv"), [[1.0, 0.3], [2.0, 0.7]], self.g2)

      self.cell = CellDynamicsStatic15_9_3_2()
      self.cell.setAdvectionClass(4, 2)
      self.cell.setAdvectionClass(10, 1)

   def testAdvect7(self):
      pop = list(range(self.cell.getNumberOfPopulations() + self.cell.getNumberOfParameters()))

      # start with a population at index=0 and see it remain there
      all_quantities = PopulationsAndParameters(self.g1, self.cell)
      all_quantities.setPopulationAndParameters(0, pop)
      spatial = SpatialDynamics(self.g1, all_quantities)
      spatial.advect(array.array('f', [1.0, 0.125, 0.5]), self.w1)
      spatial.outputCSV(os.path.join(findbin, "2D_advection_out_0_1.csv"), 4, "0", "")
      with open(os.path.join(findbin, "2D_advection_out_0_1.csv")) as f:
         data = f.readlines()
      self.assertTrue(arrayfuzzyequal([float(d) for d in data[2].strip().split(",")], [4.0, 0, 0, 0], 1E-5))
      for row in [3, 4]:
         self.assertTrue(arrayfuzzyequal([float(d) for d in data[row].strip().split(",")], [0.0] * 4, 1E-5))

   def testAdvect8(self):
      pop = list(range(1, self.cell.getNumberOfPopulations() + self.cell.getNumberOfParameters() + 1))
      pop[0] = 1
      pop[10] = 8
      pop[4] = 4
      # start with a population at index=1 and see:
      # class=0 (population_index=0, population quantity=1) all go to index=5
      # class=1 (population_index=10, population quantity=8) 0.5 go to index=5
      # class=2 (population_index=4, population quantity = 4) 0.25 go to index=5
      all_quantities = PopulationsAndParameters(self.g1, self.cell)
      all_quantities.setPopulationAndParameters(1, pop)
      spatial = SpatialDynamics(self.g1, all_quantities)
      spatial.advect(array.array('f', [1.0, 0.5, 0.25]), self.w1)

      # class=0 (index=0, initial_population=1)
      spatial.outputCSV(os.path.join(findbin, "2D_advection_out_8_0.csv"), 0, "0", "")
      with open(os.path.join(findbin, "2D_advection_out_8_0.csv")) as f:
         data = f.readlines()
      self.assertTrue(arrayfuzzyequal([float(d) for d in data[3].strip().split(",")], [0, 1.0, 0, 0], 1E-5))
      for row in [2, 4]:
         self.assertTrue(arrayfuzzyequal([float(d) for d in data[row].strip().split(",")], [0.0] * 4, 1E-5))

      # class=1 (index=10, initial population=8)
      spatial.outputCSV(os.path.join(findbin, "2D_advection_out_8_1.csv"), 10, "0", "")
      with open(os.path.join(findbin, "2D_advection_out_8_1.csv")) as f:
         data = f.readlines()
      self.assertTrue(arrayfuzzyequal([float(d) for d in data[2].strip().split(",")], [0, 4.0, 0, 0], 1E-5))
      self.assertTrue(arrayfuzzyequal([float(d) for d in data[3].strip().split(",")], [0, 4.0, 0, 0], 1E-5))
      self.assertTrue(arrayfuzzyequal([float(d) for d in data[4].strip().split(",")], [0, 0, 0, 0], 1E-5))
         
      # class=2 (index=4, initial population=4)
      spatial.outputCSV(os.path.join(findbin, "2D_advection_out_8_2.csv"), 4, "0", "")
      with open(os.path.join(findbin, "2D_advection_out_8_2.csv")) as f:
         data = f.readlines()
      self.assertTrue(arrayfuzzyequal([float(d) for d in data[2].strip().split(",")], [0, 3.0, 0, 0], 1E-5))
      self.assertTrue(arrayfuzzyequal([float(d) for d in data[3].strip().split(",")], [0, 1.0, 0, 0], 1E-5))
      self.assertTrue(arrayfuzzyequal([float(d) for d in data[4].strip().split(",")], [0, 0, 0, 0], 1E-5))
         
if __name__ == '__main__':
   unittest.main()

