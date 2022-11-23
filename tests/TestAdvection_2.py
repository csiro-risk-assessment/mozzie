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

class TestAdvection_1(unittest.TestCase):

   def setUp(self):
      # This is the grid and wind used in TestWind, so we can be sure it's correct
      self.g1 = Grid(1.0, 2.0, 3.0, 4, 3)
      self.w1 = Wind(os.path.join(findbin, "wind1.csv"), os.path.join(findbin, "wind1_processed.csv"), [[1.0, 0.3], [2.0, 0.7]], self.g1)
      self.w1.parseRawFile()
      self.g2 = Grid(1.0, 2.0, 3.0, 4, 3)
      self.g2.setActiveAndInactive(os.path.join(findbin, "inactive_active.csv"))
      self.w2 = Wind(os.path.join(findbin, "wind1.csv"), os.path.join(findbin, "wind2_processed.csv"), [[1.0, 0.3], [2.0, 0.7]], self.g2)
      self.w2.parseRawFile()

      self.cell = CellDynamicsStatic15_9_3_2()

   def testAdvect4(self):
      pop = list(range(self.cell.getNumberOfPopulations() + self.cell.getNumberOfParameters()))

      # start with a population at active index=1 and see it go to index=0
      all_quantities = PopulationsAndParameters(self.g2, self.cell)
      all_quantities.setPopulationAndParameters(1, pop)
      spatial = SpatialDynamics(self.g2, all_quantities)
      spatial.advect(array.array('f', [1.0]), self.w2)
      spatial.outputCSV(os.path.join(findbin, "2D_advection_out_4_1.csv"), 4, "0", "")
      with open(os.path.join(findbin, "2D_advection_out_4_1.csv")) as f:
         data = f.readlines()
      self.assertTrue(arrayfuzzyequal([float(d) for d in data[2].strip().split(",")], [4.0, 0, 0.0, 0], 1E-5))
      for row in [3, 4]:
         self.assertTrue(arrayfuzzyequal([float(d) for d in data[row].strip().split(",")], [0.0] * 4, 1E-5))

   def testAdvect5(self):
      pop = list(range(self.cell.getNumberOfPopulations() + self.cell.getNumberOfParameters()))

      # start with a population at active index=2 and see it disappear
      all_quantities = PopulationsAndParameters(self.g2, self.cell)
      all_quantities.setPopulationAndParameters(2, pop)
      spatial = SpatialDynamics(self.g2, all_quantities)
      spatial.advect(array.array('f', [1.0]), self.w2)
      spatial.outputCSV(os.path.join(findbin, "2D_advection_out_5_1.csv"), 4, "0", "")
      with open(os.path.join(findbin, "2D_advection_out_5_1.csv")) as f:
         data = f.readlines()
      for row in [2, 3, 4]:
         self.assertTrue(arrayfuzzyequal([float(d) for d in data[row].strip().split(",")], [0.0] * 4, 1E-5))

   def testAdvect6(self):
      pop5 = list(range(self.cell.getNumberOfPopulations() + self.cell.getNumberOfParameters()))
      pop6 = [1.0] * (self.cell.getNumberOfPopulations() + self.cell.getNumberOfParameters())

      # start with a population at active index=5 and active index = 6 and see them advect
      all_quantities = PopulationsAndParameters(self.g2, self.cell)
      all_quantities.setPopulationAndParameters(5, pop5)
      all_quantities.setPopulationAndParameters(6, pop6)
      spatial = SpatialDynamics(self.g2, all_quantities)
      spatial.advect(array.array('f', [0.6]), self.w2)
      spatial.outputCSV(os.path.join(findbin, "2D_advection_out_6_1.csv"), 4, "0", "")
      with open(os.path.join(findbin, "2D_advection_out_6_1.csv")) as f:
         data = f.readlines()
      self.assertTrue(arrayfuzzyequal([float(d) for d in data[2].strip().split(",")], [0, 0, 0, 0], 1E-5))
      self.assertTrue(arrayfuzzyequal([float(d) for d in data[3].strip().split(",")], [0, 0, 0, 0], 1E-5))
      self.assertTrue(arrayfuzzyequal([float(d) for d in data[4].strip().split(",")], [0, 0, 3.46, 1.54], 1E-5))

if __name__ == '__main__':
   unittest.main()

