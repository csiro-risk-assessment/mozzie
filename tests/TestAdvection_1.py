import os
import sys
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
   return all([(a[i] > b[i] - eps and a[i] < b[i] + eps) for i in range(0, len(a))])

class TestAdvection_1(unittest.TestCase):

   def setUp(self):
      # This is the grid and wind used in TestWind, so we can be sure it's correct
      self.g1 = Grid(1.0, 2.0, 3.0, 4, 3)
      self.w1 = Wind(os.path.join(findbin, "wind1.csv"), os.path.join(findbin, "wind1_processed.csv"), [[1.0, 0.3], [2.0, 0.7]], self.g1)
      self.w1.parseRawFile()
      self.g2 = Grid(1.0, 2.0, 3.0, 4, 3)
      self.g2.setActiveAndInactive(os.path.join(findbin, "inactive_active.csv"))
      self.w2 = Wind(os.path.join(findbin, "wind1.csv"), os.path.join(findbin, "wind2_processed.csv"), [[1.0, 0.3], [2.0, 0.7]], self.g2)

      self.cell = CellDynamicsStatic15_9_3_2()

   def testUnprocessedWind(self):
      w_not_processed = Wind(os.path.join(findbin, "wind1.csv"), os.path.join(findbin, "wind1_processed.csv"), [[1.0, 0.3], [2.0, 0.7]], self.g1)
      pop = list(range(self.cell.getNumberOfPopulations() + self.cell.getNumberOfParameters()))
      all_quantities = PopulationsAndParameters(self.g1, self.cell)
      all_quantities.setPopulationAndParameters(0, pop)
      spatial = SpatialDynamics(self.g1, all_quantities)
      with self.assertRaises(ValueError) as the_err:
         spatial.advect(1.0, w_not_processed)
      self.assertEqual(str(the_err.exception), "Wind must have been processed before being used")


   def testAdvect0(self):
      pop = list(range(self.cell.getNumberOfPopulations() + self.cell.getNumberOfParameters()))

      # start with a population at index=0 and see it remain there
      all_quantities = PopulationsAndParameters(self.g1, self.cell)
      all_quantities.setPopulationAndParameters(0, pop)
      spatial = SpatialDynamics(self.g1, all_quantities)
      spatial.advect(1.0, self.w1)
      spatial.outputCSV(os.path.join(findbin, "2D_advection_out_0_1.csv"), 4)
      with open(os.path.join(findbin, "2D_advection_out_0_1.csv")) as f:
         data = f.readlines()
      self.assertTrue(arrayfuzzyequal([float(d) for d in data[0].strip().split(",")], [4.0, 0, 0, 0], 1E-5))
      for row in [1, 2]:
         self.assertTrue(arrayfuzzyequal([float(d) for d in data[row].strip().split(",")], [0.0] * 4, 1E-5))

   def testAdvect1(self):
      pop = list(range(self.cell.getNumberOfPopulations() + self.cell.getNumberOfParameters()))

      # start with a population at index=1 and see it all go to index=5
      all_quantities = PopulationsAndParameters(self.g1, self.cell)
      all_quantities.setPopulationAndParameters(1, pop)
      spatial = SpatialDynamics(self.g1, all_quantities)
      spatial.advect(1.0, self.w1)
      spatial.outputCSV(os.path.join(findbin, "2D_advection_out_1_1.csv"), 4)
      with open(os.path.join(findbin, "2D_advection_out_1_1.csv")) as f:
         data = f.readlines()
      self.assertTrue(arrayfuzzyequal([float(d) for d in data[1].strip().split(",")], [0, 4.0, 0, 0], 1E-5))
      for row in [0, 2]:
         self.assertTrue(arrayfuzzyequal([float(d) for d in data[row].strip().split(",")], [0.0] * 4, 1E-5))

   def testAdvect2(self):
      pop = list(range(self.cell.getNumberOfPopulations() + self.cell.getNumberOfParameters()))

      # start with a population at index=6 and see 30% of it go to index 2 and 70% to index 0
      # use advection_fraction = 5 for more testing
      all_quantities = PopulationsAndParameters(self.g1, self.cell)
      all_quantities.setPopulationAndParameters(6, pop)
      spatial = SpatialDynamics(self.g1, all_quantities)
      spatial.advect(5.0, self.w1)
      spatial.outputCSV(os.path.join(findbin, "2D_advection_out_2_1.csv"), 4)
      with open(os.path.join(findbin, "2D_advection_out_2_1.csv")) as f:
         data = f.readlines()
      self.assertTrue(arrayfuzzyequal([float(d) for d in data[0].strip().split(",")], [14, 0, 6.0, 0], 1E-5))
      self.assertTrue(arrayfuzzyequal([float(d) for d in data[1].strip().split(",")], [0, 0, -16, 0], 1E-5))
      self.assertTrue(arrayfuzzyequal([float(d) for d in data[2].strip().split(",")], [0] * 4, 1E-5))

if __name__ == '__main__':
   unittest.main()

