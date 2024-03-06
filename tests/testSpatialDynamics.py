import os
import sys
import unittest

# so we can find our ../code no matter how we are called
findbin = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(findbin + "/../code")

from grid import Grid
from cellDynamics import CellDynamicsStatic15_9_3_2
from spatialDynamics import SpatialDynamics
from populationsAndParameters import PopulationsAndParameters

delete_output = True

def arrayfuzzyequal(a, b, eps):
   return (len(a) == len(b)) and all([(a[i] > b[i] - eps and a[i] < b[i] + eps) for i in range(min(len(a), len(b)))])

class TestSpatialDynamics(unittest.TestCase):

   def setUp(self):
      # grid that is 2kmx2km long with 8x8 active cells
      nx = 8
      self.grid = Grid(-0.25 * nx, -0.25 * nx, 0.5, nx, nx)

      self.cell = CellDynamicsStatic15_9_3_2()

      # all the population and parameters
      self.all_quantities = PopulationsAndParameters(self.grid, self.cell)
      # centre cell starts with nonzero population
      pop = list(range(self.cell.getNumberOfPopulations() + self.cell.getNumberOfParameters()))
      self.all_quantities.setPopulationAndParameters((nx * (nx + 1)) // 2, pop)

      self.spatial = SpatialDynamics(self.grid, self.all_quantities)


   def testExcept(self):
      fn = os.path.join(findbin, "test_sd_out_1.csv")
      with self.assertRaises(ValueError) as the_err:
         self.spatial.outputCSVsubset(1, 2, 9, 10, fn, 2, "0", "")
      self.assertEqual(str(the_err.exception), "Cell range (1,2) to (9,10) outside range of grid with x_max = 8 and y_max = 8")
      if os.path.isfile(fn) and delete_output: os.remove(fn)

   def testOutputCSVsubset(self):
      fn = os.path.join(findbin, "test_sd_out_2.csv")
      self.spatial.outputCSVsubset(3, 3, 5, 6, fn, 2, "0", "")
      with open(fn) as f:
         data = f.readlines()
      self.assertTrue(arrayfuzzyequal([float(d) for d in data[2].strip().split(",")], [0, 0], 1E-8))
      self.assertTrue(arrayfuzzyequal([float(d) for d in data[3].strip().split(",")], [0, 2], 1E-8))
      self.assertTrue(arrayfuzzyequal([float(d) for d in data[4].strip().split(",")], [0, 0], 1E-8))
      if os.path.isfile(fn) and delete_output: os.remove(fn)

if __name__ == '__main__':
   unittest.main()

