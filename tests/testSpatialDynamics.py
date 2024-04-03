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

   def testBadInit(self):
      bad_grid = Grid(0, 0, 1, 7, 7)
      with self.assertRaises(ValueError) as the_err:
         spatial = SpatialDynamics(bad_grid, self.all_quantities)
      self.assertEqual(str(the_err.exception), "Incorrect size of all_quantities: 1088 which should be the product of 49 and 17")
      
   def testExcept(self):
      fn = os.path.join(findbin, "test_sd_out_1.csv")
      with self.assertRaises(ValueError) as the_err:
         self.spatial.outputCSVsubset(1, 2, 9, 10, fn, 2, "0", "")
      self.assertEqual(str(the_err.exception), "Cell range (1,2) to (9,10) outside range of grid with x_max = 8 and y_max = 8")
      if os.path.isfile(fn) and delete_output: os.remove(fn)

   def testSetGetBirthQuantities(self):
      bqs = list(range(self.cell.getNumberOfDiffusingPopulations() * self.grid.getNumActiveCells()))
      with self.assertRaises(ValueError) as the_err:
         self.spatial.setBirthQuantities(bqs[:-1])
      self.assertEqual(str(the_err.exception), "Incorrect size of birth_quantities, 575, which should be the total number of diffusing populations, 576")

      self.spatial.setBirthQuantities(bqs)
      self.assertEqual(list(self.spatial.getBirthQuantities()), bqs)
      
   def testSetGetYypQuantities(self):
      yyp = list(range(self.cell.getNumberOfDiffusingPopulations() * self.grid.getNumActiveCells()))
      with self.assertRaises(ValueError) as the_err:
         self.spatial.setYypQuantities(yyp[:-1])
      self.assertEqual(str(the_err.exception), "Incorrect size of yyp_quantities, 575, which should be the total number of diffusing populations, 576")

      self.spatial.setYypQuantities(yyp)
      self.assertEqual(list(self.spatial.getYypQuantities()), yyp)

   def testSetGetAllQuantities(self):
      all_q = list(range((self.cell.getNumberOfPopulations() + self.cell.getNumberOfParameters()) * self.grid.getNumActiveCells()))
      with self.assertRaises(ValueError) as the_err:
         self.spatial.setAllQuantities(all_q[:-1])
      self.assertEqual(str(the_err.exception), "Incorrect size of all_quantities, 1087, which should be the total number of quantities, 1088")

      self.spatial.setAllQuantities(all_q)
      self.assertEqual(list(self.spatial.getAllQuantities()), all_q)

   def testBadOutputCSV(self):
      with self.assertRaises(ValueError) as the_err:
         self.spatial.outputCSV("fn.csv", 999, "0", "")
      self.assertEqual(str(the_err.exception), "You requested pop_or_param number 999 but there are only 17 at each cell")
      with self.assertRaises(ValueError) as the_err:
         self.spatial.outputCSV("fn.csv", -999, "0", "")
      self.assertEqual(str(the_err.exception), "You requested pop_or_param number -999 but the minimum accepted value is -18")

   def testOutputCSVsubset(self):
      fn = os.path.join(findbin, "test_sd_out_2.csv")
      self.spatial.outputCSVsubset(3, 3, 5, 6, fn, 2, "0", "")
      with open(fn) as f:
         data = f.readlines()
      self.assertTrue(arrayfuzzyequal([float(d) for d in data[2].strip().split(",")], [0, 0], 1E-8))
      self.assertTrue(arrayfuzzyequal([float(d) for d in data[3].strip().split(",")], [0, 2], 1E-8))
      self.assertTrue(arrayfuzzyequal([float(d) for d in data[4].strip().split(",")], [0, 0], 1E-8))
      if os.path.isfile(fn) and delete_output: os.remove(fn)

   def testOutputBirthQuantities(self):
      fn = os.path.join(findbin, "test_sd_out_3.csv")

      # create some inactive cells
      self.grid.setActiveAndInactive("inactive_active_8x8.csv")
      all_quantities = PopulationsAndParameters(self.grid, self.cell)
      spatial = SpatialDynamics(self.grid, all_quantities)

      # create some birth quantities and write some of them
      diff = self.cell.getNumberOfDiffusingPopulations()
      bqs = list(range(diff * self.grid.getNumActiveCells()))
      spatial.setBirthQuantities(bqs)
      pop_num = 3
      spatial.outputCSVsubset(1, 3, 4, 6, fn, -1 - pop_num, "-999", "")

      # check correct values
      with open(fn) as f:
         data = f.readlines()
      start_val = pop_num + 12 * diff
      self.assertTrue(arrayfuzzyequal([float(d) for d in data[2].strip().split(",")], [start_val, -999, -999], 1E-8))
      start_val += 2 * diff
      self.assertTrue(arrayfuzzyequal([float(d) for d in data[3].strip().split(",")], [start_val, start_val + diff, start_val + 2 * diff], 1E-8))
      start_val += 6 * diff
      self.assertTrue(arrayfuzzyequal([float(d) for d in data[4].strip().split(",")], [-999, -999, start_val], 1E-8))
      if os.path.isfile(fn) and delete_output: os.remove(fn)

   def testOutputYyp(self):
      fn = os.path.join(findbin, "test_sd_out_4.csv")

      # create some inactive cells
      self.grid.setActiveAndInactive("inactive_active_8x8.csv")
      all_quantities = PopulationsAndParameters(self.grid, self.cell)
      spatial = SpatialDynamics(self.grid, all_quantities)

      # create some Yyp and write some of them
      diff = self.cell.getNumberOfDiffusingPopulations()
      yyp = list(range(diff * self.grid.getNumActiveCells()))
      spatial.setYypQuantities(yyp)
      pop_num = 3
      spatial.outputCSVsubset(1, 3, 4, 6, fn, -1 - diff - pop_num, "-999", "")

      # check correct values
      with open(fn) as f:
         data = f.readlines()
      start_val = pop_num + 12 * diff
      self.assertTrue(arrayfuzzyequal([float(d) for d in data[2].strip().split(",")], [start_val, -999, -999], 1E-8))
      start_val += 2 * diff
      self.assertTrue(arrayfuzzyequal([float(d) for d in data[3].strip().split(",")], [start_val, start_val + diff, start_val + 2 * diff], 1E-8))
      start_val += 6 * diff
      self.assertTrue(arrayfuzzyequal([float(d) for d in data[4].strip().split(",")], [-999, -999, start_val], 1E-8))
      if os.path.isfile(fn) and delete_output: os.remove(fn)

if __name__ == '__main__':
   unittest.main()

