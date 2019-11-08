import os
import sys
import unittest

# so we can find our ../code no matter how we are called
findbin = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(findbin + "/../code")

from grid import Grid
from cellDynamics import CellDynamicsStatic15_9_3_2
from populationsAndParameters import PopulationsAndParameters

def arrayequal(a, b):
   return all([a[i] == b[i] for i in range(0, len(a))])

def arrayfuzzyequal(a, b, eps):
   return all([(a[i] > b[i] - eps and a[i] < b[i] + eps) for i in range(0, len(a))])

class TestDiffusion_1(unittest.TestCase):

   def setUp(self):
      # grid that is 2kmx2km long with 8x8 active cells
      nx = 8
      self.grid = Grid(-0.25 * nx, -0.25 * nx, 0.5, nx, nx)

      self.cell = CellDynamicsStatic15_9_3_2()

      # all the population and parameters
      self.pap = PopulationsAndParameters(self.grid, self.cell)

   def testGetCell(self):
      self.assertEqual(self.pap.getCell().getNumberOfPopulations(), 15)
      self.assertEqual(self.pap.getCell().getNumberOfDiffusingPopulations(), 9)
      self.assertEqual(self.pap.getCell().getNumberOfAdvectingPopulations(), 3)
      self.assertEqual(self.pap.getCell().getNumberOfParameters(), 2)
      self.assertTrue(arrayequal(self.pap.getCell().getDiffusingIndices(), [0, 2, 3, 4, 5, 6, 8, 10, 11]))
      self.assertTrue(arrayequal(self.pap.getCell().getAdvectingIndices(), [0, 4, 10]))

   def testBadSetPopulationAndParameters(self):
      with self.assertRaises(ValueError) as the_err:
         self.pap.setPopulationAndParameters(0, list(range(12345)))
      self.assertEqual(str(the_err.exception), "Length of pop_and_params in setPopulationAndParameters is 12345 which should be 17")
      with self.assertRaises(ValueError) as the_err:
         self.pap.setPopulationAndParameters(1234567, list(range(17)))
      self.assertEqual(str(the_err.exception), "Active cell index is 1234567 which should be less than " + str(8 * 8))

   def testSetPopulationAndParameters(self):
      result = [0.0] * 8 * 8 * 17
      self.pap.setPopulationAndParameters(0, list(range(17)))
      for i in range(17):
         result[i] = 1.0 * i
      self.assertTrue(arrayequal(self.pap.getQuantities()[0:17], result))

      pp = [1.0, -3.0, -33.0, 9.0, 17.0, 1.0, -3.0, -33.0, 9.0, 17.0, 11.0, -19.0, 13.0, -55.0, 15.0, 0.0, 1.0]
      self.pap.setPopulationAndParameters(2, pp)
      for i in range(17):
         result[17 * 2 + i] = pp
      self.assertTrue(arrayequal(self.pap.getQuantities()[0:17], result))


if __name__ == '__main__':
   unittest.main()

