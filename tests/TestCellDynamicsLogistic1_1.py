import os
import sys
import unittest
import array

# so we can find our ../code no matter how we are called
findbin = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(findbin + "/../code")

from cellDynamics import CellDynamicsLogistic1_1

def arrayequal(a, b):
   return (len(a) == len(b)) and all([a[i] == b[i] for i in range(min(len(a), len(b)))])

def arrayfuzzyequal(a, b, eps):
   return (len(a) == len(b)) and all([(a[i] > b[i] - eps and a[i] < b[i] + eps) for i in range(min(len(a), len(b)))])

class TestCellDynamicsLogistic1_1(unittest.TestCase):

   def setUp(self):
      self.c = CellDynamicsLogistic1_1()

   def testGetNumberOfPopulations(self):
      self.assertEqual(self.c.getNumberOfPopulations(), 1)

   def testGetNumberOfDiffusingPopulations(self):
      self.assertEqual(self.c.getNumberOfDiffusingPopulations(), 1)

   def testGetDiffusingIndices(self):
      self.assertTrue(arrayequal(self.c.getDiffusingIndices(), [0]))

   def testGetNumberOfAdvectingPopulations(self):
      self.assertEqual(self.c.getNumberOfAdvectingPopulations(), 1)

   def testGetAdvectingIndices(self):
      self.assertTrue(arrayequal(self.c.getAdvectingIndices(), [0]))

   def testSetGetAdvectionClass(self):
      self.assertTrue(arrayequal(self.c.getAdvectionClass(), [0]))
      with self.assertRaises(ValueError) as the_err:
         self.c.setAdvectionClass(1, 23)
      self.assertEqual(str(the_err.exception), "setAdvectionClass: population index 1 has not defined to be advecting")
      self.c.setAdvectionClass(0, 23)
      self.assertTrue(arrayequal(self.c.getAdvectionClass(), [23]))

   def testGetNumberOfParameters(self):
      self.assertEqual(self.c.getNumberOfParameters(), 1)

   def testEvolve(self):
      pop_and_params = array.array('f', [1.0, 2.0])
      self.c.evolve(40.0, pop_and_params)
      self.assertTrue(arrayfuzzyequal(pop_and_params, [1.0 + 40.0 * 0.01 * (1.0 - 1.0 / 2.0), 2.0], 1E-5))
      

if __name__ == '__main__':
   unittest.main()

