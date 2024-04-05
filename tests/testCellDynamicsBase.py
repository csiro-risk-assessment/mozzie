import os
import sys
import unittest
import array

# so we can find our ../code no matter how we are called
findbin = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(findbin + "/../code")

from cellDynamics import CellDynamicsBase

def arrayequal(a, b):
   return (len(a) == len(b)) and all([a[i] == b[i] for i in range(min(len(a), len(b)))])

class TestCellDynamicsBase(unittest.TestCase):

   def setUp(self):
      self.c = CellDynamicsBase()

   def testGetNumberOfPopulations(self):
      self.assertEqual(self.c.getNumberOfPopulations(), 0)

   def testGetNumberOfDiffusingPopulations(self):
      self.assertEqual(self.c.getNumberOfDiffusingPopulations(), 0)

   def testGetDiffusingIndices(self):
      self.assertEqual(len(self.c.getDiffusingIndices()), 0)

   def testGetNumberOfAdvectingPopulations(self):
      self.assertEqual(self.c.getNumberOfAdvectingPopulations(), 0)

   def testGetAdvectingIndices(self):
      self.assertEqual(len(self.c.getAdvectingIndices()), 0)

   def testGetAdvectionClass(self):
      self.assertEqual(len(self.c.getAdvectionClass()), 0)

   def testSetGetAdvectionClass(self):
      with self.assertRaises(ValueError) as the_err:
         self.c.setAdvectionClass(1, 23)
      self.assertEqual(str(the_err.exception), "setAdvectionClass: population index 1 has not defined to be advecting")

   def testGetNumberOfParameters(self):
      self.assertEqual(self.c.getNumberOfParameters(), 0)

   def testGetSetSmallValue(self):
      self.assertEqual(self.c.getSmallValue(), 0.0)
      self.c.setSmallValue(1.0)
      self.assertEqual(self.c.getSmallValue(), 1.0)

   def testEvolve(self):
      pop_and_params = array.array('f', [1.0])
      self.c.evolve(2.0, pop_and_params)
      self.assertTrue(arrayequal(pop_and_params, [1.0]))

   def testCalcQm(self):
      eqm_pop_and_params = array.array('f', [1.0])
      self.c.calcQm(eqm_pop_and_params)
      self.assertTrue(arrayequal(eqm_pop_and_params, [1.0]))

   def testGetNumSpecies(self):
      self.assertEqual(self.c.getNumSpecies(), 0)

if __name__ == '__main__':
   unittest.main()

