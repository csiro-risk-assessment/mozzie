import os
import sys
import unittest
import array

# so we can find our ../code no matter how we are called
findbin = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(findbin + "/../code")

from cellDynamics import CellDynamicsMosquito23

def arrayequal(a, b):
   return all([a[i] == b[i] for i in range(0, len(a))])

def arrayfuzzyequal(a, b, eps):
   return all([(a[i] > b[i] - eps and a[i] < b[i] + eps) for i in range(0, len(a))])

class TestCellDynamicsMosquito23(unittest.TestCase):

   def setUp(self):
      self.c = CellDynamicsMosquito23()

   def testGetNumberOfPopulations(self):
      self.assertEqual(self.c.getNumberOfPopulations(), 2 * 3 * 2 * 1)
      self.c.setNumAges(4)
      self.c.setNumSpecies(5)
      self.assertEqual(self.c.getNumberOfPopulations(), 2 * 3 * 4 * 5)

   def testGetNumberOfParameters(self):
      self.assertEqual(self.c.getNumberOfParameters(), 1)

   def testGetNumberOfDiffusingPopulations(self):
      self.assertEqual(self.c.getNumberOfDiffusingPopulations(), 6)
      self.c.setNumAges(4)
      self.assertEqual(self.c.getNumberOfDiffusingPopulations(), 6)

   def testGetDiffusingIndices(self):
      self.assertTrue(arrayequal(self.c.getDiffusingIndices(), [6, 7, 8, 9, 10, 11]))
      self.c.setNumAges(4)
      self.assertTrue(arrayequal(self.c.getDiffusingIndices(), [18, 19, 20, 21, 22, 23]))

   def testGetNumberOfAdvectingPopulations(self):
      self.assertEqual(self.c.getNumberOfAdvectingPopulations(), 6)
      self.c.setNumAges(4)
      self.assertEqual(self.c.getNumberOfAdvectingPopulations(), 6)

   def testGetAdvectingIndices(self):
      self.assertTrue(arrayequal(self.c.getAdvectingIndices(), [6, 7, 8, 9, 10, 11]))
      self.c.setNumAges(5)
      self.assertTrue(arrayequal(self.c.getAdvectingIndices(), [24, 25, 26, 27, 28, 29]))

   def testSetGetMuLarvae(self):
      self.c.setMuLarvae(1234.0)
      self.assertEqual(self.c.getMuLarvae(), 1234.0)

   def testSetGetMuAdult(self):
      self.c.setMuAdult(123.0)
      self.assertEqual(self.c.getMuAdult(), 123.0)

   def testSetGetFecundity(self):
      self.c.setFecundity(-123.0)
      self.assertEqual(self.c.getFecundity(), -123.0)

   def testSetGetAgingRate(self):
      self.c.setAgingRate(777.0)
      self.assertEqual(self.c.getAgingRate(), 777.0)

   def testSetGetNumAges(self):
      self.c.setNumAges(7)
      self.assertEqual(self.c.getNumAges(), 7)

   def testSetGetNumSpecies(self):
      self.c.setNumSpecies(4)
      self.assertEqual(self.c.getNumSpecies(), 4)

   def testSetGetAccuracy(self):
      self.c.setAccuracy(0.125)
      self.assertEqual(self.c.getAccuracy(), 0.125)

   def testEvolve(self):
      # THIS TEST CURRENTLY ONLY TESTS THAT evolve CAN BE CALLED WITHOUT RAISING AN ERROR
      dt = 0.1
      pap = array.array('f', list(range(13)))
      self.c.evolve(0.1, pap)

if __name__ == '__main__':
   unittest.main()

