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
      self.c.setNumSpecies(3)
      self.assertEqual(self.c.getNumberOfDiffusingPopulations(), 6 * 3)

   def testGetDiffusingIndices(self):
      self.assertTrue(arrayequal(self.c.getDiffusingIndices(), [6, 7, 8, 9, 10, 11]))
      self.c.setNumAges(4)
      self.c.setNumSpecies(2)
      self.assertTrue(arrayequal(self.c.getDiffusingIndices(), list(range(36, 48))))

   def testGetNumberOfAdvectingPopulations(self):
      self.assertEqual(self.c.getNumberOfAdvectingPopulations(), 6)
      self.c.setNumAges(4)
      self.c.setNumSpecies(5)
      self.assertEqual(self.c.getNumberOfAdvectingPopulations(), 6 * 5)

   def testGetAdvectingIndices(self):
      self.assertTrue(arrayequal(self.c.getAdvectingIndices(), [6, 7, 8, 9, 10, 11]))
      self.c.setNumAges(5)
      self.c.setNumSpecies(2)
      self.assertTrue(arrayequal(self.c.getAdvectingIndices(), list(range(48, 60))))

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

   def testEvolveZeroFecundityZeroAging(self):
      dt = 0.01
      self.c.setFecundity(0.0)
      self.c.setAgingRate(0.0)
      self.c.setMuLarvae(0.5)
      self.c.setMuAdult(0.7)

      initial_condition = list(range(13))
      pap = array.array('f', initial_condition)
      expected_answer = [x * (1 - 0.5 * dt) for x in initial_condition[:6]] + [x * (1 - 0.7 * dt) for x in initial_condition[6:12]] + [12]
      self.c.evolve(dt, pap)
      self.assertTrue(arrayfuzzyequal(pap, expected_answer, 3E-4))

      self.c.setNumSpecies(2)
      self.c.setNumAges(5)
      initial_condition = list(range(self.c.getNumberOfPopulations() + self.c.getNumberOfParameters()))
      pap = array.array('f', initial_condition)
      ### THIS CURRENTLY CRASHES DUE TO INDEXING ERROR: self.c.evolve(dt, pap)

   def testEvolveAgingOnly(self):
      dt = 0.01
      aging_rate = 0.25
      self.c.setFecundity(0.0)
      self.c.setAgingRate(aging_rate)
      self.c.setMuLarvae(0.0)
      self.c.setMuAdult(0.0)
      initial_condition = list(range(13))
      pap = array.array('f', initial_condition)
      expected_answer = [x * (1 - aging_rate * dt) for x in initial_condition[:6]] + [initial_condition[i] + initial_condition[i - 6] * aging_rate * dt for i in range(6, 12)] + [12]
      self.c.evolve(dt, pap)
      self.assertTrue(arrayfuzzyequal(pap, expected_answer, 3E-5))
      conserved = [pap[i] + pap[i + 6] - initial_condition[i] - initial_condition[i + 6] for i in range(6)]
      self.assertTrue(arrayfuzzyequal(conserved, [0.0] * 6, 6E-7))

if __name__ == '__main__':
   unittest.main()

