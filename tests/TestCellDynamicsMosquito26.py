import os
import sys
import unittest
import array
import random

# so we can find our ../code no matter how we are called
findbin = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(findbin + "/../code")

from cellDynamics import CellDynamicsMosquito26

def arrayequal(a, b):
   return all([a[i] == b[i] for i in range(0, len(a))])

def arrayfuzzyequal(a, b, eps):
   return all([(a[i] > b[i] - eps and a[i] < b[i] + eps) for i in range(0, len(a))])

class TestCellDynamicsMosquito26(unittest.TestCase):

   def setUp(self):
      self.c = CellDynamicsMosquito26()

   def testGetNumberOfPopulations(self):
      self.assertEqual(self.c.getNumberOfPopulations(), 2 * 24)
      self.c.setNumAges(4)
      self.assertEqual(self.c.getNumberOfPopulations(), 4 * 24)

   def testGetNumberOfParameters(self):
      self.assertEqual(self.c.getNumberOfParameters(), 2)

   def testGetNumberOfDiffusingPopulations(self):
      self.assertEqual(self.c.getNumberOfDiffusingPopulations(), 24)
      self.c.setNumAges(4)
      self.assertEqual(self.c.getNumberOfDiffusingPopulations(), 24)

   def testGetDiffusingIndices(self):
      self.assertTrue(arrayequal(self.c.getDiffusingIndices(), list(range(24, 48))))
      self.c.setNumAges(4)
      self.assertTrue(arrayequal(self.c.getDiffusingIndices(), list(range(24 * 3, 24 * 4))))

   def testGetNumberOfAdvectingPopulations(self):
      self.assertEqual(self.c.getNumberOfAdvectingPopulations(), 24)
      self.c.setNumAges(4)
      self.assertEqual(self.c.getNumberOfAdvectingPopulations(), 24)

   def testGetAdvectingIndices(self):
      self.assertTrue(arrayequal(self.c.getAdvectingIndices(), list(range(24, 48))))
      self.c.setNumAges(5)
      self.assertTrue(arrayequal(self.c.getAdvectingIndices(), list(range(24 * 4, 24 * 5))))

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

   def testSetGetAlpha(self):
      self.assertEqual(self.c.getAlphaComponentFromPython(0, 0), 1.0)
      self.assertTrue(arrayfuzzyequal([self.c.getAlphaComponentFromPython(0, 1)], [0.4], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getAlphaComponentFromPython(1, 0)], [0.4], 1E-6))
      self.assertEqual(self.c.getAlphaComponentFromPython(1, 1), 1.0)
      with self.assertRaises(ValueError) as the_err:
         self.c.setAlphaComponent(0, 2, 0)
      self.assertEqual(str(the_err.exception), "sp0 0 and sp1 2 must be less than the number of species, 2")
      with self.assertRaises(ValueError) as the_err:
         self.c.setAlphaComponent(2, 1, 0)
      self.assertEqual(str(the_err.exception), "sp0 2 and sp1 1 must be less than the number of species, 2")

   def testSetGetHybridisationRate(self):
      self.assertEqual(self.c.getHybridisationRateFromPython(0, 0, 0), 1.0)
      self.assertEqual(self.c.getHybridisationRateFromPython(1, 1, 1), 1.0)
      self.assertEqual(self.c.getHybridisationRateFromPython(0, 1, 0), 0.5)
      self.assertEqual(self.c.getHybridisationRateFromPython(0, 1, 1), 0.5)
      self.assertEqual(self.c.getHybridisationRateFromPython(1, 0, 0), 0.5)
      self.assertEqual(self.c.getHybridisationRateFromPython(1, 0, 1), 0.5)

   def testSetGetMating(self):
      self.assertEqual(self.c.getMatingComponentFromPython(0, 0), 1.0)
      self.assertEqual(self.c.getMatingComponentFromPython(1, 1), 1.0)
      self.assertTrue(arrayfuzzyequal([self.c.getMatingComponentFromPython(1, 0)], [0.01], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getMatingComponentFromPython(0, 1)], [0.01], 1E-6))

   def testInheritance(self):
      # father, mother, offspring
      # ww, wc, wr, cc, cr, rr
      k_c = 0.995
      k_j = 0.02
      k_ne = 1E-4
      w = 0.5 - 0.5 * k_c
      c = 0.5 + 0.5 * k_c * (1 - k_j) * (1 - k_ne)
      r = 0.5 * k_c * (k_ne + k_j * (1 - k_ne))

      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 0, 0)], [1], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 0, 1)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 0, 2)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 0, 3)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 0, 4)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 0, 5)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 1, 0)], [w], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 1, 1)], [c], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 1, 2)], [r], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 1, 3)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 1, 4)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 1, 5)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 2, 0)], [0.5], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 2, 1)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 2, 2)], [0.5], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 2, 3)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 2, 4)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 2, 5)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 3, 0)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 3, 1)], [1], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 3, 2)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 3, 3)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 3, 4)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 3, 5)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 4, 0)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 4, 1)], [0.5], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 4, 2)], [0.5], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 4, 3)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 4, 4)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 4, 5)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 5, 0)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 5, 1)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 5, 2)], [1], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 5, 3)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 5, 4)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 5, 5)], [0], 1E-6))

      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 1, 0)], [w * w], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 1, 1)], [2 * w * c], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 1, 2)], [2 * w * r], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 1, 3)], [c * c], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 1, 4)], [2 * c * r], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 1, 5)], [r * r], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 2, 0)], [0.5 * w], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 2, 1)], [0.5 * c], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 2, 2)], [0.5 * (w + r)], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 2, 3)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 2, 4)], [0.5 * c], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 2, 5)], [0.5 * r], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 3, 0)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 3, 1)], [w], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 3, 2)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 3, 3)], [c], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 3, 4)], [r], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 3, 5)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 4, 0)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 4, 1)], [0.5 * w], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 4, 2)], [0.5 * w], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 4, 3)], [0.5 * c], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 4, 4)], [0.5 * (c + r)], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 4, 5)], [0.5 * r], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 5, 0)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 5, 1)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 5, 2)], [w], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 5, 3)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 5, 4)], [c], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 5, 5)], [r], 1E-6))

      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 2, 0)], [0.25], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 2, 1)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 2, 2)], [0.5], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 2, 3)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 2, 4)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 2, 5)], [0.25], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 3, 0)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 3, 1)], [0.5], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 3, 2)], [0.5], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 3, 3)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 3, 4)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 3, 5)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 4, 0)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 4, 1)], [0.25], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 4, 2)], [0.25], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 4, 3)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 4, 4)], [0.25], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 4, 5)], [0.25], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 5, 0)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 5, 1)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 5, 2)], [0.5], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 5, 3)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 5, 4)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 5, 5)], [0.5], 1E-6))

      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(3, 3, 0)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(3, 3, 1)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(3, 3, 2)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(3, 3, 3)], [1], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(3, 3, 4)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(3, 3, 5)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(3, 4, 0)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(3, 4, 1)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(3, 4, 2)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(3, 4, 3)], [0.5], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(3, 4, 4)], [0.5], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(3, 4, 5)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(3, 5, 0)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(3, 5, 1)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(3, 5, 2)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(3, 5, 3)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(3, 5, 4)], [1], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(3, 5, 5)], [0], 1E-6))

      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(4, 4, 0)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(4, 4, 1)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(4, 4, 2)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(4, 4, 3)], [0.25], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(4, 4, 4)], [0.5], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(4, 4, 5)], [0.25], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(4, 5, 0)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(4, 5, 1)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(4, 5, 2)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(4, 5, 3)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(4, 5, 4)], [0.5], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(4, 5, 5)], [0.5], 1E-6))

      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(5, 5, 0)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(5, 5, 1)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(5, 5, 2)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(5, 5, 3)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(5, 5, 4)], [0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(5, 5, 5)], [1], 1E-6))

      for male in range(6):
         for female in range(6):
            for offspring in range(6):
               self.assertEqual(self.c.getInheritanceFromPython(male, female, offspring), self.c.getInheritanceFromPython(female, male, offspring))
            
      with self.assertRaises(ValueError) as the_err:
         self.c.getInheritanceFromPython(6, 0, 0)
      self.assertEqual(str(the_err.exception), "All genotypes, 6, 0, 0 must be less than the number of genotypes, 6")
      with self.assertRaises(ValueError) as the_err:
         self.c.getInheritanceFromPython(0, 6, 0)
      self.assertEqual(str(the_err.exception), "All genotypes, 0, 6, 0 must be less than the number of genotypes, 6")
      with self.assertRaises(ValueError) as the_err:
         self.c.getInheritanceFromPython(0, 0, 6)
      self.assertEqual(str(the_err.exception), "All genotypes, 0, 0, 6 must be less than the number of genotypes, 6")



if __name__ == '__main__':
   unittest.main()

