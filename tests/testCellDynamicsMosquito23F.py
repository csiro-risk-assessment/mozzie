import os
import sys
import unittest
import array
import random

# so we can find our ../code no matter how we are called
findbin = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(findbin + "/../code")

from cellDynamics import CellDynamicsMosquito23F

def arrayequal(a, b):
   return (len(a) == len(b)) and all([a[i] == b[i] for i in range(min(len(a), len(b)))])

def arrayfuzzyequal(a, b, eps):
   return (len(a) == len(b)) and all([(a[i] > b[i] - eps and a[i] < b[i] + eps) for i in range(min(len(a), len(b)))])

class TestCellDynamicsMosquito23F(unittest.TestCase):

   def setUp(self):
      self.c = CellDynamicsMosquito23F()

   def testGetFecundityProportion(self):
      # expected_fecundity[genotype_male][genotype_female][sex_offspring]
      expected_fecundity = [[[0.5, 0.5], [0.45, 0.55], [0.45, 0.55]], # genotype_male=ww
                            [[0.95, 0.05], [0.95, 0.05], [0.95, 0.05]], # genotype_male=Gw
                            [[0.95, 0.05], [0.95, 0.05], [0.95, 0.05]]] # genotype_male=GG
      for gM in range(3):
         for gF in range(3):
            for offspring in range(2):
               self.assertTrue(arrayfuzzyequal([self.c.getFecundityProportion(offspring, gF, gM)], [expected_fecundity[gM][gF][offspring]], 1E-6))

if __name__ == '__main__':
   unittest.main()

