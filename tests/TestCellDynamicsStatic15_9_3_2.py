import os
import sys
import unittest

# so we can find our ../code no matter how we are called
findbin = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(findbin + "/../code")

from cellDynamics import CellDynamicsStatic15_9_3_2

def arrayequal(a, b):
   return (len(a) == len(b)) and all([a[i] == b[i] for i in range(min(len(a), len(b)))])

class TestCellDynamicsStatic15_9_3_2(unittest.TestCase):

   def setUp(self):
      self.c = CellDynamicsStatic15_9_3_2()

   def testGetNumberOfPopulations(self):
      self.assertEqual(self.c.getNumberOfPopulations(), 15)

   def testGetNumberOfDiffusingPopulations(self):
      self.assertEqual(self.c.getNumberOfDiffusingPopulations(), 9)

   def testGetDiffusingIndices(self):
      self.assertTrue(arrayequal(self.c.getDiffusingIndices(), [0, 2, 3, 4, 5, 6, 8, 10, 11]))

   def testGetNumberOfAdvectingPopulations(self):
      self.assertEqual(self.c.getNumberOfAdvectingPopulations(), 3)

   def testGetAdvectingIndices(self):
      self.assertTrue(arrayequal(self.c.getAdvectingIndices(), [0, 4, 10]))

   def testGetNumberOfParameters(self):
      self.assertEqual(self.c.getNumberOfParameters(), 2)

if __name__ == '__main__':
   unittest.main()

