import os
import sys
import unittest

# so we can find our ../code no matter how we are called
findbin = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(findbin + "/../code")

from grid import Grid

def arrayequal(a, b):
   return all([a[i] == b[i] for i in range(0, len(a))])

def flattenarrayequal(a, b):
   return arrayequal(misc.common.flatten(a), misc.common.flatten(b))

def arraydiff(a, b):
   return sum([x - y for x, y in zip(a, b)])

class TestMe(unittest.TestCase):

   def setUp(self):
      self.g1 = Grid(1.0, 2.0, 3.0, 4, 3)

   def testNegativeCellSize(self):
      with self.assertRaises(ValueError) as the_err:
         Grid(1.0, 2.0, -3.0, 4, 3)
      self.assertEqual(str(the_err.exception), "Cell size must be positive")

   def testGetXmin(self):
      self.assertEqual(self.g1.getXmin(), 1.0)

   def testGetYmin(self):
      self.assertEqual(self.g1.getYmin(), 2.0)

   def testGetNx(self):
      self.assertEqual(self.g1.getNx(), 4)

   def testGetNy(self):
      self.assertEqual(self.g1.getNy(), 3)

   def testGetNumCells(self):
      self.assertEqual(self.g1.getNumCells(), 12)

   def testGetNumActiveCells(self):
      self.assertEqual(self.g1.getNumActiveCells(), 12)

   def testGetNumConnections(self):
      self.assertEqual(self.g1.getNumConnections(), 34)

   def testConnections(self):
      f = self.g1.getConnectionsFrom()
      t = self.g1.getConnectionsTo()
      r = range(len(f))

      self.assertTrue(set([t[i] for i in r if f[i] == 0]) == set([1, 4]))
      self.assertTrue(set([t[i] for i in r if f[i] == 1]) == set([0, 5, 2]))
      self.assertTrue(set([t[i] for i in r if f[i] == 2]) == set([1, 6, 3]))
      self.assertTrue(set([t[i] for i in r if f[i] == 3]) == set([2, 7]))
      self.assertTrue(set([t[i] for i in r if f[i] == 4]) == set([0, 5, 8]))
      self.assertTrue(set([t[i] for i in r if f[i] == 5]) == set([4, 1, 6, 9]))
      self.assertTrue(set([t[i] for i in r if f[i] == 6]) == set([5, 2, 7, 10]))
      self.assertTrue(set([t[i] for i in r if f[i] == 7]) == set([3, 6, 11]))
      self.assertTrue(set([t[i] for i in r if f[i] == 8]) == set([4, 9]))
      self.assertTrue(set([t[i] for i in r if f[i] == 9]) == set([8, 5, 10]))
      self.assertTrue(set([t[i] for i in r if f[i] == 10]) == set([9, 6, 11]))
      self.assertTrue(set([t[i] for i in r if f[i] == 11]) == set([10, 7]))



if __name__ == '__main__':
   unittest.main()

