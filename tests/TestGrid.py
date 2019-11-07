import os
import sys
import unittest

# so we can find our ../code no matter how we are called
findbin = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(findbin + "/../code")

from grid import Grid

def arrayequal(a, b):
   return all([a[i] == b[i] for i in range(0, len(a))])

class TestGrid(unittest.TestCase):

   def setUp(self):
      self.g1 = Grid(1.0, 2.0, 3.0, 4, 3)
      self.g2 = Grid(1.0, 2.0, 3.0, 4, 3)
      self.g2.setActiveAndInactive(os.path.join(findbin, "inactive_active.csv"))

   def testGetActiveFilename(self):
      self.assertTrue("None", self.g1.getActiveFilename())
      self.assertTrue("inactive_active.csv", self.g2.getActiveFilename())

   def testNegativeCellSize(self):
      with self.assertRaises(ValueError) as the_err:
         Grid(1.0, 2.0, -3.0, 4, 3)
      self.assertEqual(str(the_err.exception), "Cell size must be positive")

   def testGetXmin(self):
      self.assertEqual(self.g1.getXmin(), 1.0)
      self.assertEqual(self.g2.getXmin(), 1.0)

   def testGetYmin(self):
      self.assertEqual(self.g1.getYmin(), 2.0)
      self.assertEqual(self.g2.getYmin(), 2.0)

   def testGetCellSize(self):
      self.assertEqual(self.g1.getCellSize(), 3.0)
      self.assertEqual(self.g2.getCellSize(), 3.0)

   def testGetNx(self):
      self.assertEqual(self.g1.getNx(), 4)
      self.assertEqual(self.g2.getNx(), 4)

   def testGetNy(self):
      self.assertEqual(self.g1.getNy(), 3)
      self.assertEqual(self.g2.getNy(), 3)

   def testGetNumCells(self):
      self.assertEqual(self.g1.getNumCells(), 12)
      self.assertEqual(self.g2.getNumCells(), 12)

   def testGetNumActiveCells(self):
      self.assertEqual(self.g1.getNumActiveCells(), 12)
      self.assertEqual(self.g2.getNumActiveCells(), 7)

   def testGetNumConnections(self):
      self.assertEqual(self.g1.getNumConnections(), 34)
      self.assertEqual(self.g2.getNumConnections(), 10)

   def testGetGlobalIndex(self):
      self.assertTrue(arrayequal(self.g1.getGlobalIndex(), range(12)))
      self.assertTrue(arrayequal(self.g2.getGlobalIndex(), [0, 2, 3, 4, 6, 10, 11]))

   def testGetActiveIndex(self):
      self.assertTrue(arrayequal(self.g1.getActiveIndex(), range(12)))
      self.assertTrue(arrayequal(self.g2.getActiveIndex(), [0, 12, 1, 2, 3, 12, 4, 12, 12, 12, 5, 6]))

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


      f = self.g2.getConnectionsFrom()
      t = self.g2.getConnectionsTo()
      r = range(len(f))

      self.assertTrue(set([t[i] for i in r if f[i] == 0]) == set([3]))
      self.assertTrue(set([t[i] for i in r if f[i] == 1]) == set([2, 4]))
      self.assertTrue(set([t[i] for i in r if f[i] == 2]) == set([1]))
      self.assertTrue(set([t[i] for i in r if f[i] == 3]) == set([0]))
      self.assertTrue(set([t[i] for i in r if f[i] == 4]) == set([1, 5]))
      self.assertTrue(set([t[i] for i in r if f[i] == 5]) == set([4, 6]))
      self.assertTrue(set([t[i] for i in r if f[i] == 6]) == set([5]))


   def testNoActiveFile(self):
      with self.assertRaises(IOError) as the_err:
         self.g1.setActiveAndInactive("no_such_file")
      self.assertEqual(str(the_err.exception), "Cannot open or read no_such_file")

   def testBadActiveFile(self):
      for num in range(1, 13):
         num = str(num)
         with self.assertRaises(ValueError) as the_err:
            self.g1.setActiveAndInactive(os.path.join(findbin, "bad_inactive_active" + num + ".csv"))
         self.assertEqual(str(the_err.exception), "Header lines in " + os.path.join(findbin, "bad_inactive_active" + num + ".csv") + " must include #xmin=1.0,ymin=2.0,cell_size=3.0,nx=4,ny=3")

      for num in range(13, 15):
         num = str(num)
         with self.assertRaises(ValueError) as the_err:
            self.g1.setActiveAndInactive(os.path.join(findbin, "bad_inactive_active" + num + ".csv"))
         self.assertEqual(str(the_err.exception), "There must be 4 entries per line in " + os.path.join(findbin, "bad_inactive_active" + num + ".csv"))

      for num in range(15, 16):
         num = str(num)
         with self.assertRaises(ValueError) as the_err:
            self.g1.setActiveAndInactive(os.path.join(findbin, "bad_inactive_active" + num + ".csv"))
         self.assertEqual(str(the_err.exception), "There must be 3 data lines in " + os.path.join(findbin, "bad_inactive_active" + num + ".csv"))

      for num in range(16, 17):
         num = str(num)
         with self.assertRaises(ValueError) as the_err:
            self.g1.setActiveAndInactive(os.path.join(findbin, "bad_inactive_active" + num + ".csv"))
         self.assertEqual(str(the_err.exception), "The data entries in " + os.path.join(findbin, "bad_inactive_active" + num + ".csv must be either 0 or 1"))

   def testOutputActiveCSV(self):
      fn = os.path.join(findbin, "active_output.csv")
      if os.path.isfile(fn): os.remove(fn)
      self.g2.outputActiveCSV(fn)
      f = open(fn, 'r')
      data = f.readlines()
      f.close()
      self.assertTrue(data[2] == "#xmin=1.0,ymin=2.0,cell_size=3.0,nx=4,ny=3\n")
      self.assertTrue(data[3] == "1,0,1,1\n")
      self.assertTrue(data[4] == "1,0,1,0\n")
      self.assertTrue(data[5] == "0,0,1,1\n")


if __name__ == '__main__':
   unittest.main()

