import os
import sys
import unittest
import array

# so we can find our ../code no matter how we are called
findbin = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(findbin + "/../code")

from cellDynamics import CellDynamicsMosquito26Delay

def arrayequal(a, b):
   return (len(a) == len(b)) and all([a[i] == b[i] for i in range(min(len(a), len(b)))])

def arrayfuzzyequal(a, b, eps):
   return (len(a) == len(b)) and all([(a[i] > b[i] - eps and a[i] < b[i] + eps) for i in range(min(len(a), len(b)))])

class TestCellDynamicsMosquito26Delay(unittest.TestCase):

   def setUp(self):
      self.c = CellDynamicsMosquito26Delay()

   def testSetGetDelay(self):
      self.c.setParameters(22, 0, 4, [1.0] * 4 * 6)
      self.assertEqual(self.c.getDelay(), 22)

   def testGetNumberOfPopulations(self):
      self.c.setParameters(22, 123, 4, [1.0] * 4 * 6)
      self.assertEqual(self.c.getNumberOfPopulations(), 12 * 23 * 4)

   def testGetNumberOfDiffusingPopulations(self):
      self.c.setParameters(11, 123, 4, [1.0] * 4 * 6)
      self.assertEqual(self.c.getNumberOfDiffusingPopulations(), 12 * 4)

   def testGetNumberOfAdvectingPopulations(self):
      self.c.setParameters(11, 123, 4, [1.0] * 4 * 6)
      self.assertEqual(self.c.getNumberOfAdvectingPopulations(), 12 * 4)

   def testGetAdvectingIndices(self):
      self.c.setParameters(11, 7, 4, [1.0] * 4 * 6)
      gold = list(range(12 * 4 * 7, 12 * 4 * (7 + 1)))
      self.assertTrue(arrayequal(self.c.getAdvectingIndices(), gold))

   def testGetDiffusingIndices(self):
      self.c.setParameters(11, 7, 4, [1.0] * 4 * 6)
      gold = list(range(12 * 4 * 7, 12 * 4 * (7 + 1)))
      self.assertTrue(arrayequal(self.c.getDiffusingIndices(), gold))

   def testGetIncrementCurrentIndex(self):
      self.c.setParameters(5, 2, 4, [1.0] * 4 * 6)
      for i in range(13):
         current_index = self.c.getCurrentIndex()
         self.c.incrementCurrentIndex()
         current_index = (current_index + 1) % (self.c.getDelay() + 1)
         self.assertEqual(self.c.getCurrentIndex(), current_index)
         gold = list(range(12 * 4 * current_index, 12 * 4 * (current_index + 1)))
         self.assertTrue(arrayequal(self.c.getDiffusingIndices(), gold))

   def testGetNumberOfParameters(self):
      self.c.setParameters(5, 2, 14, [1.0] * 14 * 6)
      self.assertEqual(self.c.getNumberOfParameters(), 14)

   def testBadSetDeathRate(self):
      with self.assertRaises(ValueError) as the_err:
         self.c.setParameters(5, 2, 14, [3, 1, 2])
      self.assertEqual(str(the_err.exception), "size of death_rate, 3, must be equal to 6 * 14")
      with self.assertRaises(ValueError) as the_err:
         self.c.setParameters(5, 2, 22, list(range(-1, 22 * 6 - 1)))
      self.assertEqual(str(the_err.exception), "all death rates must be positive")

      self.c.setParameters(5, 2, 9, [2.0] * 6 * 9)
      with self.assertRaises(ValueError) as the_err:
         self.c.setDeathRate([2, 3, 4, 5])
      self.assertEqual(str(the_err.exception), "size of death_rate, 4, must be equal to 6 * 9")
      with self.assertRaises(ValueError) as the_err:
         self.c.setDeathRate(list(range(6 * 9)))
      self.assertEqual(str(the_err.exception), "all death rates must be positive")

   def testSetGetDeathRate(self):
      self.c.setParameters(5, 2, 2, list(range(1, 13)))
      self.assertTrue(arrayequal(self.c.getDeathRate(), list(range(1, 13))))

      
if __name__ == '__main__':
   unittest.main()

