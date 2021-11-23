import os
import sys
import unittest
import array
import random

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
      self.c.setDelayCurrentIndexNumSpecies(22, 0, 4)
      self.assertEqual(self.c.getDelay(), 22)

   def testGetIncrementCurrentIndex(self):
      current_index = self.c.getCurrentIndex()
      self.c.incrementCurrentIndex()
      current_index = (current_index + 1) % (self.c.getDelay() + 1)
      self.assertEqual(self.c.getCurrentIndex(), current_index)

if __name__ == '__main__':
   unittest.main()

