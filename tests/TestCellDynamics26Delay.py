import os
import sys
import unittest
import array
import random
from math import exp

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
      self.c.setParameters(22, 0, 4, [1.0] * 4 * 6, [0.0] * 4 * 4, [0.0] * 2 * 6 * 4, [0.0] * 4 * 4)
      self.assertEqual(self.c.getDelay(), 22)

   def testGetNumberOfPopulations(self):
      self.c.setParameters(22, 123, 4, [1.0] * 4 * 6, [0.0] * 4 * 4, [0.0] * 2 * 6 * 4, [0.0] * 4 * 4)
      self.assertEqual(self.c.getNumberOfPopulations(), 12 * 23 * 4)

   def testGetNumberOfDiffusingPopulations(self):
      self.c.setParameters(11, 123, 4, [1.0] * 4 * 6, [0.0] * 4 * 4, [0.0] * 2 * 6 * 4, [0.0] * 4 * 4)
      self.assertEqual(self.c.getNumberOfDiffusingPopulations(), 12 * 4)

   def testGetNumberOfAdvectingPopulations(self):
      self.c.setParameters(11, 123, 4, [1.0] * 4 * 6, [0.0] * 4 * 4, [0.0] * 2 * 6 * 4, [0.0] * 4 * 4)
      self.assertEqual(self.c.getNumberOfAdvectingPopulations(), 12 * 4)

   def testGetAdvectingIndices(self):
      self.c.setParameters(11, 7, 4, [1.0] * 4 * 6, [0.0] * 4 * 4, [0.0] * 2 * 6 * 4, [0.0] * 4 * 4)
      gold = list(range(12 * 4 * 7, 12 * 4 * (7 + 1)))
      self.assertTrue(arrayequal(self.c.getAdvectingIndices(), gold))

   def testGetDiffusingIndices(self):
      self.c.setParameters(11, 7, 4, [1.0] * 4 * 6, [0.0] * 4 * 4, [0.0] * 2 * 6 * 4, [0.0] * 4 * 4)
      gold = list(range(12 * 4 * 7, 12 * 4 * (7 + 1)))
      self.assertTrue(arrayequal(self.c.getDiffusingIndices(), gold))

   def testGetIncrementCurrentIndex(self):
      self.c.setParameters(5, 2, 4, [1.0] * 4 * 6, [0.0] * 4 * 4, [0.0] * 2 * 6 * 4, [0.0] * 4 * 4)
      for i in range(13):
         current_index = self.c.getCurrentIndex()
         self.c.incrementCurrentIndex()
         current_index = (current_index + 1) % (self.c.getDelay() + 1)
         self.assertEqual(self.c.getCurrentIndex(), current_index)
         gold = list(range(12 * 4 * current_index, 12 * 4 * (current_index + 1)))
         self.assertTrue(arrayequal(self.c.getDiffusingIndices(), gold))

   def testGetNumberOfParameters(self):
      self.c.setParameters(5, 2, 14, [1.0] * 14 * 6, [0.0] * 14 * 14, [0.0] * 2 * 6 * 14, [0.0] * 14 * 14)
      self.assertEqual(self.c.getNumberOfParameters(), 14)

   def testBadSetDeathRate(self):
      with self.assertRaises(ValueError) as the_err:
         self.c.setParameters(5, 2, 14, [3, 1, 2], [0.0] * 14 * 14, [0.0] * 2 * 6 * 14, [0.0] * 14 * 14)
      self.assertEqual(str(the_err.exception), "size of death_rate, 3, must be equal to 6 * 14")
      with self.assertRaises(ValueError) as the_err:
         self.c.setParameters(5, 2, 22, list(range(-1, 22 * 6 - 1)), [0.0] * 22 * 22, [0.0] * 2 * 6 * 22, [0.0] * 22 * 22)
      self.assertEqual(str(the_err.exception), "all death rates must be positive")

      self.c.setParameters(5, 2, 9, [2.0] * 6 * 9, [0.0] * 9 * 9, [0.0] * 2 * 6 * 9, [0.0] * 9 * 9)
      with self.assertRaises(ValueError) as the_err:
         self.c.setDeathRate([2, 3, 4, 5])
      self.assertEqual(str(the_err.exception), "size of death_rate, 4, must be equal to 6 * 9")
      with self.assertRaises(ValueError) as the_err:
         self.c.setDeathRate(list(range(6 * 9)))
      self.assertEqual(str(the_err.exception), "all death rates must be positive")

   def testSetGetDeathRate(self):
      self.c.setParameters(5, 2, 2, list(range(1, 13)), [0.0] * 2 * 2, [0.0] * 2 * 6 * 2, [0.0] * 2 * 2)
      self.assertTrue(arrayequal(self.c.getDeathRate(), list(range(1, 13))))

   def testEvolveTrial(self):
      # Tests dx/dt = -d * x + lambdah * x(t - delay)
      delay = 5
      current_index = 2
      death_rate = list(range(1, 7))
      self.c.setParameters(delay, current_index, 1, death_rate, [0.0], [0.0] * 2 * 6 * 1, [0.0] * 1 * 1)

      initial_condition = [random.random() for i in range(self.c.getNumberOfPopulations() + self.c.getNumberOfParameters())]
      pap = array.array('f', initial_condition)
      dt = 1.23E-2
      self.c.evolveTrial(dt, pap)
      self.c.incrementCurrentIndex() # not necessary here: just good practice to increment after evolve has been called for all grid cells

      lambdah = 1.1
      expected_answer = list(initial_condition)
      for i in range(24, 36):
         dr = death_rate[(i - 24) % 6]
         new_pop = lambdah * initial_condition[i + 12] / dr + (initial_condition[i] - lambdah * initial_condition[i + 12] / dr) * exp(- dr * dt)
         expected_answer[i + 12] = new_pop

      self.assertTrue(arrayfuzzyequal(pap, expected_answer, 1E-6))
      
   def testBadSetCompetition(self):
      with self.assertRaises(ValueError) as the_err:
         self.c.setParameters(5, 2, 2, list(range(1, 13)), list(range(3)), [0.0] * 2 * 6 * 2, [0.0] * 2 * 2)
      self.assertEqual(str(the_err.exception), "size of competition, 3, must be equal to 2 * 2")

      self.c.setParameters(5, 2, 9, [2.0] * 6 * 9, [0.0] * 9 * 9, [0.0] * 2 * 6 * 9, [0.0] * 9 * 9)
      with self.assertRaises(ValueError) as the_err:
         self.c.setCompetition([2, 3, 4, 5])
      self.assertEqual(str(the_err.exception), "size of competition, 4, must be equal to 9 * 9")

   def testSetGetCompetition(self):
      self.c.setParameters(5, 2, 2, list(range(1, 13)), list(range(4)), [0.0] * 2 * 6 * 2, [0.0] * 2 * 2)
      self.assertTrue(arrayequal(self.c.getCompetition(), list(range(4))))

   def testBadSetEmergenceRate(self):
      with self.assertRaises(ValueError) as the_err:
         self.c.setParameters(5, 2, 14, [3] * 6 * 14, [0.0] * 14 * 14, [0.0] * 4, [0.0] * 14 * 14)
      self.assertEqual(str(the_err.exception), "size of emergence_rate, 4, must be equal to 2 * 6 * 14")
      with self.assertRaises(ValueError) as the_err:
         self.c.setParameters(5, 2, 22, [3] * 6 * 22, [0.0] * 22 * 22, [-1.0] * 2 * 6 * 22, [0.0] * 22 * 22)
      self.assertEqual(str(the_err.exception), "all emergence rates must be non-negative")

      self.c.setParameters(5, 2, 9, [2.0] * 6 * 9, [0.0] * 9 * 9, [0.0] * 2 * 6 * 9, [0.0] * 9 * 9)
      with self.assertRaises(ValueError) as the_err:
         self.c.setEmergenceRate([2, 3, 4, 5])
      self.assertEqual(str(the_err.exception), "size of emergence_rate, 4, must be equal to 2 * 6 * 9")
      with self.assertRaises(ValueError) as the_err:
         self.c.setEmergenceRate(list(range(-1, 2 * 6 * 9 - 1)))
      self.assertEqual(str(the_err.exception), "all emergence rates must be non-negative")

   def testSetGetEmergenceRate(self):
      self.c.setParameters(5, 2, 2, list(range(1, 13)), [0.0] * 2 * 2, list(range(2 * 6 * 2)), [0.0] * 2 * 2)
      self.assertTrue(arrayequal(self.c.getEmergenceRate(), list(range(2 * 6 * 2))))

   def testBadSetActivity(self):
      with self.assertRaises(ValueError) as the_err:
         self.c.setParameters(5, 2, 14, [3] * 6 * 14, [0.0] * 14 * 14, [0.0] * 2 * 6 * 14, [0.0] * 4)
      self.assertEqual(str(the_err.exception), "size of activity, 4, must be equal to 14 * 14")
      with self.assertRaises(ValueError) as the_err:
         self.c.setParameters(5, 2, 22, [3] * 6 * 22, [0.0] * 22 * 22, [0.0] * 2 * 6 * 22, [-1] * 22 * 22)
      self.assertEqual(str(the_err.exception), "all activity values must be non-negative")

      self.c.setParameters(5, 2, 9, [2.0] * 6 * 9, [0.0] * 9 * 9, [0.0] * 2 * 6 * 9, [0.0] * 9 * 9)
      with self.assertRaises(ValueError) as the_err:
         self.c.setActivity([2, 3, 4, 5])
      self.assertEqual(str(the_err.exception), "size of activity, 4, must be equal to 9 * 9")
      with self.assertRaises(ValueError) as the_err:
         self.c.setActivity(list(range(-1, 9 * 9 - 1)))
      self.assertEqual(str(the_err.exception), "all activity values must be non-negative")

   def testSetGetActivity(self):
      self.c.setParameters(5, 2, 2, list(range(1, 13)), [0.0] * 2 * 2, list(range(2 * 6 * 2)), list(range(2 * 2)))
      self.assertTrue(arrayequal(self.c.getActivity(), list(range(2 * 2))))

   def testSetGetSexRatio(self):
      self.c.setFecundityP(0.625, 0.75)
      self.assertEqual(self.c.getSexRatio(), 0.625)

   def testSetGetFemaleBias(self):
      self.c.setFecundityP(0.625, 0.75)
      self.assertEqual(self.c.getFemaleBias(), 0.75)


if __name__ == '__main__':
   unittest.main()

