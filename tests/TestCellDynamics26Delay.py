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
      self.ide4 = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]
      self.red = [list(range(1, 7)), list(range(2, 8)), [1] * 6, list(range(4, 10)), [2] * 6, [1, 2, 1, 2, 1, 2]]
      self.hyb = [[[1, 2, 3, 4], [5, 4, 3, 2], [3, 3, 2, 1], [5, 6, 7, 4]], [[6, 7, 3, 2], [4, 6, 1, 8], [4, 7, 3, 2], [9, 8, 0, 7]], [[5, 7, 2, 9], [4, 5, 3, 6], [4, 3, 1, 2], [4, 2, 6, 8]], [[6, 5, 7, 4], [5, 6, 3, 7], [5, 9, 0, 1], [0, 8, 0, 4]]]
      self.d = CellDynamicsMosquito26Delay(num_species = 4, delay = 17, current_index = 7, death_rate = [[1.0] * 4] * 6, competition = self.ide4, emergence_rate = [1.0] * 4, activity = self.ide4, reduction = self.red, hybridisation = self.hyb, min_cc = 0.125)

   def testSetGetDelay(self):
      self.assertEqual(self.d.getDelay(), 17)

   def testGetNumberOfPopulations(self):
      self.assertEqual(self.d.getNumberOfPopulations(), 12 * 18 * 4)

   def testGetNumberOfDiffusingPopulations(self):
      self.assertEqual(self.d.getNumberOfDiffusingPopulations(), 12 * 4)

   def testGetNumberOfAdvectingPopulations(self):
      self.assertEqual(self.d.getNumberOfAdvectingPopulations(), 12 * 4)

   def testGetAdvectingIndices(self):
      gold = list(range(12 * 4 * 7, 12 * 4 * (7 + 1)))
      self.assertTrue(arrayequal(self.d.getAdvectingIndices(), gold))

   def testGetDiffusingIndices(self):
      gold = list(range(12 * 4 * 7, 12 * 4 * (7 + 1)))
      self.assertTrue(arrayequal(self.d.getDiffusingIndices(), gold))

   def testGetIncrementCurrentIndex(self):
      for i in range(37):
         current_index = self.d.getCurrentIndex()
         self.d.incrementCurrentIndex()
         current_index = (current_index + 1) % (self.d.getDelay() + 1)
         self.assertEqual(self.d.getCurrentIndex(), current_index)
         gold = list(range(12 * 4 * current_index, 12 * 4 * (current_index + 1)))
         self.assertTrue(arrayequal(self.d.getDiffusingIndices(), gold))

   def testGetNumberOfParameters(self):
      self.assertEqual(self.d.getNumberOfParameters(), 4)

   def testBadSetDeathRate(self):
      with self.assertRaises(ValueError) as the_err:
         e = CellDynamicsMosquito26Delay(num_species = 4, delay = 17, current_index = 7, death_rate = [1, 2, 3], competition = self.ide4, emergence_rate = [1.0] * 4, activity = self.ide4)
      self.assertEqual(str(the_err.exception), "size of death_rate, 3, must be equal to 6")
      with self.assertRaises(ValueError) as the_err:
         e = CellDynamicsMosquito26Delay(num_species = 4, delay = 17, current_index = 7, death_rate = [list(range(5))] * 6, emergence_rate = [1.0] * 4, activity = self.ide4)
      self.assertEqual(str(the_err.exception), "size of death_rate[0], 5, must be equal to 4")
      with self.assertRaises(ValueError) as the_err:
         e = CellDynamicsMosquito26Delay(num_species = 4, delay = 17, current_index = 7, death_rate = [list(range(4))] * 6, emergence_rate = [1.0] * 4, activity = self.ide4)
      self.assertEqual(str(the_err.exception), "all death rates must be positive")

      with self.assertRaises(ValueError) as the_err:
         self.d.setDeathRate([2, 3, 4, 5])
      self.assertEqual(str(the_err.exception), "size of death_rate, 4, must be equal to 6")
      with self.assertRaises(ValueError) as the_err:
         self.d.setDeathRate([list(range(-1, 3))] * 6)
      self.assertEqual(str(the_err.exception), "all death rates must be positive")

   def testSetGetDeathRate(self):
      self.assertTrue(arrayequal(self.d.getDeathRate(), [[1.0] * 4] * 6))
      dr = [[i, i + 6, i + 123] for i in range(1, 7)]
      self.c.setDeathRate(dr)
      self.assertTrue(arrayequal(self.c.getDeathRate(), dr))

   def testEvolveTrial(self):
      # Tests dx/dt = -d * x + lambdah * x(t - delay)
      delay = 5
      current_index = 2
      death_rate = [[i] for i in range(1, 7)]
      tiny = CellDynamicsMosquito26Delay(num_species = 1, delay = delay, current_index = current_index, death_rate = death_rate, competition = [[0.0]], emergence_rate = [0.0], activity = [[0.0]], hybridisation = [[[1.0]]])

      initial_condition = [random.random() for i in range(tiny.getNumberOfPopulations() + tiny.getNumberOfParameters())]
      pap = array.array('f', initial_condition)
      dt = 1.23E-2
      tiny.evolveTrial(dt, pap)
      tiny.incrementCurrentIndex() # not necessary here: just good practice to increment after evolve has been called for all grid cells

      lambdah = 1.1
      expected_answer = list(initial_condition)
      for i in range(24, 36):
         dr = death_rate[(i - 24) % 6][0]
         new_pop = lambdah * initial_condition[i + 12] / dr + (initial_condition[i] - lambdah * initial_condition[i + 12] / dr) * exp(- dr * dt)
         expected_answer[i + 12] = new_pop

      self.assertTrue(arrayfuzzyequal(pap, expected_answer, 1E-6))
      
   def testBadSetCompetition(self):
      with self.assertRaises(ValueError) as the_err:
         e = CellDynamicsMosquito26Delay(num_species = 2, delay = 17, current_index = 7, death_rate = [[1.0] * 2] * 6, competition = list(range(3)), emergence_rate = [1.0] * 2, activity = [[0, 0], [0, 0]])
      self.assertEqual(str(the_err.exception), "size of competition, 3, must be equal to 2")
      with self.assertRaises(ValueError) as the_err:
         e = CellDynamicsMosquito26Delay(num_species = 2, delay = 17, current_index = 7, death_rate = [[1.0] * 2] * 6, competition = [[1, 2], [3]], emergence_rate = [1.0] * 2, activity = [[0, 0], [0, 0]])
      self.assertEqual(str(the_err.exception), "size of competition[1], 1, must be equal to 2")

      with self.assertRaises(ValueError) as the_err:
         self.c.setCompetition([2, 3, 4, 5])
      self.assertEqual(str(the_err.exception), "size of competition, 4, must be equal to 3")
      with self.assertRaises(ValueError) as the_err:
         self.c.setCompetition([[0.0] * 7] * 3)
      self.assertEqual(str(the_err.exception), "size of competition[0], 7, must be equal to 3")

   def testSetGetCompetition(self):
      self.assertTrue(arrayequal(self.d.getCompetition(), self.ide4))
      a = [[1, 2, 3], [4, 5, 6], [6, 7, 8]]
      self.c.setCompetition(a)
      self.assertTrue(arrayequal(self.c.getCompetition(), a))

   def testBadSetEmergenceRate(self):
      with self.assertRaises(ValueError) as the_err:
         e = CellDynamicsMosquito26Delay(num_species = 14, delay = 17, current_index = 7, death_rate = [[1.0] * 14] * 6, competition = [[0] * 14] * 14, emergence_rate = [1.0] * 4, activity = [[0] * 14] * 14)
      self.assertEqual(str(the_err.exception), "size of emergence_rate, 4, must be equal to 14")
      with self.assertRaises(ValueError) as the_err:
         e = CellDynamicsMosquito26Delay(num_species = 14, delay = 17, current_index = 7, death_rate = [[1.0] * 14] * 6, competition = [[0] * 14] * 14, emergence_rate = [-1.0] * 14, activity = [[0] * 14] * 14)
      self.assertEqual(str(the_err.exception), "all emergence rates must be non-negative")

      with self.assertRaises(ValueError) as the_err:
         self.d.setEmergenceRate([2, 3, 4])
      self.assertEqual(str(the_err.exception), "size of emergence_rate, 3, must be equal to 4")
      with self.assertRaises(ValueError) as the_err:
         self.d.setEmergenceRate([2, 3, -4, 5])
      self.assertEqual(str(the_err.exception), "all emergence rates must be non-negative")

   def testSetGetEmergenceRate(self):
      self.assertTrue(arrayequal(self.d.getEmergenceRate(), [1.0] * 4))
      self.d.setEmergenceRate([1, 2, 3, 4])
      self.assertTrue(arrayequal(self.d.getEmergenceRate(), [1, 2, 3, 4]))

   def testBadSetActivity(self):
      with self.assertRaises(ValueError) as the_err:
         e = CellDynamicsMosquito26Delay(num_species = 14, delay = 17, current_index = 7, death_rate = [[1.0] * 14] * 6, competition = [[0] * 14] * 14, emergence_rate = [1.0] * 14, activity = [[0] * 14] * 4)
      self.assertEqual(str(the_err.exception), "size of activity, 4, must be equal to 14")
      with self.assertRaises(ValueError) as the_err:
         e = CellDynamicsMosquito26Delay(num_species = 2, delay = 17, current_index = 7, death_rate = [[1.0] * 2] * 6, competition = [[1, 2], [3, 4]], emergence_rate = [1.0] * 2, activity = [[1, 2], [3, 2, 1]])
      self.assertEqual(str(the_err.exception), "size of activity[1], 3, must be equal to 2")
      with self.assertRaises(ValueError) as the_err:
         e = CellDynamicsMosquito26Delay(num_species = 2, delay = 17, current_index = 7, death_rate = [[1.0] * 2] * 6, competition = [[1, 2], [3, 4]], emergence_rate = [1.0] * 2, activity = [[1, 2], [3, -2]])
      self.assertEqual(str(the_err.exception), "all activity values must be non-negative")

      with self.assertRaises(ValueError) as the_err:
         self.c.setActivity([2, 3, 4, 5])
      self.assertEqual(str(the_err.exception), "size of activity, 4, must be equal to 3")
      with self.assertRaises(ValueError) as the_err:
         self.c.setActivity([[1, 2, 3], [3, 4, 5], [5, 6]])
      self.assertEqual(str(the_err.exception), "size of activity[2], 2, must be equal to 3")
      with self.assertRaises(ValueError) as the_err:
         self.c.setActivity([[1, 2, 3], [3, -4, 5], [5, 6, 7]])
      self.assertEqual(str(the_err.exception), "all activity values must be non-negative")

   def testSetGetActivity(self):
      self.assertTrue(arrayequal(self.d.getActivity(), self.ide4))
      a = [[3, 4, 5], [1, 5, 6], [8, 8, 1]]
      self.c.setActivity(a)
      self.assertTrue(arrayequal(self.c.getActivity(), a))

   def testSetGetSexRatio(self):
      self.c.setFecundityP(0.625, 0.75)
      self.assertEqual(self.c.getSexRatio(), 0.625)

   def testSetGetFemaleBias(self):
      self.c.setFecundityP(0.625, 0.75)
      self.assertEqual(self.c.getFemaleBias(), 0.75)

   def testSetGetFecundityP(self):
      ww = 0
      wc = 1
      wr = 2
      cc = 3
      cr = 4
      rr = 5
      male = 0
      female = 1
      sex_ratio = 0.625
      female_bias = 0.75
      p = [[[0.5 for s in range(2)] for gF in range(6)] for gM in range(6)]
      for gM in range(6):
         for gF in range(6):
            for s in range(2):
               if gM == wc or gM == cc:
                  if s == male: p[gM][gF][s] = sex_ratio
                  else: p[gM][gF][s] = 1 - sex_ratio
               elif gM == ww and (gF == wc or gF == cc):
                  if s == male: p[gM][gF][s] = 1 - female_bias
                  else: p[gM][gF][s] = female_bias

      self.c.setFecundityP(sex_ratio, female_bias)
      self.assertEqual(self.c.getFecundityP(), p)
      
   def testBadSetReduction(self):
      with self.assertRaises(ValueError) as the_err:
         e = CellDynamicsMosquito26Delay(num_species = 14, delay = 17, current_index = 7, death_rate = [[1.0] * 14] * 6, competition = [[0] * 14] * 14, emergence_rate = [1.0] * 14, activity = [[0] * 14] * 14, reduction = [list(range(6))] * 5)
      self.assertEqual(str(the_err.exception), "size of reduction, 5, must be equal to 6")
      with self.assertRaises(ValueError) as the_err:
         e = CellDynamicsMosquito26Delay(num_species = 2, delay = 17, current_index = 7, death_rate = [[1.0] * 2] * 6, competition = [[1, 2], [3, 4]], emergence_rate = [1.0] * 2, activity = [[1, 2], [2, 1]], reduction = [list(range(5))] * 6)
      self.assertEqual(str(the_err.exception), "size of reduction[0], 5, must be equal to 6")

      with self.assertRaises(ValueError) as the_err:
         self.c.setReduction([2, 3, 4, 5])
      self.assertEqual(str(the_err.exception), "size of reduction, 4, must be equal to 6")
      with self.assertRaises(ValueError) as the_err:
         self.c.setReduction([[0] * 6, [1], [2], [3], [4], [5]])
      self.assertEqual(str(the_err.exception), "size of reduction[1], 1, must be equal to 6")

   def testSetGetReduction(self):
      self.assertTrue(arrayequal(self.d.getReduction(), self.red))
      a = [[3, 4, 5, 1, 2, 3], [1, 5, 6, 6, 5, 1], [8, 8, 1, 9, 3, 1], list(range(1, 7)), [4, 5, 6, 4, 4, 4], [3, 6, 1, 7, 4, 6]]
      self.c.setReduction(a)
      self.assertTrue(arrayequal(self.c.getReduction(), a))

   def testBadSetHybridisation(self):
      with self.assertRaises(ValueError) as the_err:
         e = CellDynamicsMosquito26Delay(num_species = 14, delay = 17, current_index = 7, death_rate = [[1.0] * 14] * 6, competition = [[0] * 14] * 14, emergence_rate = [1.0] * 14, activity = [[0] * 14] * 14, reduction = [list(range(6))] * 6, hybridisation = [[list(range(14))] * 14] * 12)
      self.assertEqual(str(the_err.exception), "size of hybridisation, 12, must be equal to 14")
      with self.assertRaises(ValueError) as the_err:
         e = CellDynamicsMosquito26Delay(num_species = 2, delay = 17, current_index = 7, death_rate = [[1.0] * 2] * 6, competition = [[1, 2], [3, 4]], emergence_rate = [1.0] * 2, activity = [[1, 2], [2, 1]], reduction = [list(range(6))] * 6, hybridisation = [[list(range(14))] * 11] * 2)
      self.assertEqual(str(the_err.exception), "size of hybridisation[0], 11, must be equal to 2")
      with self.assertRaises(ValueError) as the_err:
         e = CellDynamicsMosquito26Delay(num_species = 2, delay = 17, current_index = 7, death_rate = [[1.0] * 2] * 6, competition = [[1, 2], [3, 4]], emergence_rate = [1.0] * 2, activity = [[1, 2], [2, 1]], reduction = [list(range(6))] * 6, hybridisation = [[[1, 2], [1, 2]], [[3333], [1, 2]]])
      self.assertEqual(str(the_err.exception), "size of hybridisation[1][0], 1, must be equal to 2")

      with self.assertRaises(ValueError) as the_err:
         self.c.setHybridisation([2, 3, 4, 5])
      self.assertEqual(str(the_err.exception), "size of hybridisation, 4, must be equal to 3")
      with self.assertRaises(ValueError) as the_err:
         self.c.setHybridisation([[[0, 1, 2], [3, 2, 1], [2, 3, 1]], [[0, 1, 2], [3, 2, 1], [111111, 2222]], [[0, 1, 2], [3, 2, 1], [2, 3, 1]]])
      self.assertEqual(str(the_err.exception), "size of hybridisation[1][2], 2, must be equal to 3")

   def testSetGetHybridisation(self):
      self.assertTrue(arrayequal(self.d.getReduction(), self.red))
      a = [[[3, 4, 5], [1, 2, 3], [1, 5, 6]], [[6, 5, 1], [8, 8, 1], [9, 3, 1]], [[4, 5, 6], [4, 4, 4], [3, 6, 1]]]
      self.c.setHybridisation(a)
      self.assertTrue(arrayequal(self.c.getHybridisation(), a))

   def testGetInheritance(self):
      m_w = 0.125
      m_c = 0.3
      d = CellDynamicsMosquito26Delay(m_w = m_w, m_c = m_c)
      ic = self.c.getInheritance()
      ic = d.getInheritance()

      w = 0.5 * (1 - m_w)
      c = 0.5 * (1 - m_c)
      r = 0.5 * (m_w + m_c)
      gold = [[[4*w**2,0,4*(1-2*w)*w,0,0,(1-2*w)**2],
               [2*w**2,2*c*w,(1+2*r-2*w)*w,0,c-2*c*w,r-2*r*w],
               [2*w**2,0,(3-4*w)*w,0,0,(-1+w)*(-1+2*w)],
               [0,4*c*w,2*(1-2*c)*w,0,2*c*(1-2*w),(-1+2*c)*(-1+2*w)],
               [0,2*c*w,-2*(-1+c)*w,0,c-2*c*w,(-1+c)*(-1+2*w)],
               [0,0,2*w,0,0,1-2*w]],
              [[2*w**2,2*c*w,(1+2*r-2*w)*w,0,c-2*c*w,r-2*r*w],
               [w**2,2*c*w,2*r*w,c**2,2*c*r,r**2],
               [w**2,c*w,(1+r-w)*w,0,c-c*w,r-r*w],
               [0,2*c*w,w-2*c*w,2*c**2,c*(1-2*c+2*r),r-2*c*r],
               [0,c*w,w-c*w,c**2,c*(1-c+r),r-c*r],
               [0,0,w,0,c,r]],
              [[2*w**2,0,(3-4*w)*w,0,0,(-1+w)*(-1+2*w)],
               [w**2,c*w,(1+r-w)*w,0,c-c*w,r-r*w],
               [w**2,0,-2*(-1+w)*w,0,0,(-1+w)**2],
               [0,2*c*w,w-2*c*w,0,-2*c*(-1+w),(-1+2*c)*(-1+w)],
               [0,c*w,w-c*w,0,c-c*w,(-1+c)*(-1+w)],
               [0,0,w,0,0,1-w]],
              [[0,4*c*w,2*(1-2*c)*w,0,2*c*(1-2*w),(-1+2*c)*(-1+2*w)],
               [0,2*c*w,w-2*c*w,2*c**2,c*(1-2*c+2*r),r-2*c*r],
               [0,2*c*w,w-2*c*w,0,-2*c*(-1+w),(-1+2*c)*(-1+w)],
               [0,0,0,4*c**2,4*(1-2*c)*c,(1-2*c)**2],
               [0,0,0,2*c**2,(3-4*c)*c,(-1+c)*(-1+2*c)],
               [0,0,0,0,2*c,1-2*c]],
              [[0,2*c*w,-2*(-1+c)*w,0,c-2*c*w,(-1+c)*(-1+2*w)],
               [0,c*w,w-c*w,c**2,c*(1-c+r),r-c*r],
               [0,c*w,w-c*w,0,c-c*w,(-1+c)*(-1+w)],
               [0,0,0,2*c**2,(3-4*c)*c,(-1+c)*(-1+2*c)],
               [0,0,0,c**2,-2*(-1+c)*c,(-1+c)**2],
               [0,0,0,0,c,1-c]],
              [[0,0,2*w,0,0,1-2*w],
               [0,0,w,0,c,r],
               [0,0,w,0,0,1-w],
               [0,0,0,0,2*c,1-2*c],
               [0,0,0,0,c,1-c],
               [0,0,0,0,0,1]]]
      for gM in range(6):
         for gF in range(6):
            for g in range(6):
               self.assertTrue(abs(gold[gM][gF][g] - ic[gM + gF * 6 + g * 36]) < 1E-6)
 
 
   def testSetGetMinCarryingCapacity(self):
      self.assertEqual(self.d.getMinCarryingCapacity(), 0.125)
      self.c.setMinCarryingCapacity(123)
      self.assertEqual(self.c.getMinCarryingCapacity(), 123)
      
      

   def testEvolve1(self):
      # Tests evolve doesn't crash

      initial_condition = [random.random() for i in range(self.c.getNumberOfPopulations() + self.c.getNumberOfParameters())]
      pap = array.array('f', initial_condition)
      dt = 1.23E-2
      self.c.evolve(dt, pap)
      self.c.incrementCurrentIndex() # not necessary here: just good practice to increment after evolve has been called for all grid cells

      self.assertTrue(True)
      

   def testEvolve2(self):
      # Test carrying_capacity = 0 case
      delay = 5
      cd = CellDynamicsMosquito26Delay(delay = delay, current_index = 2)

      # initialise populations and carrying capacities
      initial_condition = [random.random() for i in range(cd.getNumberOfPopulations())] + [0, 0, 0] # last 3 are the carrying capacities for the 3-species case
      pap = array.array('f', initial_condition)
      # set the death rates
      dr = [[random.random() for species in range(3)] for genotype in range(6)]
      cd.setDeathRate(dr)

      dt = 1.23E-2
      for timestep in range(13):
         # form the expected answer, and put in "gold"
         gold = list(pap) # result from previous timestep
         current_index = cd.getCurrentIndex()
         current_base = current_index * 3 * 6 * 2
         new_adults = [0 for i in range(3 * 6 * 2)]
         for sex in range(2):
            for genotype in range(6):
               for species in range(3):
                  ind = species + genotype * 3 + sex * 3 * 6
                  new_adults[ind] = gold[current_base + ind] * exp(- dr[genotype][species] * dt)
         for sex in range(2):
            for genotype in range(6):
               for species in range(3):
                  ind = species + genotype * 3 + sex * 3 * 6
                  gold[(current_index + 1) % (delay + 1) * 3 * 6 * 2 + ind] = new_adults[ind]

         # get the code to evolve and check answer is gold
         cd.evolve(dt, pap)
         cd.incrementCurrentIndex() # note: current_index must be incremented after evolve has been called for all grid cells (in this case, there is no spatial structure, ie, no grid cells)

         self.assertTrue(arrayfuzzyequal(pap, gold, 1E-6))
      
if __name__ == '__main__':
   unittest.main()

