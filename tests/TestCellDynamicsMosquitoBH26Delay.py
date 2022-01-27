import os
import sys
import unittest
import array
import random
from math import exp

# so we can find our ../code no matter how we are called
findbin = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(findbin + "/../code")

from cellDynamics import CellDynamicsMosquitoBH26Delay

def arrayequal(a, b):
   return (len(a) == len(b)) and all([a[i] == b[i] for i in range(min(len(a), len(b)))])

def arrayfuzzyequal(a, b, eps):
   return (len(a) == len(b)) and all([(a[i] > b[i] - eps and a[i] < b[i] + eps) for i in range(min(len(a), len(b)))])

class TestCellDynamicsMosquitoBH26Delay(unittest.TestCase):

   def setUp(self):
      self.c = CellDynamicsMosquitoBH26Delay()
      self.ide4 = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]
      self.red = [list(range(1, 7)), list(range(2, 8)), [1] * 6, list(range(4, 10)), [2] * 6, [1, 2, 1, 2, 1, 2]]
      self.hyb = [[[1, 2, 3, 4], [5, 4, 3, 2], [3, 3, 2, 1], [5, 6, 7, 4]], [[6, 7, 3, 2], [4, 6, 1, 8], [4, 7, 3, 2], [9, 8, 0, 7]], [[5, 7, 2, 9], [4, 5, 3, 6], [4, 3, 1, 2], [4, 2, 6, 8]], [[6, 5, 7, 4], [5, 6, 3, 7], [5, 9, 0, 1], [0, 8, 0, 4]]]
      self.d = CellDynamicsMosquitoBH26Delay(num_species = 4, delay = 17, current_index = 7, death_rate = [[1.0] * 4] * 6, competition = self.ide4, emergence_rate = [1.0] * 4, activity = self.ide4, reduction = self.red, hybridisation = self.hyb)

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

   def testSetGetDeathRate(self):
      self.assertTrue(arrayequal(self.d.getDeathRate(), [[1.0] * 4] * 6))
      dr = [[i, i + 6, i + 123] for i in range(1, 7)]
      self.c.setDeathRate(dr)
      self.assertTrue(arrayequal(self.c.getDeathRate(), dr))

   def testSetGetCompetition(self):
      self.assertTrue(arrayequal(self.d.getCompetition(), self.ide4))
      a = [[1, 2, 3], [4, 5, 6], [6, 7, 8]]
      self.c.setCompetition(a)
      self.assertTrue(arrayequal(self.c.getCompetition(), a))

   def testSetGetEmergenceRate(self):
      self.assertTrue(arrayequal(self.d.getEmergenceRate(), [1.0] * 4))
      self.d.setEmergenceRate([1, 2, 3, 4])
      self.assertTrue(arrayequal(self.d.getEmergenceRate(), [1, 2, 3, 4]))

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
      
   def testSetGetReduction(self):
      self.assertTrue(arrayequal(self.d.getReduction(), self.red))
      a = [[3, 4, 5, 1, 2, 3], [1, 5, 6, 6, 5, 1], [8, 8, 1, 9, 3, 1], list(range(1, 7)), [4, 5, 6, 4, 4, 4], [3, 6, 1, 7, 4, 6]]
      self.c.setReduction(a)
      self.assertTrue(arrayequal(self.c.getReduction(), a))

   def testSetGetHybridisation(self):
      self.assertTrue(arrayequal(self.d.getReduction(), self.red))
      a = [[[3, 4, 5], [1, 2, 3], [1, 5, 6]], [[6, 5, 1], [8, 8, 1], [9, 3, 1]], [[4, 5, 6], [4, 4, 4], [3, 6, 1]]]
      self.c.setHybridisation(a)
      self.assertTrue(arrayequal(self.c.getHybridisation(), a))

   def testGetInheritance(self):
      m_w = 0.125
      m_c = 0.3
      d = CellDynamicsMosquitoBH26Delay(m_w = m_w, m_c = m_c)
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
 
   def testEvolve1(self):
      # Tests evolve doesn't crash

      initial_condition = [random.random() for i in range(self.c.getNumberOfPopulations() + self.c.getNumberOfParameters())]
      pap = array.array('f', initial_condition)
      dt = 1.23
      self.c.evolve(dt, pap)
      self.c.incrementCurrentIndex() # not necessary here: just good practice to increment after evolve has been called for all grid cells

      self.assertTrue(True)
      

   def testEvolve2(self):
      # Test qm = 0 case
      delay = 5
      cd = CellDynamicsMosquitoBH26Delay(delay = delay, current_index = 2)

      # initialise populations and carrying capacities
      initial_condition = [random.random() for i in range(cd.getNumberOfPopulations())] + [0, 0, 0] # last 3 are the "qm" for the 3-species case
      pap = array.array('f', initial_condition)
      # set the death rates
      dr = [[random.random() for species in range(3)] for genotype in range(6)]
      cd.setDeathRate(dr)

      dt = 1.23
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
      
   def testEvolve3(self):
      # Test emergence = 0 case, where the ODE is just deaths-only
      delay = 3
      cd = CellDynamicsMosquitoBH26Delay(delay = delay, emergence_rate = [0, 0, 0])

      # initialise populations and carrying capacities
      initial_condition = [random.random() for i in range(cd.getNumberOfPopulations())] + [1E4, 2E4, 3E4] # last 3 are "qm"
      pap = array.array('f', initial_condition)
      # set the death rates
      dr = [[random.random() for species in range(3)] for genotype in range(6)]
      cd.setDeathRate(dr)

      dt = 1.23
      for timestep in range(4):
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
      
   def testEvolve_zeroMales(self):
      # Test case where the are zero males, where the ODE reduces to just deaths-only
      num_species = 4
      dr = [[random.random() for species in range(num_species)] for genotype in range(6)]
      self.d.setDeathRate(dr)

      # initialise populations and carrying capacities
      initial_condition = [random.random() for i in range(self.d.getNumberOfPopulations())] + [1E4, 2E4, 3E4, 4E4] # last num_species are the "qm"
      small_value = 0.004
      self.d.setSmallValue(small_value) # any population less than this will be zeroed
      for delay in range(self.d.getDelay() + 1):
         for species in range(num_species):
            for genotype in range(6):
               for sex in [0]: # males only
                  ind = species + genotype * num_species + sex * num_species * 6 + delay * num_species * 6 * 2
                  initial_condition[ind] = small_value * 1E-5  # this is zero males because of SmallValue
      pap = array.array('f', initial_condition)

      dt = 1.23
      for timestep in range(2):
         # form the expected answer, and put in "gold"
         gold = list(pap) # result from previous timestep
         current_index = self.d.getCurrentIndex()
         current_base = current_index * num_species * 6 * 2
         new_adults = [0 for i in range(num_species * 6 * 2)]
         for sex in range(2):
            for genotype in range(6):
               for species in range(num_species):
                  ind = species + genotype * num_species + sex * num_species * 6
                  new_adults[ind] = gold[current_base + ind] * exp(- dr[genotype][species] * dt)
         for sex in range(2):
            for genotype in range(6):
               for species in range(num_species):
                  ind = species + genotype * num_species + sex * num_species * 6
                  if new_adults[ind] > small_value:
                     gold[(current_index + 1) % (delay + 1) * num_species * 6 * 2 + ind] = new_adults[ind]
                  else:
                     gold[(current_index + 1) % (delay + 1) * num_species * 6 * 2 + ind] = 0.0

         # get the code to evolve and check answer is gold
         self.d.evolve(dt, pap)
         self.d.incrementCurrentIndex() # note: current_index must be incremented after evolve has been called for all grid cells (in this case, there is no spatial structure, ie, no grid cells)

         self.assertTrue(arrayfuzzyequal(pap, gold, 1E-6))
      
   def testEvolve_zeroFemales(self):
      # Test case where the are zero females, where the ODE reduces to just deaths-only
      num_species = 4
      dr = [[random.random() for species in range(num_species)] for genotype in range(6)]
      self.d.setDeathRate(dr)

      # initialise populations and carrying capacities
      initial_condition = [random.random() for i in range(self.d.getNumberOfPopulations())] + [1E4, 2E4, 3E4, 4E4] # last num_species are the "qm"
      small_value = 1E-3
      self.d.setSmallValue(small_value) # any population much less than this gets zeroed
      for delay in range(self.d.getDelay() + 1):
         for species in range(num_species):
            for genotype in range(6):
               for sex in [1]: # females only
                  ind = species + genotype * num_species + sex * num_species * 6 + delay * num_species * 6 * 2
                  initial_condition[ind] = 1E-10 # zero females (value is less than small_value)
      pap = array.array('f', initial_condition)

      dt = 1.23
      for timestep in range(11):
         # form the expected answer, and put in "gold"
         gold = list(pap) # result from previous timestep
         current_index = self.d.getCurrentIndex()
         current_base = current_index * num_species * 6 * 2
         new_adults = [0 for i in range(num_species * 6 * 2)]
         for sex in range(2):
            for genotype in range(6):
               for species in range(num_species):
                  ind = species + genotype * num_species + sex * num_species * 6
                  new_adults[ind] = gold[current_base + ind] * exp(- dr[genotype][species] * dt)
         for sex in range(2):
            for genotype in range(6):
               for species in range(num_species):
                  ind = species + genotype * num_species + sex * num_species * 6
                  if new_adults[ind] > small_value:
                     gold[(current_index + 1) % (delay + 1) * num_species * 6 * 2 + ind] = new_adults[ind]
                  else:
                     gold[(current_index + 1) % (delay + 1) * num_species * 6 * 2 + ind] = 0.0

         # get the code to evolve and check answer is gold
         self.d.evolve(dt, pap)
         self.d.incrementCurrentIndex() # note: current_index must be incremented after evolve has been called for all grid cells (in this case, there is no spatial structure, ie, no grid cells)

         self.assertTrue(arrayfuzzyequal(pap, gold, 1E-6))
      
   def testEvolve_oneSpecies(self):
      # Test case where there is just one species, and random parameters and random initial conditions
      num_species = 1
      delay = 5
      current_index = 2
      death_rate = [[random.random()] for i in range(6)]
      competition = random.random()
      emergence_rate = random.random()
      activity = random.random()
      reduction = [[random.random() for i in range(6)] for j in range(6)]
      hybridisation = random.random()
      sex_ratio = random.random()
      female_bias = random.random()
      m_w = random.random()
      m_c = random.random()
      tiny = CellDynamicsMosquitoBH26Delay(num_species = num_species, delay = delay, current_index = current_index, death_rate = death_rate, competition = [[competition]], emergence_rate = [emergence_rate], activity = [[activity]], reduction = reduction, hybridisation = [[[hybridisation]]], sex_ratio = sex_ratio, female_bias = female_bias, m_w = m_w, m_c = m_c)

      # precalculate O * R.  Since other tests have shown inheritance, fecundity and reduction are OK, just use them here
      ic = [[[tiny.getInheritance()[gM + gF * 6 + g * 36] for g in range(6)] for gF in range(6)] for gM in range(6)]
      pp = tiny.getFecundityP()
      rr = tiny.getReduction()
      o_times_r = [[[[ic[gM][gF][g] * pp[gM][gF][s] * rr[gM][gF] for s in range(2)] for g in range(6)] for gF in range(6)] for gM in range(6)]

      # initialise populations and carrying capacities
      qm = 10 + random.random()
      initial_condition = [random.random() for i in range(tiny.getNumberOfPopulations())] + [qm]
      pap = array.array('f', initial_condition)

      dt = 1.23
      for timestep in range(11):
         # form the expected answer, and put in "gold"
         gold = list(pap) # result from previous timestep
         current_index = tiny.getCurrentIndex()
         current_base = current_index * num_species * 6 * 2
         delayed_base = (current_index + 1) % (delay + 1) * num_species * 6 * 2

         # for num_species = 1, the proportinate-mixing is just:
         denom = sum([gold[delayed_base + gprime] for gprime in range(6)])
         xprimeM = [gold[delayed_base + g] / denom for g in range(6)]

         # compute y[s][g] and C
         y = [[sum([sum([hybridisation * o_times_r[gM][gF][g][s] * gold[delayed_base + gF + 1 * 6] * xprimeM[gM] for gF in range(6)]) for gM in range(6)]) for g in range(6)] for s in range(2)]
         cc = competition * sum([sum(y[s]) for s in range(2)])

         # compute B[s][g]
         bb = [[qm * emergence_rate * y[s][g] / (qm + cc) for g in range(6)] for s in range(2)]

         # compute the new adult populations
         new_adults = [0 for i in range(num_species * 6 * 2)]
         for sex in range(2):
            for genotype in range(6):
               for species in range(num_species):
                  ind = species + genotype * num_species + sex * num_species * 6
                  new_adults[ind] = bb[sex][genotype] / death_rate[genotype][species] + (gold[current_base + ind] - bb[sex][genotype] / death_rate[genotype][species]) * exp(- death_rate[genotype][species] * dt)

         # put the result into the correct slots in pops_and_params
         for sex in range(2):
            for genotype in range(6):
               for species in range(num_species):
                  ind = species + genotype * num_species + sex * num_species * 6
                  gold[delayed_base + ind] = new_adults[ind]

         # get the code to evolve and check answer is gold
         tiny.evolve(dt, pap)
         tiny.incrementCurrentIndex() # note: current_index must be incremented after evolve has been called for all grid cells (in this case, there is no spatial structure, ie, no grid cells)

         self.assertTrue(arrayfuzzyequal(pap, gold, 1E-6))
      
   def testEvolve_onlyWild(self):
      # Test that when there are only wild-types and m_w=0, no wild-types are produced
      num_species = 4
      wild = CellDynamicsMosquitoBH26Delay(num_species = num_species, delay = 7, current_index = 2, death_rate = [[1.0] * 4] * 6, competition = self.ide4, emergence_rate = [1.0] * 4, activity = self.ide4, reduction = self.red, hybridisation = self.hyb, m_w = 0)

      # initialise populations and qm values
      initial_condition = [random.random() for i in range(wild.getNumberOfPopulations())] + [10 + random.random() for i in range(wild.getNumberOfParameters())]
      for delay in range(wild.getDelay() + 1):
         for species in range(num_species):
            for genotype in range(6):
               for sex in range(2):
                  ind = species + genotype * num_species + sex * num_species * 6 + delay * num_species * 6 * 2
                  if genotype != 0:
                     initial_condition[ind] = 0.0 # zero non-wildtype
      pap = array.array('f', initial_condition)
      
      dt = 1.23
      for timestep in range(11):
         wild.evolve(dt, pap)
         wild.incrementCurrentIndex() # note: current_index must be incremented after evolve has been called for all grid cells (in this case, there is no spatial structure, ie, no grid cells)
         for delay in range(wild.getDelay() + 1):
            for species in range(num_species):
               for genotype in range(6):
                  for sex in range(2):
                     ind = species + genotype * num_species + sex * num_species * 6 + delay * num_species * 6 * 2
                     if genotype != 0:
                        self.assertEqual(pap[ind], 0.0)
      
      
   def testEvolve_oneSpeciesWW(self):
      # Test long-term behaviour where there is just one species and wild-type
      num_species = 1
      delay = 3
      current_index = 0
      death_rate = [[random.random()] for i in range(6)]
      #death_rate[0][0] = 0.75
      competition = random.random()
      emergence_rate = 2.01 + random.random() # to ensure a nonzero steady-state
      activity = random.random()
      reduction = [[random.random() for i in range(6)] for j in range(6)]
      reduction[0][0] = 1.0 + random.random() # to ensure a nonzero steady-state
      hybridisation = 1.0 + random.random() # to ensure a nonzero steady-state
      sex_ratio = random.random()
      female_bias = random.random()
      m_w = 0.0 # so only ww mosquitos
      m_c = random.random()
      tiny = CellDynamicsMosquitoBH26Delay(num_species = num_species, delay = delay, current_index = current_index, death_rate = death_rate, competition = [[competition]], emergence_rate = [emergence_rate], activity = [[activity]], reduction = reduction, hybridisation = [[[hybridisation]]], sex_ratio = sex_ratio, female_bias = female_bias, m_w = m_w, m_c = m_c)
      qm = 1.0 + random.random()

      # calculate steady-state adult population
      ic = 1
      pp = 0.5
      rr = reduction[0][0] # >= 1
      hh = hybridisation # >= 1
      HORXprimeM = hh * ic * pp * rr # >= 1/2
      lambda_m = emergence_rate # > 2
      # note that because lambda_m > 2 and death_rate <= 1, the steady_state is positive
      num_sexes = 2
      steady_state = (1.0 / (num_sexes * competition) / HORXprimeM) * qm * (HORXprimeM * lambda_m / death_rate[0][0] - 1)


      # initialise very small populations
      initial_condition = [1E-6 * random.random() * steady_state for i in range(tiny.getNumberOfPopulations())] + [qm]
      for delay in range(tiny.getDelay() + 1):
         for species in range(num_species):
            for genotype in range(6):
               for sex in range(2):
                  ind = species + genotype * num_species + sex * num_species * 6 + delay * num_species * 6 * 2
                  if genotype != 0:
                     initial_condition[ind] = 0.0 # zero non-wildtype
      pap = array.array('f', initial_condition)

      dt = 1.23
      for timestep in range(10000):
         # get the code to evolve and check answer is gold
         tiny.evolve(dt, pap)
         tiny.incrementCurrentIndex() # note: current_index must be incremented after evolve has been called for all grid cells (in this case, there is no spatial structure, ie, no grid cells)
         error = abs(1 - pap[6] / steady_state)
         if error < 1E-3:
            break
      self.assertTrue(error < 1E-3)

      # initialise very large populations
      initial_condition = [1E6 * (1 + random.random()) * steady_state for i in range(tiny.getNumberOfPopulations())] + [qm]
      for delay in range(tiny.getDelay() + 1):
         for species in range(num_species):
            for genotype in range(6):
               for sex in range(2):
                  ind = species + genotype * num_species + sex * num_species * 6 + delay * num_species * 6 * 2
                  if genotype != 0:
                     initial_condition[ind] = 0.0 # zero non-wildtype
      pap = array.array('f', initial_condition)

      dt = 1.23
      for timestep in range(10000):
         # get the code to evolve and check answer is gold
         tiny.evolve(dt, pap)
         tiny.incrementCurrentIndex() # note: current_index must be incremented after evolve has been called for all grid cells (in this case, there is no spatial structure, ie, no grid cells)
         error = abs(1 - pap[6] / steady_state)
         if error < 1E-3:
            break
      self.assertTrue(error < 1E-3)

   def testGetSetSmallValue(self):
      self.assertEqual(self.c.getSmallValue(), 0.0)
      self.c.setSmallValue(1.0)
      self.assertEqual(self.c.getSmallValue(), 1.0)

if __name__ == '__main__':
   unittest.main()

