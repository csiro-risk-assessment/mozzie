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
   return (len(a) == len(b)) and all([a[i] == b[i] for i in range(min(len(a), len(b)))])

def arrayfuzzyequal(a, b, eps):
   if len(a) != len(b):
      return False
   abs_tol = 1E-4 * eps
   for i in range(len(a)):
      if (abs(b[i]) < abs_tol): # if very small then look at absolute tolerance
         if (a[i] < b[i] - abs_tol) or (a[i] > b[i] + abs_tol):
            return False
      else: # do relative tolerance
         if abs(a[i] / b[i] - 1.0) > eps:
            return False
   return True

def h(m, mm, mf): # m = offspring mosquito_type, mm = Male mosquito_type, mf = Female mosquito_type
   if mm == mf:
      if m == mm:
         return 1 # all same species
      else:
         return 0 # cannot create offspring of different type
   else:
      if m == mm:
         return 0.5
      elif m == mf:
         return 0.5
   return 0.0

def Wm(matingcomp, mm, mf): # mm = Male mosquito_type, mf = Female mosquito_type
   if mm == mf: # offspring is the male mosquito type
      return 1
   return matingcomp # 0.01 value in "26" constructor

def Wg(g): # g = genotype of male
   h_e = h_n = 0.5
   s_e = 0.1
   s_n = 0.05
   if g == 0: #ww
      return 1.0
   elif g == 1: #wc
      return (1.0 - h_e * s_e) * (1.0 - h_n * s_n)
   elif g == 2: #wr
      return 1.0
   elif g == 3: #cc
      return (1.0 - s_e) * (1.0 - s_n)
   elif g == 4: #cr
      return (1.0 - h_e * s_e) * (1.0 - h_n * s_n)
   else: #rr
      return 1.0

def mm(matingcomp, m, g, mM, mF, gM, X): # m=offspring mosquito_type, g=offspring geonotype, mM = Male mosquito_type, mF = Female mosquito_type, gM = Male genotype, X = pap
   denom = 0.0
   for gp in range(6):
      for mp in range(2):
         ind = mp + 2 * (gp + 6 * (0 + 2 * 1))
         denom += Wm(matingcomp, mp, mF) * Wg(gp) * X[ind]
   ind = mM + 2 * (gM + 6 * (0 + 2 * 1))
   num = h(m, mM, mF) * Wm(matingcomp, mM, mF) * Wg(gM) * X[ind]
   if denom == 0:
      if num == 0:
         return 0.0
      else:
         sys.stderr.write("denom == 0 but num != 0\n")
         sys.exit(1)
   return num / denom


def bb_no_fecundity(carrying, matingcomp, s, g, m, X): # carrying[m]=carrying_cap, matingcomp=w, s=sex, g=genotype, m=mosquito_type, X=pap
   term = 0.0
   for gM in range(6): # male genotype
      for gF in range(6): # female genotype
         for mM in range(2): # male mosquito_type
            for mF in range(2): # female mosquito_type
               ind = mF + 2 * (gF + 6 * (1 + 2 * 1))
               term += mm(matingcomp, m, g, mM, mF, gM, X) * OO(g, gM, gF) * X[ind]
   return max(0.0, 1 - cc(m, X) / carrying[m]) * term

def OO(g, gM, gF): # g=genotype, gM=Male genotype, gF=Female genotype
   p = 0.5 # no sex bias in "26"
   return p * inher(g, gM, gF)
def inher(g, gM, gF): # g=genotype, gM=Male genotype, gF=Female genotype
   p = 0.5
   kc = 0.995
   kj = 0.02
   kne = 0.0001
   w = 0.5 - 0.5 * kc
   c = 0.5 + kc * (1 - kj) * (1 - kne) * 0.5
   r = kc * (kne + kj * (1 - kne)) * 0.5

   # since matings are symmetrical, just code the 'lower-diagonal' part
   if gM <= gF:
      parent0 = gM
      parent1 = gF
   else:
      parent0 = gF
      parent1 = gM
   # now we know parent1 >= parent0

   if parent0 == 0: # Parent0 ww
      if parent1 == 0: # Parent1 ww
         if g == 0: # offspring ww
            return 1.0
         return 0.0
      elif parent1 == 1: # Parent1 wc
         if g == 0: # offspring ww
            return w
         elif g == 1: # offspring wc
            return c
         elif g == 2: # offpsring wr
            return r
         return 0.0
      elif parent1 == 2: # Parent1 wr
         if g == 0: # offspring ww
            return 0.5
         elif g == 1: # offspring wc
            return 0
         elif g == 2: # offpsring wr
            return 0.5
         return 0.0
      elif parent1 == 3: # Parent1 cc
         if g == 1: # offspring wc
            return 1
         return 0.0
      elif parent1 == 4: # Parent1 cr
         if g == 1: # offspring wc
            return 0.5
         elif g == 2: # offpsring wr
            return 0.5
         return 0.0
      else: # Parent1 rr
         if g == 2: # offspring wr
            return 1
         return 0.0

   elif parent0 == 1: # Parent0 wc
      if parent1 == 1: # Parent1 wc
         if g == 0: # offspring ww
            return w * w
         elif g == 1: # offspring wc
            return 2 * w * c
         elif g == 2: # offpsring wr
            return 2 * w * r
         elif g == 3: # offpsring cc
            return c * c
         elif g == 4: # offpsring cr
            return 2 * c * r
         else: # offspring rr
            return r * r
      elif parent1 == 2: # Parent1 wr
         if g == 0: # offspring ww
            return 0.5 * w
         elif g == 1: # offspring wc
            return 0.5 * c
         elif g == 2: # offpsring wr
            return 0.5 * (w + r)
         elif g == 3: # offpsring cc
            return 0.0
         elif g == 4: # offpsring cr
            return 0.5 * c
         else: # offspring rr
            return 0.5 * r
      elif parent1 == 3: # Parent1 cc
         if g == 0: # offspring ww
            return 0
         elif g == 1: # offspring wc
            return w
         elif g == 2: # offpsring wr
            return 0
         elif g == 3: # offpsring cc
            return c
         elif g == 4: # offpsring cr
            return r
         else: # offspring rr
            return 0
      elif parent1 == 4: # Parent1 cr
         if g == 0: # offspring ww
            return 0
         elif g == 1: # offspring wc
            return 0.5 * w
         elif g == 2: # offpsring wr
            return 0.5 * w
         elif g == 3: # offpsring cc
            return 0.5 * c
         elif g == 4: # offpsring cr
            return 0.5 * (c + r)
         else: # offspring rr
            return 0.5 * r
      else: # Parent1 rr
         if g == 0: # offspring ww
            return 0
         elif g == 1: # offspring wc
            return 0
         elif g == 2: # offpsring wr
            return w
         elif g == 3: # offpsring cc
            return 0
         elif g == 4: # offpsring cr
            return c
         else: # offspring rr
            return r

   elif parent0 == 2: # Parent0 wr
      if parent1 == 2: # Parent1 wr
         if g == 0: # offspring ww
            return 0.25
         elif g == 1: # offspring wc
            return 0
         elif g == 2: # offpsring wr
            return 0.5
         elif g == 3: # offpsring cc
            return 0
         elif g == 4: # offpsring cr
            return 0
         else: # offspring rr
            return 0.25
      elif parent1 == 3: # Parent1 cc
         if g == 0: # offspring ww
            return 0
         elif g == 1: # offspring wc
            return 0.5
         elif g == 2: # offpsring wr
            return 0.5
         elif g == 3: # offpsring cc
            return 0
         elif g == 4: # offpsring cr
            return 0
         else: # offspring rr
            return 0
      elif parent1 == 4: # Parent1 cr
         if g == 0: # offspring ww
            return 0
         elif g == 1: # offspring wc
            return 0.25
         elif g == 2: # offpsring wr
            return 0.25
         elif g == 3: # offpsring cc
            return 0
         elif g == 4: # offpsring cr
            return 0.25
         else: # offspring rr
            return 0.25
      else: # Parent1 rr
         if g == 0: # offspring ww
            return 0
         elif g == 1: # offspring wc
            return 0
         elif g == 2: # offpsring wr
            return 0.5
         elif g == 3: # offpsring cc
            return 0
         elif g == 4: # offpsring cr
            return 0
         else: # offspring rr
            return 0.5

   elif parent0 == 3: # Parent0 cc
      if parent1 == 3: # Parent1 cc
         if g == 0: # offspring ww
            return 0
         elif g == 1: # offspring wc
            return 0
         elif g == 2: # offpsring wr
            return 0
         elif g == 3: # offpsring cc
            return 1
         elif g == 4: # offpsring cr
            return 0
         else: # offspring rr
            return 0
      elif parent1 == 4: # Parent1 cr
         if g == 0: # offspring ww
            return 0
         elif g == 1: # offspring wc
            return 0
         elif g == 2: # offpsring wr
            return 0
         elif g == 3: # offpsring cc
            return 0.5
         elif g == 4: # offpsring cr
            return 0.5
         else: # offspring rr
            return 0
      else: # Parent1 rr
         if g == 0: # offspring ww
            return 0
         elif g == 1: # offspring wc
            return 0
         elif g == 2: # offpsring wr
            return 0
         elif g == 3: # offpsring cc
            return 0
         elif g == 4: # offpsring cr
            return 1
         else: # offspring rr
            return 0

   elif parent0 == 4: # Parent0 cr
      if parent1 == 4: # Parent1 cr
         if g == 0: # offspring ww
            return 0
         elif g == 1: # offspring wc
            return 0
         elif g == 2: # offpsring wr
            return 0
         elif g == 3: # offpsring cc
            return 0.25
         elif g == 4: # offpsring cr
            return 0.5
         else: # offspring rr
            return 0.25
      else: # Parent1 rr
         if g == 0: # offspring ww
            return 0
         elif g == 1: # offspring wc
            return 0
         elif g == 2: # offpsring wr
            return 0
         elif g == 3: # offpsring cc
            return 0
         elif g == 4: # offpsring cr
            return 0.5
         else: # offspring rr
            return 0.5

   else: # Parent0 rr
      if parent1 == 5: # Parent1 rr
         if g == 0: # offspring ww
            return 0
         elif g == 1: # offspring wc
            return 0
         elif g == 2: # offpsring wr
            return 0
         elif g == 3: # offpsring cc
            return 0
         elif g == 4: # offpsring cr
            return 0
         else: # offspring rr
            return 1
   sys.stderr.write("Ooops")
   sys.exit(1)

def lotka(m, mp): # m=mosquito_type, mp=mosquito_type
   if m == mp:
      return 1.0
   else:
      return 0.4 # setAlphaComponent in mosquito26 constructor
   
def cc(m, X): # alpha=lotkavolterra, m=mosquito_type, X=pap
   term = 0.0
   for mp in range(2): # mosquito_type
      a = 0 # juveniles only
      for s in range(2): # sex
         for g in range(6): # genotype
            ind = mp + 2 * (g + 6 * (s + 2 * a))
            term += lotka(m, mp) * X[ind]
   return term

def predict_answer(carrying, matingcomp, fecundity, mu_adult, mu_juv, aging, dt, pap):
   prediction = [0] * 2 * 6 * 2 * 2
   # juveniles
   for s in range(2): # sex
      for g in range(6): # genotype
         for m in range(2): # mosquito_type
            ind = m + 2 * (g + 6 * s)
            prediction[ind] = pap[ind] * (1 - dt * mu_juv) - dt * aging * pap[ind] + dt * bb_no_fecundity(carrying, matingcomp, s, g, m, pap) * fecundity
   # adults
   for s in range(2): # sex
      for g in range(6): # genotype
         for m in range(2): # mosquito_type
            ind = m + 2 * (g + 6 * (s + 2))
            juvenile_ind = m + 2 * (g + 6 * (s + 0))
            prediction[ind] = pap[ind] * (1 - dt * mu_adult) + dt * aging * pap[juvenile_ind]
   return prediction

ep_tol = 1E-5

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

   def testGetAdvectionClass(self):
      gold = [1] * 24
      self.assertTrue(arrayequal(self.c.getAdvectionClass(), gold))
      self.c.setAdvectionClass(25, 26)
      gold[1] = 26
      self.assertTrue(arrayequal(self.c.getAdvectionClass(), gold))
      self.c.setNumAges(5)
      gold[1] = 1
      self.assertTrue(arrayequal(self.c.getAdvectionClass(), gold))

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
      self.assertTrue(arrayfuzzyequal([self.c.getAlphaComponentFromPython(0, 1)], [0.4], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getAlphaComponentFromPython(1, 0)], [0.4], ep_tol))
      self.assertEqual(self.c.getAlphaComponentFromPython(1, 1), 1.0)
      with self.assertRaises(ValueError) as the_err:
         self.c.setAlphaComponent(0, 2, 0)
      self.assertEqual(str(the_err.exception), "sp0 0 and sp1 2 must be less than the number of species, 2")
      with self.assertRaises(ValueError) as the_err:
         self.c.setAlphaComponent(2, 1, 0)
      self.assertEqual(str(the_err.exception), "sp0 2 and sp1 1 must be less than the number of species, 2")

   def testSetGetFitnessComponents(self):
      h_e = 0.5
      h_n = 0.5
      s_e = 0.1
      s_n = 0.05
      self.assertTrue(arrayfuzzyequal([self.c.getFitnessComponentFromPython(0)], [1], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getFitnessComponentFromPython(1)], [(1 - h_e * s_e) * (1 - h_n * s_n)], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getFitnessComponentFromPython(2)], [1], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getFitnessComponentFromPython(3)], [(1 - s_e) * (1 - s_n)], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getFitnessComponentFromPython(4)], [(1 - h_e * s_e) * (1 - h_n * s_n)], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getFitnessComponentFromPython(5)], [1], ep_tol))

      with self.assertRaises(ValueError) as the_err:
         self.c.getFitnessComponentFromPython(6)
      self.assertEqual(str(the_err.exception), "Genotype 6 must be less than the number of genotypes, 6")

      h_e = 0.25
      h_n = 0.625
      s_e = 0.75
      s_n = 0.125
      self.c.setFitnessComponents26(h_e, h_n, s_e, s_n)
      self.assertTrue(arrayfuzzyequal([self.c.getFitnessComponentFromPython(0)], [1], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getFitnessComponentFromPython(1)], [(1 - h_e * s_e) * (1 - h_n * s_n)], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getFitnessComponentFromPython(2)], [1], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getFitnessComponentFromPython(3)], [(1 - s_e) * (1 - s_n)], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getFitnessComponentFromPython(4)], [(1 - h_e * s_e) * (1 - h_n * s_n)], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getFitnessComponentFromPython(5)], [1], ep_tol))

      self.c.setNumGenotypes(2, 5)
      with self.assertRaises(ValueError) as the_err:
         self.c.setFitnessComponents26(h_e, h_n, s_e, s_n)
      self.assertEqual(str(the_err.exception), "setFitnessComponents26 can only be used if the number of genotypes is 6")
      

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
      self.assertTrue(arrayfuzzyequal([self.c.getMatingComponentFromPython(1, 0)], [0.01], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getMatingComponentFromPython(0, 1)], [0.01], ep_tol))

   def testSetGetInheritance(self):

      # default values in constructor
      k_c = 0.995
      k_j = 0.02
      k_ne = 1E-4
      w = 0.5 - 0.5 * k_c
      c = 0.5 + 0.5 * k_c * (1 - k_j) * (1 - k_ne)
      r = 0.5 * k_c * (k_ne + k_j * (1 - k_ne))

      # father, mother, offspring
      # ww, wc, wr, cc, cr, rr
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 0, 0)], [1], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 0, 1)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 0, 2)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 0, 3)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 0, 4)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 0, 5)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 1, 0)], [w], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 1, 1)], [c], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 1, 2)], [r], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 1, 3)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 1, 4)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 1, 5)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 2, 0)], [0.5], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 2, 1)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 2, 2)], [0.5], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 2, 3)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 2, 4)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 2, 5)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 3, 0)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 3, 1)], [1], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 3, 2)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 3, 3)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 3, 4)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 3, 5)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 4, 0)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 4, 1)], [0.5], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 4, 2)], [0.5], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 4, 3)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 4, 4)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 4, 5)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 5, 0)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 5, 1)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 5, 2)], [1], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 5, 3)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 5, 4)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 5, 5)], [0], ep_tol))

      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 1, 0)], [w * w], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 1, 1)], [2 * w * c], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 1, 2)], [2 * w * r], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 1, 3)], [c * c], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 1, 4)], [2 * c * r], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 1, 5)], [r * r], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 2, 0)], [0.5 * w], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 2, 1)], [0.5 * c], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 2, 2)], [0.5 * (w + r)], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 2, 3)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 2, 4)], [0.5 * c], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 2, 5)], [0.5 * r], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 3, 0)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 3, 1)], [w], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 3, 2)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 3, 3)], [c], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 3, 4)], [r], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 3, 5)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 4, 0)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 4, 1)], [0.5 * w], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 4, 2)], [0.5 * w], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 4, 3)], [0.5 * c], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 4, 4)], [0.5 * (c + r)], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 4, 5)], [0.5 * r], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 5, 0)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 5, 1)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 5, 2)], [w], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 5, 3)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 5, 4)], [c], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 5, 5)], [r], ep_tol))

      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 2, 0)], [0.25], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 2, 1)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 2, 2)], [0.5], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 2, 3)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 2, 4)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 2, 5)], [0.25], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 3, 0)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 3, 1)], [0.5], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 3, 2)], [0.5], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 3, 3)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 3, 4)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 3, 5)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 4, 0)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 4, 1)], [0.25], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 4, 2)], [0.25], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 4, 3)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 4, 4)], [0.25], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 4, 5)], [0.25], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 5, 0)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 5, 1)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 5, 2)], [0.5], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 5, 3)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 5, 4)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 5, 5)], [0.5], ep_tol))

      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(3, 3, 0)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(3, 3, 1)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(3, 3, 2)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(3, 3, 3)], [1], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(3, 3, 4)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(3, 3, 5)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(3, 4, 0)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(3, 4, 1)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(3, 4, 2)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(3, 4, 3)], [0.5], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(3, 4, 4)], [0.5], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(3, 4, 5)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(3, 5, 0)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(3, 5, 1)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(3, 5, 2)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(3, 5, 3)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(3, 5, 4)], [1], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(3, 5, 5)], [0], ep_tol))

      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(4, 4, 0)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(4, 4, 1)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(4, 4, 2)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(4, 4, 3)], [0.25], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(4, 4, 4)], [0.5], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(4, 4, 5)], [0.25], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(4, 5, 0)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(4, 5, 1)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(4, 5, 2)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(4, 5, 3)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(4, 5, 4)], [0.5], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(4, 5, 5)], [0.5], ep_tol))

      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(5, 5, 0)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(5, 5, 1)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(5, 5, 2)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(5, 5, 3)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(5, 5, 4)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(5, 5, 5)], [1], ep_tol))

      # symmetry
      for male in range(6):
         for female in range(6):
            for offspring in range(6):
               self.assertEqual(self.c.getInheritanceFromPython(male, female, offspring), self.c.getInheritanceFromPython(female, male, offspring))

      # some new values
      k_c = 0.875
      k_j = 0.0625
      k_ne = 0.03125
      w = 0.5 - 0.5 * k_c
      c = 0.5 + 0.5 * k_c * (1 - k_j) * (1 - k_ne)
      r = 0.5 * k_c * (k_ne + k_j * (1 - k_ne))
      self.c.setInheritance26(k_c, k_j, k_ne)

      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 0, 0)], [1], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 0, 1)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 0, 2)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 0, 3)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 0, 4)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 0, 5)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 1, 0)], [w], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 1, 1)], [c], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 1, 2)], [r], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 1, 3)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 1, 4)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 1, 5)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 2, 0)], [0.5], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 2, 1)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 2, 2)], [0.5], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 2, 3)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 2, 4)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 2, 5)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 3, 0)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 3, 1)], [1], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 3, 2)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 3, 3)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 3, 4)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 3, 5)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 4, 0)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 4, 1)], [0.5], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 4, 2)], [0.5], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 4, 3)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 4, 4)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 4, 5)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 5, 0)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 5, 1)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 5, 2)], [1], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 5, 3)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 5, 4)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 5, 5)], [0], ep_tol))

      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 1, 0)], [w * w], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 1, 1)], [2 * w * c], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 1, 2)], [2 * w * r], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 1, 3)], [c * c], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 1, 4)], [2 * c * r], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 1, 5)], [r * r], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 2, 0)], [0.5 * w], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 2, 1)], [0.5 * c], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 2, 2)], [0.5 * (w + r)], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 2, 3)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 2, 4)], [0.5 * c], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 2, 5)], [0.5 * r], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 3, 0)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 3, 1)], [w], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 3, 2)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 3, 3)], [c], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 3, 4)], [r], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 3, 5)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 4, 0)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 4, 1)], [0.5 * w], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 4, 2)], [0.5 * w], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 4, 3)], [0.5 * c], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 4, 4)], [0.5 * (c + r)], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 4, 5)], [0.5 * r], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 5, 0)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 5, 1)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 5, 2)], [w], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 5, 3)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 5, 4)], [c], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 5, 5)], [r], ep_tol))

      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 2, 0)], [0.25], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 2, 1)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 2, 2)], [0.5], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 2, 3)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 2, 4)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 2, 5)], [0.25], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 3, 0)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 3, 1)], [0.5], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 3, 2)], [0.5], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 3, 3)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 3, 4)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 3, 5)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 4, 0)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 4, 1)], [0.25], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 4, 2)], [0.25], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 4, 3)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 4, 4)], [0.25], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 4, 5)], [0.25], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 5, 0)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 5, 1)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 5, 2)], [0.5], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 5, 3)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 5, 4)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 5, 5)], [0.5], ep_tol))

      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(3, 3, 0)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(3, 3, 1)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(3, 3, 2)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(3, 3, 3)], [1], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(3, 3, 4)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(3, 3, 5)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(3, 4, 0)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(3, 4, 1)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(3, 4, 2)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(3, 4, 3)], [0.5], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(3, 4, 4)], [0.5], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(3, 4, 5)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(3, 5, 0)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(3, 5, 1)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(3, 5, 2)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(3, 5, 3)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(3, 5, 4)], [1], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(3, 5, 5)], [0], ep_tol))

      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(4, 4, 0)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(4, 4, 1)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(4, 4, 2)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(4, 4, 3)], [0.25], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(4, 4, 4)], [0.5], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(4, 4, 5)], [0.25], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(4, 5, 0)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(4, 5, 1)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(4, 5, 2)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(4, 5, 3)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(4, 5, 4)], [0.5], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(4, 5, 5)], [0.5], ep_tol))

      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(5, 5, 0)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(5, 5, 1)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(5, 5, 2)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(5, 5, 3)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(5, 5, 4)], [0], ep_tol))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(5, 5, 5)], [1], ep_tol))

      # symmetry
      for male in range(6):
         for female in range(6):
            for offspring in range(6):
               self.assertEqual(self.c.getInheritanceFromPython(male, female, offspring), self.c.getInheritanceFromPython(female, male, offspring))

      # exception testing
      with self.assertRaises(ValueError) as the_err:
         self.c.getInheritanceFromPython(6, 0, 0)
      self.assertEqual(str(the_err.exception), "All genotypes, 6, 0, 0 must be less than the number of genotypes, 6")
      with self.assertRaises(ValueError) as the_err:
         self.c.getInheritanceFromPython(0, 6, 0)
      self.assertEqual(str(the_err.exception), "All genotypes, 0, 6, 0 must be less than the number of genotypes, 6")
      with self.assertRaises(ValueError) as the_err:
         self.c.getInheritanceFromPython(0, 0, 6)
      self.assertEqual(str(the_err.exception), "All genotypes, 0, 0, 6 must be less than the number of genotypes, 6")
      self.c.setNumGenotypes(2, 5)
      with self.assertRaises(ValueError) as the_err:
         self.c.setInheritance26(k_c, k_j, k_ne)
      self.assertEqual(str(the_err.exception), "setInheritance26 can only be used if the number of genotypes is 6")


   def testEvolveSpecies1WW(self):
      dt = 3.0
      mu_adult = 0.125
      self.c.setMuAdult(mu_adult)
      mu_juv = 0.25
      self.c.setMuLarvae(mu_juv)
      fecundity = 8
      self.c.setFecundity(fecundity)
      carrying_cap = 60.0 # because we start with zero juveniles, this is actually irrelevant
      aging = 0.0 # because we start with zero juveniles, this is actually irrelevant

      # start with zero juveniles
      male_juv_ww = male_juv_wc = male_juv_wr = male_juv_cc = male_juv_cr = male_juv_rr = female_juv_ww = female_juv_wc = female_juv_wr = female_juv_cc = female_juv_cr = female_juv_rr = 0.0

      # species0 adults:
      sp0_male_ad_ww = sp0_male_ad_wc = sp0_male_ad_wr = sp0_male_ad_cc = sp0_male_ad_cr = sp0_male_ad_rr = sp0_female_ad_ww = sp0_female_ad_wc = sp0_female_ad_wr = sp0_female_ad_cc = sp0_female_ad_cr = sp0_female_ad_rr = 0.0
      #sp0_male_ad_ww = sp0_female_ad_ww = 20.0
      #sp0_male_ad_wc = sp0_female_ad_wc = 10.0
      ic_species0 = [male_juv_ww, male_juv_wc, male_juv_wr, male_juv_cc, male_juv_cr, male_juv_rr, female_juv_ww, female_juv_wc, female_juv_wr, female_juv_cc, female_juv_cr, female_juv_rr, sp0_male_ad_ww, sp0_male_ad_wc, sp0_male_ad_wr, sp0_male_ad_cc, sp0_male_ad_cr, sp0_male_ad_rr, sp0_female_ad_ww, sp0_female_ad_wc, sp0_female_ad_wr, sp0_female_ad_cc, sp0_female_ad_cr, sp0_female_ad_rr]

      # species1 adults:
      sp1_male_ad_ww = sp1_male_ad_wc = sp1_male_ad_wr = sp1_male_ad_cc = sp1_male_ad_cr = sp1_male_ad_rr = sp1_female_ad_ww = sp1_female_ad_wc = sp1_female_ad_wr = sp1_female_ad_cc = sp1_female_ad_cr = sp1_female_ad_rr = 0.0
      sp1_male_ad_ww = 20.0
      sp1_female_ad_ww = 40.0
      ic_species1 = [male_juv_ww, male_juv_wc, male_juv_wr, male_juv_cc, male_juv_cr, male_juv_rr, female_juv_ww, female_juv_wc, female_juv_wr, female_juv_cc, female_juv_cr, female_juv_rr, sp1_male_ad_ww, sp1_male_ad_wc, sp1_male_ad_wr, sp1_male_ad_cc, sp1_male_ad_cr, sp1_male_ad_rr, sp1_female_ad_ww, sp1_female_ad_wc, sp1_female_ad_wr, sp1_female_ad_cc, sp1_female_ad_cr, sp1_female_ad_rr]

      # build initial condition in the standard form
      initial_condition = []
      for i in range(24):
         initial_condition.append(ic_species0[i])
         initial_condition.append(ic_species1[i])
      initial_condition += [carrying_cap, carrying_cap]
      pap = array.array('f', initial_condition)

      # predict the answer
      prediction = predict_answer([carrying_cap, carrying_cap], 0.0, fecundity, mu_adult, mu_juv, aging, dt, pap)

      # evolve according to the cython code
      self.c.evolve(dt, pap)

      # check answer
      self.assertTrue(arrayfuzzyequal(prediction, pap[:48], ep_tol))


   def testEvolveSpecies0WW(self):
      dt = 3.0
      mu_adult = 0.125
      self.c.setMuAdult(mu_adult)
      mu_juv = 0.25
      self.c.setMuLarvae(mu_juv)
      fecundity = 8
      self.c.setFecundity(fecundity)
      carrying_cap = 60.0 # because we start with zero juveniles, this is actually irrelevant
      aging = 0.0 # because we start with zero juveniles, this is actually irrelevant

      # start with zero juveniles
      male_juv_ww = male_juv_wc = male_juv_wr = male_juv_cc = male_juv_cr = male_juv_rr = female_juv_ww = female_juv_wc = female_juv_wr = female_juv_cc = female_juv_cr = female_juv_rr = 0.0

      # species0 adults:
      sp0_male_ad_ww = sp0_male_ad_wc = sp0_male_ad_wr = sp0_male_ad_cc = sp0_male_ad_cr = sp0_male_ad_rr = sp0_female_ad_ww = sp0_female_ad_wc = sp0_female_ad_wr = sp0_female_ad_cc = sp0_female_ad_cr = sp0_female_ad_rr = 0.0
      sp0_male_ad_ww = 20.0
      sp0_female_ad_ww = 30.0
      ic_species0 = [male_juv_ww, male_juv_wc, male_juv_wr, male_juv_cc, male_juv_cr, male_juv_rr, female_juv_ww, female_juv_wc, female_juv_wr, female_juv_cc, female_juv_cr, female_juv_rr, sp0_male_ad_ww, sp0_male_ad_wc, sp0_male_ad_wr, sp0_male_ad_cc, sp0_male_ad_cr, sp0_male_ad_rr, sp0_female_ad_ww, sp0_female_ad_wc, sp0_female_ad_wr, sp0_female_ad_cc, sp0_female_ad_cr, sp0_female_ad_rr]

      # species1 adults:
      sp1_male_ad_ww = sp1_male_ad_wc = sp1_male_ad_wr = sp1_male_ad_cc = sp1_male_ad_cr = sp1_male_ad_rr = sp1_female_ad_ww = sp1_female_ad_wc = sp1_female_ad_wr = sp1_female_ad_cc = sp1_female_ad_cr = sp1_female_ad_rr = 0.0
      ic_species1 = [male_juv_ww, male_juv_wc, male_juv_wr, male_juv_cc, male_juv_cr, male_juv_rr, female_juv_ww, female_juv_wc, female_juv_wr, female_juv_cc, female_juv_cr, female_juv_rr, sp1_male_ad_ww, sp1_male_ad_wc, sp1_male_ad_wr, sp1_male_ad_cc, sp1_male_ad_cr, sp1_male_ad_rr, sp1_female_ad_ww, sp1_female_ad_wc, sp1_female_ad_wr, sp1_female_ad_cc, sp1_female_ad_cr, sp1_female_ad_rr]

      # build initial condition in the standard form
      initial_condition = []
      for i in range(24):
         initial_condition.append(ic_species0[i])
         initial_condition.append(ic_species1[i])
      initial_condition += [carrying_cap, carrying_cap]
      pap = array.array('f', initial_condition)

      # predict the answer
      prediction = predict_answer([carrying_cap, carrying_cap], 0.0, fecundity, mu_adult, mu_juv, aging, dt, pap)

      # evolve according to the cython code
      self.c.evolve(dt, pap)

      # check answer
      self.assertTrue(arrayfuzzyequal(prediction, pap[:48], ep_tol))


   def testEvolveSpecies_bothWW(self):
      dt = 3.0
      mu_adult = 0.125
      self.c.setMuAdult(mu_adult)
      mu_juv = 0.25
      self.c.setMuLarvae(mu_juv)
      fecundity = 8
      self.c.setFecundity(fecundity)
      carrying_cap = 60.0 # because we start with zero juveniles, this is actually irrelevant
      aging = 0.0 # because we start with zero juveniles, this is actually irrelevant
      matingcomp = 0.0
      self.c.setMatingComponent(0, 1, matingcomp)
      self.c.setMatingComponent(1, 0, matingcomp)

      # start with zero juveniles
      male_juv_ww = male_juv_wc = male_juv_wr = male_juv_cc = male_juv_cr = male_juv_rr = female_juv_ww = female_juv_wc = female_juv_wr = female_juv_cc = female_juv_cr = female_juv_rr = 0.0

      # species0 adults:
      sp0_male_ad_ww = sp0_male_ad_wc = sp0_male_ad_wr = sp0_male_ad_cc = sp0_male_ad_cr = sp0_male_ad_rr = sp0_female_ad_ww = sp0_female_ad_wc = sp0_female_ad_wr = sp0_female_ad_cc = sp0_female_ad_cr = sp0_female_ad_rr = 0.0
      sp0_male_ad_ww = 20.0
      sp0_female_ad_ww = 30.0
      ic_species0 = [male_juv_ww, male_juv_wc, male_juv_wr, male_juv_cc, male_juv_cr, male_juv_rr, female_juv_ww, female_juv_wc, female_juv_wr, female_juv_cc, female_juv_cr, female_juv_rr, sp0_male_ad_ww, sp0_male_ad_wc, sp0_male_ad_wr, sp0_male_ad_cc, sp0_male_ad_cr, sp0_male_ad_rr, sp0_female_ad_ww, sp0_female_ad_wc, sp0_female_ad_wr, sp0_female_ad_cc, sp0_female_ad_cr, sp0_female_ad_rr]

      # species1 adults:
      sp1_male_ad_ww = sp1_male_ad_wc = sp1_male_ad_wr = sp1_male_ad_cc = sp1_male_ad_cr = sp1_male_ad_rr = sp1_female_ad_ww = sp1_female_ad_wc = sp1_female_ad_wr = sp1_female_ad_cc = sp1_female_ad_cr = sp1_female_ad_rr = 0.0
      sp1_male_ad_ww = 10.0
      sp1_female_ad_ww = 10.0
      ic_species1 = [male_juv_ww, male_juv_wc, male_juv_wr, male_juv_cc, male_juv_cr, male_juv_rr, female_juv_ww, female_juv_wc, female_juv_wr, female_juv_cc, female_juv_cr, female_juv_rr, sp1_male_ad_ww, sp1_male_ad_wc, sp1_male_ad_wr, sp1_male_ad_cc, sp1_male_ad_cr, sp1_male_ad_rr, sp1_female_ad_ww, sp1_female_ad_wc, sp1_female_ad_wr, sp1_female_ad_cc, sp1_female_ad_cr, sp1_female_ad_rr]

      # build initial condition in the standard form
      initial_condition = []
      for i in range(24):
         initial_condition.append(ic_species0[i])
         initial_condition.append(ic_species1[i])
      initial_condition += [carrying_cap, carrying_cap]
      pap = array.array('f', initial_condition)

      # predict the answer
      prediction = predict_answer([carrying_cap, carrying_cap], matingcomp, fecundity, mu_adult, mu_juv, aging, dt, pap)

      # evolve according to the cython code
      self.c.evolve(dt, pap)

      # check answer
      self.assertTrue(arrayfuzzyequal(prediction, pap[:48], ep_tol))

   def testEvolveSpecies_bothWW_matingcomp(self):
      dt = 3.0
      mu_adult = 0.125
      self.c.setMuAdult(mu_adult)
      mu_juv = 0.25
      self.c.setMuLarvae(mu_juv)
      fecundity = 8
      self.c.setFecundity(fecundity)
      carrying_cap = 60.0 # because we start with zero juveniles, this is actually irrelevant
      aging = 0.0 # because we start with zero juveniles, this is actually irrelevant
      matingcomp = 0.5
      self.c.setMatingComponent(0, 1, matingcomp)
      self.c.setMatingComponent(1, 0, matingcomp)

      # start with zero juveniles
      male_juv_ww = male_juv_wc = male_juv_wr = male_juv_cc = male_juv_cr = male_juv_rr = female_juv_ww = female_juv_wc = female_juv_wr = female_juv_cc = female_juv_cr = female_juv_rr = 0.0

      # species0 adults:
      sp0_male_ad_ww = sp0_male_ad_wc = sp0_male_ad_wr = sp0_male_ad_cc = sp0_male_ad_cr = sp0_male_ad_rr = sp0_female_ad_ww = sp0_female_ad_wc = sp0_female_ad_wr = sp0_female_ad_cc = sp0_female_ad_cr = sp0_female_ad_rr = 0.0
      sp0_male_ad_ww = 20.0
      sp0_female_ad_ww = 30.0
      ic_species0 = [male_juv_ww, male_juv_wc, male_juv_wr, male_juv_cc, male_juv_cr, male_juv_rr, female_juv_ww, female_juv_wc, female_juv_wr, female_juv_cc, female_juv_cr, female_juv_rr, sp0_male_ad_ww, sp0_male_ad_wc, sp0_male_ad_wr, sp0_male_ad_cc, sp0_male_ad_cr, sp0_male_ad_rr, sp0_female_ad_ww, sp0_female_ad_wc, sp0_female_ad_wr, sp0_female_ad_cc, sp0_female_ad_cr, sp0_female_ad_rr]

      # species1 adults:
      sp1_male_ad_ww = sp1_male_ad_wc = sp1_male_ad_wr = sp1_male_ad_cc = sp1_male_ad_cr = sp1_male_ad_rr = sp1_female_ad_ww = sp1_female_ad_wc = sp1_female_ad_wr = sp1_female_ad_cc = sp1_female_ad_cr = sp1_female_ad_rr = 0.0
      sp1_male_ad_ww = 10.0
      sp1_female_ad_ww = 10.0
      ic_species1 = [male_juv_ww, male_juv_wc, male_juv_wr, male_juv_cc, male_juv_cr, male_juv_rr, female_juv_ww, female_juv_wc, female_juv_wr, female_juv_cc, female_juv_cr, female_juv_rr, sp1_male_ad_ww, sp1_male_ad_wc, sp1_male_ad_wr, sp1_male_ad_cc, sp1_male_ad_cr, sp1_male_ad_rr, sp1_female_ad_ww, sp1_female_ad_wc, sp1_female_ad_wr, sp1_female_ad_cc, sp1_female_ad_cr, sp1_female_ad_rr]

      # build initial condition in the standard form
      initial_condition = []
      for i in range(24):
         initial_condition.append(ic_species0[i])
         initial_condition.append(ic_species1[i])
      initial_condition += [carrying_cap, carrying_cap]
      pap = array.array('f', initial_condition)

      # predict the answer
      prediction = predict_answer([carrying_cap, carrying_cap], matingcomp, fecundity, mu_adult, mu_juv, aging, dt, pap)

      # evolve according to the cython code
      self.c.evolve(dt, pap)

      # check answer
      self.assertTrue(arrayfuzzyequal(prediction, pap[:48], ep_tol))

   def testEvolveSpecies_bothWC(self):
      dt = 3.0
      mu_adult = 0.125
      self.c.setMuAdult(mu_adult)
      mu_juv = 0.25
      self.c.setMuLarvae(mu_juv)
      fecundity = 8
      self.c.setFecundity(fecundity)
      carrying_cap = 60.0 # because we start with zero juveniles, this is actually irrelevant
      aging = 0.0 # because we start with zero juveniles, this is actually irrelevant
      matingcomp = 0.0
      self.c.setMatingComponent(0, 1, matingcomp)
      self.c.setMatingComponent(1, 0, matingcomp)

      # start with zero juveniles
      male_juv_ww = male_juv_wc = male_juv_wr = male_juv_cc = male_juv_cr = male_juv_rr = female_juv_ww = female_juv_wc = female_juv_wr = female_juv_cc = female_juv_cr = female_juv_rr = 0.0

      # species0 adults:
      sp0_male_ad_ww = sp0_male_ad_wc = sp0_male_ad_wr = sp0_male_ad_cc = sp0_male_ad_cr = sp0_male_ad_rr = sp0_female_ad_ww = sp0_female_ad_wc = sp0_female_ad_wr = sp0_female_ad_cc = sp0_female_ad_cr = sp0_female_ad_rr = 0.0
      sp0_male_ad_wc = 20.0
      sp0_female_ad_wc = 30.0
      ic_species0 = [male_juv_ww, male_juv_wc, male_juv_wr, male_juv_cc, male_juv_cr, male_juv_rr, female_juv_ww, female_juv_wc, female_juv_wr, female_juv_cc, female_juv_cr, female_juv_rr, sp0_male_ad_ww, sp0_male_ad_wc, sp0_male_ad_wr, sp0_male_ad_cc, sp0_male_ad_cr, sp0_male_ad_rr, sp0_female_ad_ww, sp0_female_ad_wc, sp0_female_ad_wr, sp0_female_ad_cc, sp0_female_ad_cr, sp0_female_ad_rr]

      # species1 adults:
      sp1_male_ad_ww = sp1_male_ad_wc = sp1_male_ad_wr = sp1_male_ad_cc = sp1_male_ad_cr = sp1_male_ad_rr = sp1_female_ad_ww = sp1_female_ad_wc = sp1_female_ad_wr = sp1_female_ad_cc = sp1_female_ad_cr = sp1_female_ad_rr = 0.0
      sp1_male_ad_wc = 10.0
      sp1_female_ad_wc = 10.0
      ic_species1 = [male_juv_ww, male_juv_wc, male_juv_wr, male_juv_cc, male_juv_cr, male_juv_rr, female_juv_ww, female_juv_wc, female_juv_wr, female_juv_cc, female_juv_cr, female_juv_rr, sp1_male_ad_ww, sp1_male_ad_wc, sp1_male_ad_wr, sp1_male_ad_cc, sp1_male_ad_cr, sp1_male_ad_rr, sp1_female_ad_ww, sp1_female_ad_wc, sp1_female_ad_wr, sp1_female_ad_cc, sp1_female_ad_cr, sp1_female_ad_rr]

      # build initial condition in the standard form
      initial_condition = []
      for i in range(24):
         initial_condition.append(ic_species0[i])
         initial_condition.append(ic_species1[i])
      initial_condition += [carrying_cap, carrying_cap]
      pap = array.array('f', initial_condition)

      # predict the answer
      prediction = predict_answer([carrying_cap, carrying_cap], matingcomp, fecundity, mu_adult, mu_juv, aging, dt, pap)

      # evolve according to the cython code
      self.c.evolve(dt, pap)

      # check answer
      self.assertTrue(arrayfuzzyequal(prediction, pap[:48], ep_tol))

   def testEvolveSpecies_bothWC_matingcomp(self):
      dt = 3.0
      mu_adult = 0.125
      self.c.setMuAdult(mu_adult)
      mu_juv = 0.25
      self.c.setMuLarvae(mu_juv)
      fecundity = 8
      self.c.setFecundity(fecundity)
      carrying_cap = 60.0 # because we start with zero juveniles, this is actually irrelevant
      aging = 0.0 # because we start with zero juveniles, this is actually irrelevant
      matingcomp = 0.5
      self.c.setMatingComponent(0, 1, matingcomp)
      self.c.setMatingComponent(1, 0, matingcomp)

      # start with zero juveniles
      male_juv_ww = male_juv_wc = male_juv_wr = male_juv_cc = male_juv_cr = male_juv_rr = female_juv_ww = female_juv_wc = female_juv_wr = female_juv_cc = female_juv_cr = female_juv_rr = 0.0

      # species0 adults:
      sp0_male_ad_ww = sp0_male_ad_wc = sp0_male_ad_wr = sp0_male_ad_cc = sp0_male_ad_cr = sp0_male_ad_rr = sp0_female_ad_ww = sp0_female_ad_wc = sp0_female_ad_wr = sp0_female_ad_cc = sp0_female_ad_cr = sp0_female_ad_rr = 0.0
      sp0_male_ad_wc = 20.0
      sp0_female_ad_wc = 30.0
      ic_species0 = [male_juv_ww, male_juv_wc, male_juv_wr, male_juv_cc, male_juv_cr, male_juv_rr, female_juv_ww, female_juv_wc, female_juv_wr, female_juv_cc, female_juv_cr, female_juv_rr, sp0_male_ad_ww, sp0_male_ad_wc, sp0_male_ad_wr, sp0_male_ad_cc, sp0_male_ad_cr, sp0_male_ad_rr, sp0_female_ad_ww, sp0_female_ad_wc, sp0_female_ad_wr, sp0_female_ad_cc, sp0_female_ad_cr, sp0_female_ad_rr]

      # species1 adults:
      sp1_male_ad_ww = sp1_male_ad_wc = sp1_male_ad_wr = sp1_male_ad_cc = sp1_male_ad_cr = sp1_male_ad_rr = sp1_female_ad_ww = sp1_female_ad_wc = sp1_female_ad_wr = sp1_female_ad_cc = sp1_female_ad_cr = sp1_female_ad_rr = 0.0
      sp1_male_ad_wc = 10.0
      sp1_female_ad_wc = 10.0
      ic_species1 = [male_juv_ww, male_juv_wc, male_juv_wr, male_juv_cc, male_juv_cr, male_juv_rr, female_juv_ww, female_juv_wc, female_juv_wr, female_juv_cc, female_juv_cr, female_juv_rr, sp1_male_ad_ww, sp1_male_ad_wc, sp1_male_ad_wr, sp1_male_ad_cc, sp1_male_ad_cr, sp1_male_ad_rr, sp1_female_ad_ww, sp1_female_ad_wc, sp1_female_ad_wr, sp1_female_ad_cc, sp1_female_ad_cr, sp1_female_ad_rr]

      # build initial condition in the standard form
      initial_condition = []
      for i in range(24):
         initial_condition.append(ic_species0[i])
         initial_condition.append(ic_species1[i])
      initial_condition += [carrying_cap, carrying_cap]
      pap = array.array('f', initial_condition)

      # predict the answer
      prediction = predict_answer([carrying_cap, carrying_cap], matingcomp, fecundity, mu_adult, mu_juv, aging, dt, pap)

      # evolve according to the cython code
      self.c.evolve(dt, pap)

      # check answer
      self.assertTrue(arrayfuzzyequal(prediction, pap[:48], ep_tol))

   def testEvolveSpecies_bothWC_matingcomp_somejuvs(self):
      dt = 3.0E-1
      self.c.setMinimumDt(3.0E-1)
      mu_adult = 0.125
      self.c.setMuAdult(mu_adult)
      mu_juv = 0.25
      self.c.setMuLarvae(mu_juv)
      fecundity = 8
      self.c.setFecundity(fecundity)
      carrying_cap = [90.0, 60.0]
      aging = 0.1 # default of constructor in mosquito23
      matingcomp = 0.5
      self.c.setMatingComponent(0, 1, matingcomp)
      self.c.setMatingComponent(1, 0, matingcomp)

      # species0: start with nonzero juveniles
      sp0_male_juv_ww = 0.2
      sp0_male_juv_wc = 0.3
      sp0_male_juv_wr = 0.4
      sp0_male_juv_cc = 0.5
      sp0_male_juv_cr = 0.6
      sp0_male_juv_rr = 0.7
      sp0_female_juv_ww = 0.8
      sp0_female_juv_wc = 0.0
      sp0_female_juv_wr = 0.15
      sp0_female_juv_cc = 0.5
      sp0_female_juv_cr = 0.15
      sp0_female_juv_rr = 0.0
      
      # species0 adults:
      sp0_male_ad_ww = sp0_male_ad_wc = sp0_male_ad_wr = sp0_male_ad_cc = sp0_male_ad_cr = sp0_male_ad_rr = sp0_female_ad_ww = sp0_female_ad_wc = sp0_female_ad_wr = sp0_female_ad_cc = sp0_female_ad_cr = sp0_female_ad_rr = 0.0
      sp0_male_ad_wc = 20.0
      sp0_female_ad_wc = 30.0
      ic_species0 = [sp0_male_juv_ww, sp0_male_juv_wc, sp0_male_juv_wr, sp0_male_juv_cc, sp0_male_juv_cr, sp0_male_juv_rr, sp0_female_juv_ww, sp0_female_juv_wc, sp0_female_juv_wr, sp0_female_juv_cc, sp0_female_juv_cr, sp0_female_juv_rr, sp0_male_ad_ww, sp0_male_ad_wc, sp0_male_ad_wr, sp0_male_ad_cc, sp0_male_ad_cr, sp0_male_ad_rr, sp0_female_ad_ww, sp0_female_ad_wc, sp0_female_ad_wr, sp0_female_ad_cc, sp0_female_ad_cr, sp0_female_ad_rr]

      # species1: start with nonzero juveniles
      sp1_male_juv_ww = 0.9
      sp1_male_juv_wc = 0.7
      sp1_male_juv_wr = 0.5
      sp1_male_juv_cc = 0.3
      sp1_male_juv_cr = 0.1
      sp1_male_juv_rr = 0.3
      sp1_female_juv_ww = 0.3
      sp1_female_juv_wc = 0.4
      sp1_female_juv_wr = 0.6
      sp1_female_juv_cc = 0.8
      sp1_female_juv_cr = 0.3
      sp1_female_juv_rr = 0.2

      # species1 adults:
      sp1_male_ad_ww = sp1_male_ad_wc = sp1_male_ad_wr = sp1_male_ad_cc = sp1_male_ad_cr = sp1_male_ad_rr = sp1_female_ad_ww = sp1_female_ad_wc = sp1_female_ad_wr = sp1_female_ad_cc = sp1_female_ad_cr = sp1_female_ad_rr = 0.0
      sp1_male_ad_wc = 10.0
      sp1_female_ad_wc = 10.0
      ic_species1 = [sp1_male_juv_ww, sp1_male_juv_wc, sp1_male_juv_wr, sp1_male_juv_cc, sp1_male_juv_cr, sp1_male_juv_rr, sp1_female_juv_ww, sp1_female_juv_wc, sp1_female_juv_wr, sp1_female_juv_cc, sp1_female_juv_cr, sp1_female_juv_rr, sp1_male_ad_ww, sp1_male_ad_wc, sp1_male_ad_wr, sp1_male_ad_cc, sp1_male_ad_cr, sp1_male_ad_rr, sp1_female_ad_ww, sp1_female_ad_wc, sp1_female_ad_wr, sp1_female_ad_cc, sp1_female_ad_cr, sp1_female_ad_rr]

      # build initial condition in the standard form
      initial_condition = []
      for i in range(24):
         initial_condition.append(ic_species0[i])
         initial_condition.append(ic_species1[i])
      initial_condition += carrying_cap
      pap = array.array('f', initial_condition)

      # predict the answer
      prediction = predict_answer(carrying_cap, matingcomp, fecundity, mu_adult, mu_juv, aging, dt, pap)

      # evolve according to the cython code
      self.c.evolve(dt, pap)

      # check answer
      self.assertTrue(arrayfuzzyequal(prediction, pap[:48], ep_tol))

   def testEvolveSpecies_allGeno(self):
      dt = 3.0
      mu_adult = 0.125
      self.c.setMuAdult(mu_adult)
      mu_juv = 0.25
      self.c.setMuLarvae(mu_juv)
      fecundity = 8
      self.c.setFecundity(fecundity)
      carrying_cap = 60.0 # because we start with zero juveniles, this is actually irrelevant
      aging = 0.0 # because we start with zero juveniles, this is actually irrelevant
      matingcomp = 0.0
      self.c.setMatingComponent(0, 1, matingcomp)
      self.c.setMatingComponent(1, 0, matingcomp)

      # start with zero juveniles
      male_juv_ww = male_juv_wc = male_juv_wr = male_juv_cc = male_juv_cr = male_juv_rr = female_juv_ww = female_juv_wc = female_juv_wr = female_juv_cc = female_juv_cr = female_juv_rr = 0.0

      # species0 adults:
      sp0_male_ad_ww = 10.0
      sp0_male_ad_wc = 12.0
      sp0_male_ad_wr = 16.0
      sp0_male_ad_cc = 20.0
      sp0_male_ad_cr = 4.0
      sp0_male_ad_rr = 9.0
      sp0_female_ad_ww = 8.0
      sp0_female_ad_wc = 12.0
      sp0_female_ad_wr = 4.0
      sp0_female_ad_cc = 6.0
      sp0_female_ad_cr = 13.0
      sp0_female_ad_rr = 2.0
      ic_species0 = [male_juv_ww, male_juv_wc, male_juv_wr, male_juv_cc, male_juv_cr, male_juv_rr, female_juv_ww, female_juv_wc, female_juv_wr, female_juv_cc, female_juv_cr, female_juv_rr, sp0_male_ad_ww, sp0_male_ad_wc, sp0_male_ad_wr, sp0_male_ad_cc, sp0_male_ad_cr, sp0_male_ad_rr, sp0_female_ad_ww, sp0_female_ad_wc, sp0_female_ad_wr, sp0_female_ad_cc, sp0_female_ad_cr, sp0_female_ad_rr]

      # species1 adults:
      sp1_male_ad_ww = 1.0
      sp1_male_ad_wc = 1.2
      sp1_male_ad_wr = 1.6
      sp1_male_ad_cc = 22.0
      sp1_male_ad_cr = 4.4
      sp1_male_ad_rr = 9.9
      sp1_female_ad_ww = 13.0
      sp1_female_ad_wc = 12.0
      sp1_female_ad_wr = 14.0
      sp1_female_ad_cc = 16.0
      sp1_female_ad_cr = 1.3
      sp1_female_ad_rr = 0.2
      ic_species1 = [male_juv_ww, male_juv_wc, male_juv_wr, male_juv_cc, male_juv_cr, male_juv_rr, female_juv_ww, female_juv_wc, female_juv_wr, female_juv_cc, female_juv_cr, female_juv_rr, sp1_male_ad_ww, sp1_male_ad_wc, sp1_male_ad_wr, sp1_male_ad_cc, sp1_male_ad_cr, sp1_male_ad_rr, sp1_female_ad_ww, sp1_female_ad_wc, sp1_female_ad_wr, sp1_female_ad_cc, sp1_female_ad_cr, sp1_female_ad_rr]

      # build initial condition in the standard form
      initial_condition = []
      for i in range(24):
         initial_condition.append(ic_species0[i])
         initial_condition.append(ic_species1[i])
      initial_condition += [carrying_cap, carrying_cap]
      pap = array.array('f', initial_condition)

      # predict the answer
      prediction = predict_answer([carrying_cap, carrying_cap], matingcomp, fecundity, mu_adult, mu_juv, aging, dt, pap)

      # evolve according to the cython code
      self.c.evolve(dt, pap)

      # check answer
      self.assertTrue(arrayfuzzyequal(prediction, pap[:48], ep_tol))

   def testEvolveSpecies_allGeno_matingcomp(self):
      dt = 3.0
      mu_adult = 0.125
      self.c.setMuAdult(mu_adult)
      mu_juv = 0.25
      self.c.setMuLarvae(mu_juv)
      fecundity = 8
      self.c.setFecundity(fecundity)
      carrying_cap = 60.0 # because we start with zero juveniles, this is actually irrelevant
      aging = 0.0 # because we start with zero juveniles, this is actually irrelevant
      matingcomp = 0.5
      self.c.setMatingComponent(0, 1, matingcomp)
      self.c.setMatingComponent(1, 0, matingcomp)

      # start with zero juveniles
      male_juv_ww = male_juv_wc = male_juv_wr = male_juv_cc = male_juv_cr = male_juv_rr = female_juv_ww = female_juv_wc = female_juv_wr = female_juv_cc = female_juv_cr = female_juv_rr = 0.0

      # species0 adults:
      sp0_male_ad_ww = 10.0
      sp0_male_ad_wc = 12.0
      sp0_male_ad_wr = 16.0
      sp0_male_ad_cc = 20.0
      sp0_male_ad_cr = 4.0
      sp0_male_ad_rr = 9.0
      sp0_female_ad_ww = 8.0
      sp0_female_ad_wc = 12.0
      sp0_female_ad_wr = 4.0
      sp0_female_ad_cc = 6.0
      sp0_female_ad_cr = 13.0
      sp0_female_ad_rr = 2.0
      ic_species0 = [male_juv_ww, male_juv_wc, male_juv_wr, male_juv_cc, male_juv_cr, male_juv_rr, female_juv_ww, female_juv_wc, female_juv_wr, female_juv_cc, female_juv_cr, female_juv_rr, sp0_male_ad_ww, sp0_male_ad_wc, sp0_male_ad_wr, sp0_male_ad_cc, sp0_male_ad_cr, sp0_male_ad_rr, sp0_female_ad_ww, sp0_female_ad_wc, sp0_female_ad_wr, sp0_female_ad_cc, sp0_female_ad_cr, sp0_female_ad_rr]

      # species1 adults:
      sp1_male_ad_ww = 1.0
      sp1_male_ad_wc = 1.2
      sp1_male_ad_wr = 1.6
      sp1_male_ad_cc = 22.0
      sp1_male_ad_cr = 4.4
      sp1_male_ad_rr = 9.9
      sp1_female_ad_ww = 13.0
      sp1_female_ad_wc = 12.0
      sp1_female_ad_wr = 14.0
      sp1_female_ad_cc = 16.0
      sp1_female_ad_cr = 1.3
      sp1_female_ad_rr = 0.2
      ic_species1 = [male_juv_ww, male_juv_wc, male_juv_wr, male_juv_cc, male_juv_cr, male_juv_rr, female_juv_ww, female_juv_wc, female_juv_wr, female_juv_cc, female_juv_cr, female_juv_rr, sp1_male_ad_ww, sp1_male_ad_wc, sp1_male_ad_wr, sp1_male_ad_cc, sp1_male_ad_cr, sp1_male_ad_rr, sp1_female_ad_ww, sp1_female_ad_wc, sp1_female_ad_wr, sp1_female_ad_cc, sp1_female_ad_cr, sp1_female_ad_rr]

      # build initial condition in the standard form
      initial_condition = []
      for i in range(24):
         initial_condition.append(ic_species0[i])
         initial_condition.append(ic_species1[i])
      initial_condition += [carrying_cap, carrying_cap]
      pap = array.array('f', initial_condition)

      # predict the answer
      prediction = predict_answer([carrying_cap, carrying_cap], matingcomp, fecundity, mu_adult, mu_juv, aging, dt, pap)

      # evolve according to the cython code
      self.c.evolve(dt, pap)

      # check answer
      self.assertTrue(arrayfuzzyequal(prediction, pap[:48], ep_tol))

   def testEvolveSpecies_allGeno_matingcomp_somejuvs(self):
      dt = 3.0E-1
      self.c.setMinimumDt(3.0E-1)
      mu_adult = 0.125
      self.c.setMuAdult(mu_adult)
      mu_juv = 0.25
      self.c.setMuLarvae(mu_juv)
      fecundity = 8
      self.c.setFecundity(fecundity)
      carrying_cap = [90.0, 60.0]
      aging = 0.1 # default of constructor in mosquito23
      matingcomp = 0.5
      self.c.setMatingComponent(0, 1, matingcomp)
      self.c.setMatingComponent(1, 0, matingcomp)

      # species0: start with nonzero juveniles
      sp0_male_juv_ww = 0.2
      sp0_male_juv_wc = 0.3
      sp0_male_juv_wr = 0.4
      sp0_male_juv_cc = 0.5
      sp0_male_juv_cr = 0.6
      sp0_male_juv_rr = 0.7
      sp0_female_juv_ww = 0.8
      sp0_female_juv_wc = 0.0
      sp0_female_juv_wr = 0.15
      sp0_female_juv_cc = 0.5
      sp0_female_juv_cr = 0.15
      sp0_female_juv_rr = 0.0

      # species0 adults:
      sp0_male_ad_ww = 10.0
      sp0_male_ad_wc = 12.0
      sp0_male_ad_wr = 16.0
      sp0_male_ad_cc = 20.0
      sp0_male_ad_cr = 4.0
      sp0_male_ad_rr = 9.0
      sp0_female_ad_ww = 8.0
      sp0_female_ad_wc = 12.0
      sp0_female_ad_wr = 4.0
      sp0_female_ad_cc = 6.0
      sp0_female_ad_cr = 13.0
      sp0_female_ad_rr = 2.0
      ic_species0 = [sp0_male_juv_ww, sp0_male_juv_wc, sp0_male_juv_wr, sp0_male_juv_cc, sp0_male_juv_cr, sp0_male_juv_rr, sp0_female_juv_ww, sp0_female_juv_wc, sp0_female_juv_wr, sp0_female_juv_cc, sp0_female_juv_cr, sp0_female_juv_rr, sp0_male_ad_ww, sp0_male_ad_wc, sp0_male_ad_wr, sp0_male_ad_cc, sp0_male_ad_cr, sp0_male_ad_rr, sp0_female_ad_ww, sp0_female_ad_wc, sp0_female_ad_wr, sp0_female_ad_cc, sp0_female_ad_cr, sp0_female_ad_rr]

      # species1: start with nonzero juveniles
      sp1_male_juv_ww = 0.9
      sp1_male_juv_wc = 0.7
      sp1_male_juv_wr = 0.5
      sp1_male_juv_cc = 0.3
      sp1_male_juv_cr = 0.1
      sp1_male_juv_rr = 0.3
      sp1_female_juv_ww = 0.3
      sp1_female_juv_wc = 0.4
      sp1_female_juv_wr = 0.6
      sp1_female_juv_cc = 0.8
      sp1_female_juv_cr = 0.3
      sp1_female_juv_rr = 0.2

      # species1 adults:
      sp1_male_ad_ww = 1.0
      sp1_male_ad_wc = 1.2
      sp1_male_ad_wr = 1.6
      sp1_male_ad_cc = 22.0
      sp1_male_ad_cr = 4.4
      sp1_male_ad_rr = 9.9
      sp1_female_ad_ww = 13.0
      sp1_female_ad_wc = 12.0
      sp1_female_ad_wr = 14.0
      sp1_female_ad_cc = 16.0
      sp1_female_ad_cr = 1.3
      sp1_female_ad_rr = 0.2
      ic_species1 = [sp1_male_juv_ww, sp1_male_juv_wc, sp1_male_juv_wr, sp1_male_juv_cc, sp1_male_juv_cr, sp1_male_juv_rr, sp1_female_juv_ww, sp1_female_juv_wc, sp1_female_juv_wr, sp1_female_juv_cc, sp1_female_juv_cr, sp1_female_juv_rr, sp1_male_ad_ww, sp1_male_ad_wc, sp1_male_ad_wr, sp1_male_ad_cc, sp1_male_ad_cr, sp1_male_ad_rr, sp1_female_ad_ww, sp1_female_ad_wc, sp1_female_ad_wr, sp1_female_ad_cc, sp1_female_ad_cr, sp1_female_ad_rr]

      # build initial condition in the standard form
      initial_condition = []
      for i in range(24):
         initial_condition.append(ic_species0[i])
         initial_condition.append(ic_species1[i])
      initial_condition += carrying_cap
      pap = array.array('f', initial_condition)

      # predict the answer
      prediction = predict_answer(carrying_cap, matingcomp, fecundity, mu_adult, mu_juv, aging, dt, pap)

      # evolve according to the cython code
      self.c.evolve(dt, pap)

      # check answer
      self.assertTrue(arrayfuzzyequal(prediction, pap[:48], ep_tol))

   def testGetSetSmallValue(self):
      self.assertEqual(self.c.getSmallValue(), 0.0)
      self.c.setSmallValue(1.0)
      self.assertEqual(self.c.getSmallValue(), 1.0)

if __name__ == '__main__':
   unittest.main()

