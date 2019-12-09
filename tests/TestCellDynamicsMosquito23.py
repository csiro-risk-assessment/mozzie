import os
import sys
import unittest
import array
import random

# so we can find our ../code no matter how we are called
findbin = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(findbin + "/../code")

from cellDynamics import CellDynamicsMosquito23

def arrayequal(a, b):
   return all([a[i] == b[i] for i in range(0, len(a))])

def arrayfuzzyequal(a, b, eps):
   return all([(a[i] > b[i] - eps and a[i] < b[i] + eps) for i in range(0, len(a))])

class TestCellDynamicsMosquito23(unittest.TestCase):

   def setUp(self):
      self.c = CellDynamicsMosquito23()

   def testGetNumberOfPopulations(self):
      self.assertEqual(self.c.getNumberOfPopulations(), 2 * 3 * 2 * 1)
      self.c.setNumAges(4)
      self.c.setNumSpecies(5)
      self.assertEqual(self.c.getNumberOfPopulations(), 2 * 3 * 4 * 5)

   def testGetNumberOfParameters(self):
      self.assertEqual(self.c.getNumberOfParameters(), 1)

   def testGetNumberOfDiffusingPopulations(self):
      self.assertEqual(self.c.getNumberOfDiffusingPopulations(), 6)
      self.c.setNumAges(4)
      self.c.setNumSpecies(3)
      self.assertEqual(self.c.getNumberOfDiffusingPopulations(), 6 * 3)

   def testGetDiffusingIndices(self):
      self.assertTrue(arrayequal(self.c.getDiffusingIndices(), [6, 7, 8, 9, 10, 11]))
      self.c.setNumAges(4)
      self.c.setNumSpecies(2)
      self.assertTrue(arrayequal(self.c.getDiffusingIndices(), list(range(36, 48))))

   def testGetNumberOfAdvectingPopulations(self):
      self.assertEqual(self.c.getNumberOfAdvectingPopulations(), 6)
      self.c.setNumAges(4)
      self.c.setNumSpecies(5)
      self.assertEqual(self.c.getNumberOfAdvectingPopulations(), 6 * 5)

   def testGetAdvectingIndices(self):
      self.assertTrue(arrayequal(self.c.getAdvectingIndices(), [6, 7, 8, 9, 10, 11]))
      self.c.setNumAges(5)
      self.c.setNumSpecies(2)
      self.assertTrue(arrayequal(self.c.getAdvectingIndices(), list(range(48, 60))))

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
      self.c.setNumSpecies(3)
      self.c.setAlphaComponent(0, 0, 1.23)
      self.c.setAlphaComponent(0, 1, -1.23)
      self.c.setAlphaComponent(0, 2, 3.33)
      self.c.setAlphaComponent(1, 0, 1.4)
      self.c.setAlphaComponent(1, 2, -0.5)
      self.c.setAlphaComponent(2, 0, 0.8)
      self.assertTrue(arrayfuzzyequal([self.c.getAlphaComponentFromPython(0, 0)], [1.23], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getAlphaComponentFromPython(0, 1)], [-1.23], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getAlphaComponentFromPython(0, 2)], [3.33], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getAlphaComponentFromPython(1, 0)], [1.4], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getAlphaComponentFromPython(1, 1)], [1], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getAlphaComponentFromPython(1, 2)], [-0.5], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getAlphaComponentFromPython(2, 0)], [0.8], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getAlphaComponentFromPython(2, 1)], [0.0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getAlphaComponentFromPython(2, 2)], [1.0], 1E-6))
      with self.assertRaises(ValueError) as the_err:
         self.c.setAlphaComponent(0, 3, 0)
      self.assertEqual(str(the_err.exception), "sp0 0 and sp1 3 must be less than the number of species, 3")
      with self.assertRaises(ValueError) as the_err:
         self.c.setAlphaComponent(3, 1, 0)
      self.assertEqual(str(the_err.exception), "sp0 3 and sp1 1 must be less than the number of species, 3")

   def testSetGetHybridisationRate(self):
      self.assertEqual(self.c.getHybridisationRateFromPython(0, 0, 0), 1.0)
      self.c.setNumSpecies(3)
      for i in range(3):
         for j in range(3):
            for k in range(3):
               self.assertTrue(arrayfuzzyequal([self.c.getHybridisationRateFromPython(i, j, k)], [1 if (i == j and j == k) else 0], 1E-8))
      self.c.setHybridisationRate(0, 0, 0, 1.23)
      self.c.setHybridisationRate(0, 1, 1, -1.23)
      self.c.setHybridisationRate(0, 2, 2, 3.33)
      self.c.setHybridisationRate(1, 0, 0, 1.4)
      self.c.setHybridisationRate(1, 2, 1, -0.5)
      self.c.setHybridisationRate(2, 0, 2, 0.8)
      self.assertTrue(arrayfuzzyequal([self.c.getHybridisationRateFromPython(0, 0, 0)], [1.23], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getHybridisationRateFromPython(0, 1, 1)], [-1.23], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getHybridisationRateFromPython(0, 2, 2)], [3.33], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getHybridisationRateFromPython(1, 0, 0)], [1.4], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getHybridisationRateFromPython(1, 2, 1)], [-0.5], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getHybridisationRateFromPython(2, 0, 2)], [0.8], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getHybridisationRateFromPython(1, 1, 1)], [1.0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getHybridisationRateFromPython(2, 2, 2)], [1.0], 1E-6))
      with self.assertRaises(ValueError) as the_err:
         self.c.setHybridisationRate(3, 1, 0, 0)
      self.assertEqual(str(the_err.exception), "All species numbers, 3, 1, 0 must be less than the number of species, 3")
      with self.assertRaises(ValueError) as the_err:
         self.c.setHybridisationRate(0, 3, 0, 0)
      self.assertEqual(str(the_err.exception), "All species numbers, 0, 3, 0 must be less than the number of species, 3")
      with self.assertRaises(ValueError) as the_err:
         self.c.setHybridisationRate(1, 1, 3, 0)
      self.assertEqual(str(the_err.exception), "All species numbers, 1, 1, 3 must be less than the number of species, 3")

   def testSetGetMating(self):
      self.assertEqual(self.c.getMatingComponentFromPython(0, 0), 1.0)
      self.c.setNumSpecies(3)
      self.c.setMatingComponent(0, 0, 1.23)
      self.c.setMatingComponent(0, 1, -1.23)
      self.c.setMatingComponent(0, 2, 3.33)
      self.c.setMatingComponent(1, 0, 1.4)
      self.c.setMatingComponent(1, 2, -0.5)
      self.c.setMatingComponent(2, 0, 0.8)
      self.assertTrue(arrayfuzzyequal([self.c.getMatingComponentFromPython(0, 0)], [1.23], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getMatingComponentFromPython(0, 1)], [-1.23], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getMatingComponentFromPython(0, 2)], [3.33], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getMatingComponentFromPython(1, 0)], [1.4], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getMatingComponentFromPython(1, 1)], [1], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getMatingComponentFromPython(1, 2)], [-0.5], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getMatingComponentFromPython(2, 0)], [0.8], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getMatingComponentFromPython(2, 1)], [0.0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getMatingComponentFromPython(2, 2)], [1.0], 1E-6))
      with self.assertRaises(ValueError) as the_err:
         self.c.setMatingComponent(0, 3, 0)
      self.assertEqual(str(the_err.exception), "species_father 0 and species_mother 3 must be less than the number of species, 3")
      with self.assertRaises(ValueError) as the_err:
         self.c.setMatingComponent(3, 1, 0)
      self.assertEqual(str(the_err.exception), "species_father 3 and species_mother 1 must be less than the number of species, 3")

   def testInheritance(self):
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 0, 0)], [1], 1E-8))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 0, 1)], [0], 1E-8))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 0, 2)], [0], 1E-8))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 1, 0)], [0.5], 1E-8))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 1, 1)], [0.5], 1E-8))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 1, 2)], [0], 1E-8))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 2, 0)], [0], 1E-8))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 2, 1)], [1], 1E-8))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(0, 2, 2)], [0], 1E-8))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 0, 0)], [0.5], 1E-8))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 0, 1)], [0.5], 1E-8))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 0, 2)], [0], 1E-8))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 1, 0)], [0.25], 1E-8))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 1, 1)], [0.5], 1E-8))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 1, 2)], [0.25], 1E-8))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 2, 0)], [0], 1E-8))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 2, 1)], [0.5], 1E-8))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(1, 2, 2)], [0.5], 1E-8))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 0, 0)], [0], 1E-8))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 0, 1)], [1], 1E-8))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 0, 2)], [0], 1E-8))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 1, 0)], [0], 1E-8))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 1, 1)], [0.5], 1E-8))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 1, 2)], [0.5], 1E-8))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 2, 0)], [0], 1E-8))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 2, 1)], [0], 1E-8))
      self.assertTrue(arrayfuzzyequal([self.c.getInheritanceFromPython(2, 2, 2)], [1], 1E-8))
      with self.assertRaises(ValueError) as the_err:
         self.c.getInheritanceFromPython(3, 0, 0)
      self.assertEqual(str(the_err.exception), "All genotypes, 3, 0, 0 must be less than the number of genotypes, 3")
      with self.assertRaises(ValueError) as the_err:
         self.c.getInheritanceFromPython(0, 3, 0)
      self.assertEqual(str(the_err.exception), "All genotypes, 0, 3, 0 must be less than the number of genotypes, 3")
      with self.assertRaises(ValueError) as the_err:
         self.c.getInheritanceFromPython(0, 0, 3)
      self.assertEqual(str(the_err.exception), "All genotypes, 0, 0, 3 must be less than the number of genotypes, 3")

   def testSetTimeIntegrationMethod(self):
      self.c.setTimeIntegrationMethod("explicit_euler")
      self.c.setTimeIntegrationMethod("solve_ivp")
      self.c.setTimeIntegrationMethod("runge_kutta4")
      with self.assertRaises(ValueError) as the_err:
         self.c.setTimeIntegrationMethod("crazy_method")
      self.assertEqual(str(the_err.exception), "Time integration method crazy_method not supported")

   def testSetGetMinimumDt(self):
      self.c.setMinimumDt(123.0)
      self.assertEqual(self.c.getMinimumDt(), 123.0)

   def testSetGetAdaptive(self):
      self.c.setAdaptive(12)
      self.assertEqual(self.c.getAdaptive(), 12)

   def testEvolveAdaptiveDt1(self):
      dt = 1.0
      mu_adult = 1.5
      self.c.setNumAges(1)
      self.c.setFecundity(0.0)
      self.c.setAgingRate(0.0)
      self.c.setMuLarvae(2.0) # irrelevant here because num_ages = 1
      self.c.setMuAdult(mu_adult)
      random.seed(1)

      initial_condition = [random.random() for i in range(self.c.getNumberOfPopulations() + self.c.getNumberOfParameters())]
      pap = array.array('f', initial_condition)
      # Explicit-Euler will first give population = initial - dt * mu * initial, which is negative.
      # Adaptive timestepping will cut timestep to 0.9 * dt * initial / (dt * mu * initial) = 0.9 / mu = 0.6.
      # Then explicit will give population = initial - 0.6 * mu * initial = 0.1 * initial
      # The second substep will then try dt = min(1.1 * 0.6, 1 - 0.6) = 0.4, to give population = 0.1 * initial * (1 - dt * mu) = 0.1 * initial * 0.4
      expected_answer = [x * 0.1 * 0.4 for x in initial_condition[:6]] + [initial_condition[-1]]
      self.c.evolve(dt, pap)
      self.assertTrue(arrayfuzzyequal(pap[:-1], expected_answer[:-1], 4E-8))


   def testEvolveAdaptiveDt2(self):
      dt = 2.0
      mu_adult = 1.5
      self.c.setNumAges(1)
      self.c.setFecundity(0.0)
      self.c.setAgingRate(0.0)
      self.c.setMuLarvae(2.0) # irrelevant here because num_ages = 1
      self.c.setMuAdult(mu_adult)
      random.seed(1)

      initial_condition = [random.random() for i in range(self.c.getNumberOfPopulations() + self.c.getNumberOfParameters())]
      pap = array.array('f', initial_condition)
      # Explicit-Euler will first give population = initial - dt * mu * initial, which is negative.
      # Adaptive timestepping will cut timestep to 0.9 * dt * initial / (dt * mu * initial) = 0.9 / mu = 0.6.
      # Then explicit will give population = initial - 0.6 * mu * initial = 0.1 * initial
      # The second substep will then try dt = min(1.1 * 0.6, 2 - 0.6) = 0.66, to give population = 0.1 * initial * (1 - dt * mu) = 0.1 * initial * 0.01
      # The third substep will then try dt = min(1.1 * 0.66, 2 - (0.6 + 0.66)) = 0.726, to give population = 0.1 * 0.01 * initial * (1 - dt * mu), which is negative
      # Adaptive timestepping will cut timestep to 0.9 * dt * 0.1 * 0.01 * initial / (dt * mu * 0.1 * 0.01 * initial) = 0.6
      # then explicit will give population = 0.1 * 0.01 * initial * (1 - 0.6 * mu) = 0.1 * 0.01 * initial * 0.1
      # The next substep will then try dt = min(1.1 * 0.6, 2 - (0.6 + 0.66 + 0.6)) = 0.14
      # Then explicit will give population = 0.1 * 0.1 * 0.01 * initial * (1 - 0.14 * mu) = 0.79 * 0.1 * 0.1 * 0.01 * initial
      expected_answer = [x * 0.79 * 0.1 * 0.1 * 0.01 for x in initial_condition[:6]] + [initial_condition[-1]]
      self.c.evolve(dt, pap)
      self.assertTrue(arrayfuzzyequal(pap[:-1], expected_answer[:-1], 2E-9))

   def testSetGetZeroCutoff(self):
      self.c.setZeroCutoff(123.25)
      self.assertEqual(self.c.getZeroCutoff(), 123.25)

   def testSetGetMinCarryingCapacity(self):
      self.c.setMinCarryingCapacity(45.5)
      self.assertEqual(self.c.getMinCarryingCapacity(), 45.5)


   def testEvolveZeroCutoff(self):
      dt = 1.0
      mu_adult = 1.5
      cutoff = 0.1
      self.c.setNumAges(1)
      self.c.setFecundity(0.0)
      self.c.setAgingRate(0.0)
      self.c.setMuLarvae(2.0) # irrelevant here because num_ages = 1
      self.c.setMuAdult(mu_adult)
      self.c.setZeroCutoff(cutoff)

      initial_condition = list(range(self.c.getNumberOfPopulations() + self.c.getNumberOfParameters()))
      pap = array.array('f', initial_condition)
      # Explicit-Euler will first give population = initial - dt * mu * initial, which is negative.
      # Adaptive timestepping will cut timestep to 0.9 * dt * initial / (dt * mu * initial) = 0.9 / mu = 0.6.
      # Then explicit will give population = initial - 0.6 * mu * initial = 0.1 * initial
      # The second substep will then try dt = min(1.1 * 0.6, 1 - 0.6) = 0.4, to give population = 0.1 * initial * (1 - dt * mu) = 0.1 * initial * 0.4
      expected_answer = [x * 0.1 * 0.4 if x * 0.1 * 0.4 > cutoff else 0 for x in initial_condition[:6]] + [initial_condition[-1]]
      self.c.evolve(dt, pap)
      self.assertTrue(arrayfuzzyequal(pap[:-1], expected_answer[:-1], 1E-7))

   def testEvolveZeroFecundityZeroAging(self):
      dt = 0.01
      self.c.setFecundity(0.0)
      self.c.setAgingRate(0.0)
      self.c.setMuLarvae(0.5)
      self.c.setMuAdult(0.7)
      random.seed(1)

      initial_condition = [random.random() for i in range(self.c.getNumberOfPopulations() + self.c.getNumberOfParameters())]
      pap = array.array('f', initial_condition)
      expected_answer = [x * (1 - 0.5 * dt) for x in initial_condition[:6]] + [x * (1 - 0.7 * dt) for x in initial_condition[6:12]] + [initial_condition[-1]]
      self.c.evolve(dt, pap)
      self.assertTrue(arrayfuzzyequal(pap, expected_answer, 4E-5))

      num_ages = 5
      num_gt = 3
      num_species = 4
      num_sexes = 2
      self.c.setNumSpecies(4)
      self.c.setNumAges(5)
      initial_condition = [random.random() for i in range(self.c.getNumberOfPopulations() + self.c.getNumberOfParameters())]
      expected_answer = [x * (1 - 0.5 * dt) for x in initial_condition[:(num_ages - 1) * num_species * num_gt * num_sexes]] + [x * (1 - 0.7 * dt) for x in initial_condition[(num_ages - 1) * num_species * num_gt * num_sexes : num_ages * num_species * num_gt * num_sexes]] + [initial_condition[-1]]
      pap = array.array('f', initial_condition)
      self.c.evolve(dt, pap)
      self.assertTrue(arrayfuzzyequal(pap, expected_answer, 4E-5))

   def testEvolveZeroFecundityZeroAgingRK(self):
      dt = 0.5
      mu_larvae = 0.5
      mu_adult = 0.7
      self.c.setFecundity(0.0)
      self.c.setAgingRate(0.0)
      self.c.setMuLarvae(mu_larvae)
      self.c.setMuAdult(mu_adult)
      self.c.setTimeIntegrationMethod("runge_kutta4")
      random.seed(1)

      initial_condition = [random.random() for i in range(self.c.getNumberOfPopulations() + self.c.getNumberOfParameters())]
      rhs = [- mu_larvae * x for x in initial_condition[:6]] + [- mu_adult * x for x in initial_condition[6:12]]
      k1 = [dt * r for r in rhs]
      rky = [initial_condition[i] + 0.5 * k1[i] for i in range(12)]
      rhs = [- mu_larvae * y for y in rky[:6]] + [- mu_adult * y for y in rky[6:12]]
      k2 = [dt * r for r in rhs]
      rky = [initial_condition[i] + 0.5 * k2[i] for i in range(12)]
      rhs = [- mu_larvae * y for y in rky[:6]] + [- mu_adult * y for y in rky[6:12]]
      k3 = [dt * r for r in rhs]
      rky = [initial_condition[i] + k3[i] for i in range(12)]
      rhs = [- mu_larvae * y for y in rky[:6]] + [- mu_adult * y for y in rky[6:12]]
      k4 = [dt * r for r in rhs]
      expected_answer = [initial_condition[i] + (1.0 / 6.0) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) for i in range(12)] + [initial_condition[-1]]

      pap = array.array('f', initial_condition)
      self.c.evolve(dt, pap)
      self.assertTrue(arrayfuzzyequal(pap, expected_answer, 1E-7))

   def testEvolveAgingOnly(self):
      dt = 0.01
      aging_rate = 0.25
      self.c.setFecundity(0.0)
      self.c.setAgingRate(aging_rate)
      self.c.setMuLarvae(0.0)
      self.c.setMuAdult(0.0)
      random.seed(1)

      initial_condition = [random.random() for i in range(self.c.getNumberOfPopulations() + self.c.getNumberOfParameters())]
      pap = array.array('f', initial_condition)
      expected_answer = [x * (1 - aging_rate * dt) for x in initial_condition[:6]] + [initial_condition[i] + initial_condition[i - 6] * aging_rate * dt for i in range(6, 12)] + [initial_condition[-1]]
      self.c.evolve(dt, pap)
      self.assertTrue(arrayfuzzyequal(pap, expected_answer, 5E-6))
      conserved = [pap[i] + pap[i + 6] - initial_condition[i] - initial_condition[i + 6] for i in range(6)]
      self.assertTrue(arrayfuzzyequal(conserved, [0.0] * 6, 6E-8))

      num_ages = 5
      num_gt = 3
      num_species = 4
      num_sexes = 2
      self.c.setNumSpecies(num_species)
      self.c.setNumAges(num_ages)
      initial_condition = [random.random() for i in range(self.c.getNumberOfPopulations() + self.c.getNumberOfParameters())]
      expected_answer = [x * (1 - aging_rate * dt) for x in initial_condition[:1 * num_species * num_gt * num_sexes]] + [initial_condition[i] + initial_condition[i - num_species * num_gt * num_sexes] * aging_rate * dt for i in range(1 * num_species * num_gt * num_sexes, num_ages * num_species * num_gt * num_sexes)] + [initial_condition[-1]]
      pap = array.array('f', initial_condition)
      self.c.evolve(dt, pap)
      self.assertTrue(arrayfuzzyequal(pap, expected_answer, 4E-3))
      for sp in range(num_species):
         for sex in range(num_sexes):
            for gt in range(num_gt):
               conserved = 0
               for age in range(num_ages):
                  ind = sp + gt * num_species + sex * num_species * num_gt + age * num_species * num_gt * num_sexes
                  conserved += pap[ind] - initial_condition[ind]
               self.assertTrue((conserved < 3E-7 and conserved > -3E-7))


   def testEvolveSingleAge2(self):
      dt = 5.0
      self.c.setMuLarvae(0.6) # irrelevant here because num_ages=1
      self.c.setMuAdult(0.1)
      self.c.setFecundity(2.0)
      self.c.setAgingRate(0.1234)  # irrelevant here because num_ages=1
      self.c.setNumAges(1)
      self.c.setAccuracy(0.25)

      initial_condition = list(range(6)) + [60.0]
      pap = array.array('f', initial_condition)
      self.c.evolve(dt, pap)
      # this expected answer was worked out laboriously by hand!
      mat = [[-0.1, 0, 0, 0.125, 0.0625, 0],
             [0, -0.1, 0, 0.625, 0.375, 0.125],
             [0, 0, -0.1, 0, 0.3125, 0.625],
             [0, 0, 0, 0.375-0.1, 0.1875, 0],
             [0, 0, 0, 1.875, 1.125-0.1, 0.375],
             [0, 0, 0, 0, 0.9375, 1.875-0.1]]
      expected_result = [initial_condition[i] + dt * sum([mat[i][j] * initial_condition[j] for j in range(6)]) for i in range(6)] + [60.0]
      self.assertTrue(arrayfuzzyequal(pap, expected_result, 1E-8))

   def testEvolveSingleAge1(self):
      dt = 1.0
      self.c.setMuLarvae(0.1)
      self.c.setMuAdult(0.1)
      self.c.setFecundity(0.9)
      self.c.setAgingRate(0.1)
      self.c.setNumAges(1)
      self.c.setAccuracy(0.95)

      initial_condition = [0.1] * self.c.getNumberOfPopulations() + [9.0 / 7.0]
      pap = array.array('f', initial_condition)
      self.c.evolve(dt, pap)
      expected_result_from_nick = [0.10800000, 0.12600000, 0.10800000, 0.10231579, 0.10326316, 0.09094737] + [initial_condition[-1]]
      self.assertTrue(arrayfuzzyequal(pap, expected_result_from_nick, 4E-8))

   def testEvolveTwoAges(self):
      dt = 3.0
      self.c.setMuLarvae(0.0625)
      self.c.setMuAdult(0.09375)
      self.c.setFecundity(2.0)
      self.c.setAgingRate(0.125)
      self.c.setNumAges(2)
      self.c.setAccuracy(0.25)

      initial_condition = [0.25, 1.25, 0.375, 0.625, 1.125, 0.125, 0, 1, 2, 3, 4, 5, 60.0]
      pap = array.array('f', initial_condition)
      self.c.evolve(dt, pap)
      # this expected answer was worked out laboriously by hand!
      mat = [[-0.0625 - 0.125, 0, 0, 0, 0, 0, 0, 0, 0, 0.125 * 1.25, 0.0625 * 1.25, 0],
             [0, -0.0625 - 0.125, 0, 0, 0, 0, 0, 0, 0, 0.625 * 1.25, 0.375 * 1.25, 0.125 * 1.25],
             [0, 0, -0.0625 - 0.125, 0, 0, 0, 0, 0, 0, 0, 0.3125 * 1.25, 0.625 * 1.25],
             [0, 0, 0, -0.0625 - 0.125, 0, 0, 0, 0, 0, 0.375 * 1.25, 0.1875 * 1.25, 0],
             [0, 0, 0, 0, -0.0625 - 0.125, 0, 0, 0, 0, 1.875 * 1.25, 1.125 * 1.25, 0.375 * 1.25],
             [0, 0, 0, 0, 0, -0.0625 - 0.125, 0, 0, 0, 0, 0.9375 * 1.25, 1.875 * 1.25],
             [0.125, 0, 0, 0, 0, 0, -0.09375, 0, 0, 0, 0, 0],
             [0, 0.125, 0, 0, 0, 0, 0, -0.09375, 0, 0, 0, 0],
             [0, 0, 0.125, 0, 0, 0, 0, 0, -0.09375, 0, 0, 0],
             [0, 0, 0, 0.125, 0, 0, 0, 0, 0, -0.09375, 0, 0],
             [0, 0, 0, 0, 0.125, 0, 0, 0, 0, 0, -0.09375, 0],
             [0, 0, 0, 0, 0, 0.125, 0, 0, 0, 0, 0, -0.09375]]

      expected_result = [initial_condition[i] + dt * sum([mat[i][j] * initial_condition[j] for j in range(12)]) for i in range(12)] + [60.0]
      self.assertTrue(arrayfuzzyequal(pap, expected_result, 1E-8))

   def testEvolveThreeAges(self):
      dt = 3.0
      self.c.setMuLarvae(0.0625)
      self.c.setMuAdult(0.09375)
      self.c.setFecundity(2.0)
      self.c.setAgingRate(0.125)
      self.c.setNumAges(3)
      self.c.setAccuracy(0.25)

      initial_condition = [0.125, 0.625, 0.1875, 0.3125, 0.5625, 0.0625,  0.625, 0.125, 0.15625, 0.34375, 0.53125, 0.09375,  0, 1, 2, 3, 4, 5, 60.0]
      pap = array.array('f', initial_condition)
      self.c.evolve(dt, pap)
      # this expected answer was worked out laboriously by hand for TwoAges
      mat = [[-0.0625 - 0.125, 0, 0, 0, 0, 0,  0, 0, 0, 0 ,0, 0,  0, 0, 0, 0.125 * 1.25, 0.0625 * 1.25, 0],
             [0, -0.0625 - 0.125, 0, 0, 0, 0,  0, 0, 0, 0 ,0, 0,  0, 0, 0, 0.625 * 1.25, 0.375 * 1.25, 0.125 * 1.25],
             [0, 0, -0.0625 - 0.125, 0, 0, 0,  0, 0, 0, 0 ,0, 0,  0, 0, 0, 0, 0.3125 * 1.25, 0.625 * 1.25],
             [0, 0, 0, -0.0625 - 0.125, 0, 0,  0, 0, 0, 0 ,0, 0,  0, 0, 0, 0.375 * 1.25, 0.1875 * 1.25, 0],
             [0, 0, 0, 0, -0.0625 - 0.125, 0,  0, 0, 0, 0 ,0, 0,  0, 0, 0, 1.875 * 1.25, 1.125 * 1.25, 0.375 * 1.25],
             [0, 0, 0, 0, 0, -0.0625 - 0.125,  0, 0, 0, 0 ,0, 0,  0, 0, 0, 0, 0.9375 * 1.25, 1.875 * 1.25],
             [0.125, 0, 0, 0, 0, 0,  -0.0625 - 0.125, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0],
             [0, 0.125, 0, 0, 0, 0,  0, -0.0625 - 0.125, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0],
             [0, 0, 0.125, 0, 0, 0,  0, 0, -0.0625 - 0.125, 0, 0, 0,  0, 0, 0, 0, 0, 0],
             [0, 0, 0, 0.125, 0, 0,  0, 0, 0, -0.0625 - 0.125, 0, 0,  0, 0, 0, 0, 0, 0],
             [0, 0, 0, 0, 0.125, 0,  0, 0, 0, 0, -0.0625 - 0.125, 0,  0, 0, 0, 0, 0, 0],
             [0, 0, 0, 0, 0, 0.125,  0, 0, 0, 0, 0, -0.0625 - 0.125,  0, 0, 0, 0, 0, 0],
             [0, 0, 0, 0 ,0, 0,  0.125, 0, 0, 0, 0, 0,  -0.09375, 0, 0, 0, 0, 0],
             [0, 0, 0, 0 ,0, 0,  0, 0.125, 0, 0, 0, 0,  0, -0.09375, 0, 0, 0, 0],
             [0, 0, 0, 0 ,0, 0,  0, 0, 0.125, 0, 0, 0,  0, 0, -0.09375, 0, 0, 0],
             [0, 0, 0, 0 ,0, 0,  0, 0, 0, 0.125, 0, 0,  0, 0, 0, -0.09375, 0, 0],
             [0, 0, 0, 0 ,0, 0,  0, 0, 0, 0, 0.125, 0,  0, 0, 0, 0, -0.09375, 0],
             [0, 0, 0, 0 ,0, 0,  0, 0, 0, 0, 0, 0.125,  0, 0, 0, 0, 0, -0.09375]]

      expected_result = [initial_condition[i] + dt * sum([mat[i][j] * initial_condition[j] for j in range(18)]) for i in range(18)] + [60.0]
      self.assertTrue(arrayfuzzyequal(pap, expected_result, 1E-8))

   def testEvolveTwoAgesMinCC(self):
      dt = 5.0
      self.c.setMuLarvae(0.0625)
      self.c.setMuAdult(0.09375)
      self.c.setAgingRate(0.125)
      self.c.setNumAges(2)
      self.c.setMinCarryingCapacity(100.0) # large minimum carrying capacity

      initial_condition = list(range(12)) + [60.0] # carrying capacity < min_cc, so no new larvae
      pap = array.array('f', initial_condition)
      self.c.evolve(dt, pap)
      mat = [[-0.0625 - 0.125, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
             [0, -0.0625 - 0.125, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
             [0, 0, -0.0625 - 0.125, 0, 0, 0, 0, 0, 0, 0, 0, 0],
             [0, 0, 0, -0.0625 - 0.125, 0, 0, 0, 0, 0, 0, 0, 0],
             [0, 0, 0, 0, -0.0625 - 0.125, 0, 0, 0, 0, 0, 0, 0],
             [0, 0, 0, 0, 0, -0.0625 - 0.125, 0, 0, 0, 0, 0, 0],
             [0.125, 0, 0, 0, 0, 0, -0.09375, 0, 0, 0, 0, 0],
             [0, 0.125, 0, 0, 0, 0, 0, -0.09375, 0, 0, 0, 0],
             [0, 0, 0.125, 0, 0, 0, 0, 0, -0.09375, 0, 0, 0],
             [0, 0, 0, 0.125, 0, 0, 0, 0, 0, -0.09375, 0, 0],
             [0, 0, 0, 0, 0.125, 0, 0, 0, 0, 0, -0.09375, 0],
             [0, 0, 0, 0, 0, 0.125, 0, 0, 0, 0, 0, -0.09375]]
      
      expected_result = [initial_condition[i] + dt * sum([mat[i][j] * initial_condition[j] for j in range(12)]) for i in range(12)] + [60.0]
      self.assertTrue(arrayfuzzyequal(pap, expected_result, 1E-8))


if __name__ == '__main__':
   unittest.main()

