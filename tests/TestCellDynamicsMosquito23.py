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
      self.assertEqual(self.c.getHybridisationRate(0, 0, 0), 1.0)
      self.c.setNumSpecies(3)
      for i in range(3):
         for j in range(3):
            for k in range(3):
               self.assertTrue(arrayfuzzyequal([self.c.getHybridisationRate(i, j, k)], [1 if (i == j and j == k) else 0], 1E-8))
      self.c.setHybridisationRate(0, 0, 0, 1.23)
      self.c.setHybridisationRate(0, 1, 1, -1.23)
      self.c.setHybridisationRate(0, 2, 2, 3.33)
      self.c.setHybridisationRate(1, 0, 0, 1.4)
      self.c.setHybridisationRate(1, 2, 1, -0.5)
      self.c.setHybridisationRate(2, 0, 2, 0.8)
      self.assertTrue(arrayfuzzyequal([self.c.getHybridisationRate(0, 0, 0)], [1.23], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getHybridisationRate(0, 1, 1)], [-1.23], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getHybridisationRate(0, 2, 2)], [3.33], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getHybridisationRate(1, 0, 0)], [1.4], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getHybridisationRate(1, 2, 1)], [-0.5], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getHybridisationRate(2, 0, 2)], [0.8], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getHybridisationRate(1, 1, 1)], [1.0], 1E-6))
      self.assertTrue(arrayfuzzyequal([self.c.getHybridisationRate(2, 2, 2)], [1.0], 1E-6))
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


   def testEvolveSingleAge(self):
      sys.stderr.write("SKIPPING")
      return
      dt = 50.0
      self.c.setMuLarvae(0.1)
      self.c.setMuAdult(0.1)
      self.c.setFecundity(0.9)
      self.c.setAgingRate(0.1)
      self.c.setNumAges(1)
      self.c.setAccuracy(0.95)

      initial_condition = [0.1] * self.c.getNumberOfPopulations() + [9.0 / 7.0]
      pap = array.array('f', initial_condition)
      self.c.evolve(dt, pap)
      expected_result_from_nick = [0.32054381, 0.14635508, 0.01695156, 0.24774694, 0.03637915, 0.00153114] + [9.0 / 7.0]
      self.assertTrue(arrayfuzzyequal(pap, expected_result_from_nick, 1E-4))


if __name__ == '__main__':
   unittest.main()

