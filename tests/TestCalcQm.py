import os
import sys
import unittest
import random


# so we can find our ../code no matter how we are called
findbin = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(findbin + "/../code")

from grid import Grid
from cellDynamics import CellDynamicsMosquitoBH26Delay
from spatialDynamics import SpatialDynamics
from populationsAndParameters import PopulationsAndParameters


def arrayfuzzyequal(a, b, eps):
   return (len(a) == len(b)) and all([(a[i] > b[i] - eps and a[i] < b[i] + eps) for i in range(min(len(a), len(b)))])

class TestEvolveCells(unittest.TestCase):

   def setUp(self):
      self.num_active_cells = 6
      self.grid = Grid(0, 0, 1, self.num_active_cells, 1)

      self.num_species = 3
      delay_days = 0 # Note: in principal can have > 0 delay days, but the code below is slightly slower and more complicated, since it is necessary to put the equilibrium populations into the correct slots
      death_rate_ww = [0.1, 0.2, 0.3] # death rate of wild-types of each species, measured in 1/day.  Below we assume the other genotypes have the same death rate: if a bad assumption then just modify death_rate variable below.  Must be changed if num_species changes from 3
      competition = [[1, 0.4, 0.2], [0.4, 1, 0.3], [0.2, 0.3, 1]] # competition between subspecies.  Must be changed if num_species changes from 3
      emergence_rate = [9.0, 7.0, 6.0] # emergence rate for each species.  Must be changed if num_species changes from 3
      activity = [[1, 0.1, 0.4], [0.1, 1, 0.2], [0.4, 0.2, 1]] # activity[female_of_species1][male_of_species2] is activity level in the proportionate mixing.  Must be changed if num_species changes from 3
      hM = 0.125 # reduction R = (1 - s * h) * (1 - s * h) for genotype = wc, for instance: see the reduct function below.  Using that function, it is assumed that hM = hF
      sM = 0.25 # reduction R = (1 - s * h) * (1 - s * h) for genotype = wc, for instance: see the reduct function below.  Using that function, it is assumed that sM = sF
      hybridisation = [[[1, 0, 0], [0, 1, 0], [0, 0, 1]], [[1, 0, 0], [0, 1, 0], [0, 0, 1]], [[1, 0, 0], [0, 1, 0], [0, 0, 1]]] # hybridisation[mM][mF][m] = prob that offspring of species m results from male of species mM and female of species mF.  The current value means offspring is always same as mF.  Must be changed if num_species changes from 3
      offspring_modifier = [[[0.5, 1, 1], [1, 1.5, 1], [1, 1, 0.7]], [[0.5, 1, 1], [1, 1.5, 1], [1, 1, 0.7]]] # offspring_modifier[s][mM][mF] = suppression (if <1, or increased vigour if >1) of sppring of sex s that arises from male of species mM and female of species mF
      sex_ratio = 0.7 # probability that offspring of wc or cc fathers are male
      female_bias = 0.6 # probability that offspring of (mother wc or cc + father ww) is female
      m_w = 0
      m_c = 0
      small_value = 1E-10 # when populations, carrying capacity, etc get smaller than this, set them to zero
      death_rate = [[death_rate_ww] * 6] * 2
      def reduct(genotype):
         if genotype == 0 or genotype == 2 or genotype == 5:
            return 1.0
         elif genotype == 1 or genotype == 4:
            return 1.0 - sM * hM
         elif genotype == 3:
            return 1.0 - sM
         raise ValueError("invalid genotype")
      reduction = [[reduct(gM) * reduct(gF) for gM in range(6)] for gF in range(6)]
      self.cell = CellDynamicsMosquitoBH26Delay(num_species = self.num_species, delay = delay_days, current_index = 0, death_rate = death_rate, competition = competition, emergence_rate = emergence_rate, activity = activity, reduction = reduction, hybridisation = hybridisation, offspring_modifier = offspring_modifier, sex_ratio = sex_ratio, female_bias = female_bias, m_w = m_w, m_c = m_c)

      # all the population and parameters
      self.all_quantities = PopulationsAndParameters(self.grid, self.cell)
      # include some interesting equilibrium population numbers and thus carrying capacities
      self.cell_cc = [[0] * self.num_species * 6 * 2 + [0] * self.num_species for i in range(self.num_active_cells)]
      for active_cell in range(self.num_active_cells):
         for m in range(self.num_species):
            males_plus_females = 1 + random.random()
            for g in range(1): # only ww genotype
               for s in range(2):
                  ind = m + g * self.num_species + s * self.num_species * 6
                  self.cell_cc[active_cell][ind] = 0.5 * males_plus_females
         self.all_quantities.setPopulationAndParameters(active_cell, self.cell_cc[active_cell])

      self.spatial = SpatialDynamics(self.grid, self.all_quantities)


   def testEvolveCells1(self):
      # calculate qm
      self.spatial.calcQm()
      # now qm should be in the correct slots in self.all_quantities

      # evolve with large timestep
      self.spatial.evolveCells(1E6)

      # check that the equilibrium in self.cc is preserved
      num_params_pops = self.cell.getNumberOfPopulations() + self.cell.getNumberOfParameters()
      for active_cell in range(self.num_active_cells):
         for m in range(self.num_species):
            for g in range(1): # only ww genotype
               for s in range(2):
                  ind = m + g * self.num_species + s * self.num_species * 6
                  self.assertTrue(abs(self.all_quantities.getQuantities()[active_cell * num_params_pops + ind] - self.cell_cc[active_cell][ind]), 1E-6)
      
if __name__ == '__main__':
   unittest.main()

