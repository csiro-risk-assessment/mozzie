import os
import sys
import unittest

# so we can find our ../code no matter how we are called
findbin = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(findbin + "/../code")

from grid import Grid
from cell import Cell
from spatial import Spatial
from populations import Populations

def arrayfuzzyequal(a, b, eps):
   return all([(a[i] > b[i] - eps and a[i] < b[i] + eps) for i in range(0, len(a))])

class TestDiffusion_2(unittest.TestCase):

   def setUp(self):
      # grid that is 2kmx2km long with 8x8 active cells
      nx = 8
      self.grid = Grid(-0.25 * nx, -0.25 * nx, 0.5, nx, nx)
      # now make some inactive
      self.grid.setActiveAndInactive(os.path.join(findbin, "inactive_active_8x8.csv"))
      
      # all the populations
      self.all_pops = Populations(self.grid)
      # centre cell starts with nonzero population
      pop = list(range(Cell().getNumberOfPopulations()))
      self.all_pops.setPopulation(17, pop)

      # initialise the spatial structure with timestep = 0.5 and diffusion coefficient = 0.075
      # hence diffusion_d = 0.6
      self.spatial = Spatial(0.5, 0.075, self.grid, self.all_pops)


   def testDiffuse1(self):
      self.spatial.diffuse()
      self.spatial.outputCSV(os.path.join(findbin, "2D_diffusion_out_3.csv"), 2)
      with open(os.path.join(findbin, "2D_diffusion_out_3.csv")) as f:
         data = f.readlines()
      for row in [0, 1, 2, 3, 6, 7]:
         self.assertTrue(arrayfuzzyequal([float(d) for d in data[row].strip().split(",")], [0.0] * 8, 1E-5))
      self.assertTrue(arrayfuzzyequal([float(d) for d in data[4].strip().split(",")], [0, 0, 0, 0.3, 0.8, 0.3, 0, 0], 1E-5))
      self.assertTrue(arrayfuzzyequal([float(d) for d in data[5].strip().split(",")], [0, 0, 0, 0, 0.3, 0, 0, 0], 1E-5))
      
      self.spatial.diffuse()
      self.spatial.outputCSV(os.path.join(findbin, "2D_diffusion_out_4.csv"), 2)
      with open(os.path.join(findbin, "2D_diffusion_out_4.csv")) as f:
         data = f.readlines()
      for row in [0, 1, 2, 3, 7]:
         self.assertTrue(arrayfuzzyequal([float(d) for d in data[row].strip().split(",")], [0.0] * 8, 1E-5))
      self.assertTrue(arrayfuzzyequal([float(d) for d in data[4].strip().split(",")], [0, 0, 0.045, 0.24, 0.455, 0.24, 0, 0], 1E-5))
      self.assertTrue(arrayfuzzyequal([float(d) for d in data[5].strip().split(",")], [0, 0, 0, 0.09, 0.24, 0.09, 0, 0], 1E-5))
      self.assertTrue(arrayfuzzyequal([float(d) for d in data[6].strip().split(",")], [0, 0, 0, 0, 0.045, 0, 0, 0], 1E-5))

if __name__ == '__main__':
   unittest.main()

