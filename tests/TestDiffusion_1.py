import os
import sys
import unittest

# so we can find our ../code no matter how we are called
findbin = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(findbin + "/../code")

import array
from grid import Grid
from cellDynamics import CellDynamicsStatic15_9_3_2
from spatialDynamics import SpatialDynamics
from populationsAndParameters import PopulationsAndParameters

def arrayfuzzyequal(a, b, eps):
   return (len(a) == len(b)) and all([(a[i] > b[i] - eps and a[i] < b[i] + eps) for i in range(min(len(a), len(b)))])

class TestDiffusion_1(unittest.TestCase):

   def setUp(self):
      # grid that is 2kmx2km long with 8x8 active cells
      nx = 8
      self.grid = Grid(-0.25 * nx, -0.25 * nx, 0.5, nx, nx)

      self.cell = CellDynamicsStatic15_9_3_2()

      # all the population and parameters
      self.all_quantities = PopulationsAndParameters(self.grid, self.cell)
      # centre cell starts with nonzero population
      pop = list(range(self.cell.getNumberOfPopulations() + self.cell.getNumberOfParameters()))
      self.all_quantities.setPopulationAndParameters((nx * (nx + 1)) // 2, pop)

      self.spatial = SpatialDynamics(self.grid, self.all_quantities)

   def testExcept(self):
      with self.assertRaises(ValueError) as the_err:
         self.spatial.diffuseVarying(1, array.array('f', [1.0] * 10))
      self.assertEqual(str(the_err.exception), "diffusion_coeffs incorrectly sized")

   def testDiffuse1(self):
      # diffuse with timestep = 0.5 and diffusion coefficient = 0.075, hence diffusion_d = 0.6
      self.spatial.diffuse(0.5, 0.075)
      self.spatial.outputCSV(os.path.join(findbin, "2D_diffusion_out_1.csv"), 2, "0", "")
      with open(os.path.join(findbin, "2D_diffusion_out_1.csv")) as f:
         data = f.readlines()
      for row in [2, 3, 4, 8, 9]:
         self.assertTrue(arrayfuzzyequal([float(d) for d in data[row].strip().split(",")], [0.0] * 8, 1E-5))
      self.assertTrue(arrayfuzzyequal([float(d) for d in data[5].strip().split(",")], [0, 0, 0, 0, 0.3, 0, 0, 0], 1E-5))
      self.assertTrue(arrayfuzzyequal([float(d) for d in data[6].strip().split(",")], [0, 0, 0, 0.3, 0.8, 0.3, 0, 0], 1E-5))
      self.assertTrue(arrayfuzzyequal([float(d) for d in data[7].strip().split(",")], [0, 0, 0, 0, 0.3, 0, 0, 0], 1E-5))
      
      self.spatial.diffuse(0.5, 0.075)
      self.spatial.outputCSV(os.path.join(findbin, "2D_diffusion_out_2.csv"), 2, "0", "")
      with open(os.path.join(findbin, "2D_diffusion_out_2.csv")) as f:
         data = f.readlines()
      for row in [2, 3, 9]:
         self.assertTrue(arrayfuzzyequal([float(d) for d in data[row].strip().split(",")], [0.0] * 8, 1E-5))
      self.assertTrue(arrayfuzzyequal([float(d) for d in data[4].strip().split(",")], [0, 0, 0, 0, 0.045, 0, 0, 0], 1E-5))
      self.assertTrue(arrayfuzzyequal([float(d) for d in data[5].strip().split(",")], [0, 0, 0, 0.09, 0.24, 0.09, 0, 0], 1E-5))
      self.assertTrue(arrayfuzzyequal([float(d) for d in data[6].strip().split(",")], [0, 0, 0.045, 0.24, 0.5, 0.24, 0.045, 0], 1E-5))
      self.assertTrue(arrayfuzzyequal([float(d) for d in data[7].strip().split(",")], [0, 0, 0, 0.09, 0.24, 0.09, 0, 0], 1E-5))
      self.assertTrue(arrayfuzzyequal([float(d) for d in data[8].strip().split(",")], [0, 0, 0, 0, 0.045, 0, 0, 0], 1E-5))

   def testExcept(self):
      with self.assertRaises(ValueError) as the_err:
         self.spatial.diffuseVarying(1, array.array('f', [1.0] * 10))
      self.assertEqual(str(the_err.exception), "diffusion_coeffs incorrectly sized")

   def testDiffuseVarying(self):
      # diffuse with timestep = 0.5 and varying diffusion coefficient = 0.075 (so diffusion_d = 0.6), or diffusion coefficient = 0.1 (so diffusion_d = 0.8)
      self.spatial.diffuseVarying(0.5, array.array('f', [0.01 * i for i in range(9)]))
      diffusion_d = [0.01 * i * 4 * 0.5 / (0.5 * 0.5) for i in range(9)]
      diffusing_indices = [0, 2, 3, 4, 5, 6, 8, 10, 11]

      for i in range(9):
         self.spatial.outputCSV(os.path.join(findbin, "2D_diffusion_out_" + str(i) + ".csv"), diffusing_indices[i], "0", "")
         with open(os.path.join(findbin, "2D_diffusion_out_" + str(i) + ".csv")) as f:
            data = f.readlines()
            for row in [2, 3, 4, 8, 9]:
               self.assertTrue(arrayfuzzyequal([float(d) for d in data[row].strip().split(",")], [0.0] * 8, 1E-5))
            r = diffusing_indices[i] * diffusion_d[i] / 4.0
            f = diffusing_indices[i] * (1.0 - diffusion_d[i])
            self.assertTrue(arrayfuzzyequal([float(d) for d in data[5].strip().split(",")], [0, 0, 0, 0, r, 0, 0, 0], 1E-5))
            self.assertTrue(arrayfuzzyequal([float(d) for d in data[6].strip().split(",")], [0, 0, 0, r, f, r, 0, 0], 1E-5))
            self.assertTrue(arrayfuzzyequal([float(d) for d in data[7].strip().split(",")], [0, 0, 0, 0, r, 0, 0, 0], 1E-5))

if __name__ == '__main__':
   unittest.main()

