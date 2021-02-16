import os
import sys
import unittest
import array

# so we can find our ../code no matter how we are called
findbin = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(findbin + "/../code")

from cellDynamics import CellDynamicsBeeton2_2

def arrayequal(a, b):
   return all([a[i] == b[i] for i in range(0, len(a))])

def arrayfuzzyequal(a, b, eps):
   return all([(a[i] > b[i] - eps and a[i] < b[i] + eps) for i in range(0, len(a))])

class TestCellDynamicsBeeton2_2(unittest.TestCase):

   def setUp(self):
      self.c = CellDynamicsBeeton2_2()

   def testSetGetMuX(self):
      self.c.setMuX(1234.0)
      self.assertEqual(self.c.getMuX(), 1234.0)
      
   def testSetGetMuY(self):
      self.c.setMuY(12345.0)
      self.assertEqual(self.c.getMuY(), 12345.0)
      
   def testSetGetGaX(self):
      self.c.setGaX(123.0)
      self.assertEqual(self.c.getGaX(), 123.0)
      
   def testSetGetGaY(self):
      self.c.setGaY(123.0)
      self.assertEqual(self.c.getGaY(), 123.0)
      
   def testSetGetAxy(self):
      self.c.setAxy(12.0)
      self.assertEqual(self.c.getAxy(), 12.0)
      
   def testSetGetAyx(self):
      self.c.setAyx(-12.0)
      self.assertEqual(self.c.getAyx(), -12.0)
      
   def testSetGetW(self):
      self.c.setW(-123.0)
      self.assertEqual(self.c.getW(), -123.0)

   def testSetGetSmall(self):
      self.c.setSmall(1234.0)
      self.assertEqual(self.c.getSmall(), 1234.0)

   def testEvolve(self):
      x = 1.0
      y = 2.0
      kx = 3.0
      ky = 4.0
      mux = 5.0
      axy = 6.0
      w = 7.0
      gax = 8.0
      muy = 9.0
      ayx = 10.0
      gay = 11.0
      small = 0.2
      pap = array.array('f', [x, y, kx, ky])
      self.c.setMuX(mux)
      self.c.setMuY(muy)
      self.c.setGaX(gax)
      self.c.setGaY(gay)
      self.c.setAxy(axy)
      self.c.setAyx(ayx)
      self.c.setW(w)
      self.c.setSmall(small)
      dt = 0.02
      self.c.evolve(dt, pap)

      # this should be the solution:
      f = -mux + (1 - (x + axy * y) / (kx + small)) * (x / (x + w * y + small)) * gax
      g = -muy + (1 - (ayx * x + y) / (ky + small)) * (gay + w * x * gax / (x + w * y + small))
      rhsx = f * x
      rhsy = g * y
      xnew = x + dt * rhsx
      ynew = y + dt * rhsy

      self.assertTrue(arrayfuzzyequal(pap, [xnew, ynew, kx, ky], 1E-6))

if __name__ == '__main__':
   unittest.main()

