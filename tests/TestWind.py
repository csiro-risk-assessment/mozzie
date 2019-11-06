import os
import sys
import unittest

# so we can find our ../code no matter how we are called
findbin = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(findbin + "/../code")

from wind import Wind
from grid import Grid

def arrayequal(a, b):
   return all([a[i] == b[i] for i in range(0, len(a))])

def arrayfuzzyequal(a, b, eps):
   return all([(a[i] > b[i] - eps and a[i] < b[i] + eps) for i in range(0, len(a))])

class TestWind(unittest.TestCase):

   def setUp(self):
      self.g1 = Grid(1.0, 2.0, 3.0, 4, 3)
      self.w1 = Wind(os.path.join(findbin, "wind1.csv"), os.path.join(findbin, "wind1_processed.csv"), [[1.0, 0.3], [2.0, 0.7]], self.g1)
      self.g2 = Grid(1.0, 2.0, 3.0, 4, 3)
      self.g2.setActiveAndInactive(os.path.join(findbin, "inactive_active.csv"))
      self.w2 = Wind(os.path.join(findbin, "wind1.csv"), os.path.join(findbin, "wind1_processed.csv"), [[1.0, 0.3], [2.0, 0.7]], self.g2)
      

   def testBadPDF(self):
      pdf_error_string = "PDF defining the advecting velocities in terms of the raw velocity must be of the form [[t0, p0], [t1, p1], [t2, p2]...], where p is the probability that an advecting mosquito will stay in the wind-stream for time t.  The t's must be positive and monotonically increasing.  The p's must sum to unity"

      with self.assertRaises(ValueError) as the_err:
         Wind("no_file", "processed_no_file", [[0.0, 1.0]], self.g1)
      self.assertEqual(str(the_err.exception), pdf_error_string)

      with self.assertRaises(ValueError) as the_err:
         Wind("no_file", "processed_no_file", [[-1.0, 1.0]], self.g1)
      self.assertEqual(str(the_err.exception), pdf_error_string)

      with self.assertRaises(ValueError) as the_err:
         Wind("no_file", "processed_no_file", [["one", 1.0]], self.g1)
      self.assertEqual(str(the_err.exception), pdf_error_string)

      with self.assertRaises(ValueError) as the_err:
         Wind("no_file", "processed_no_file", [[0.6, 0.3], [0.1, 0.7]], self.g1)
      self.assertEqual(str(the_err.exception), pdf_error_string)

      with self.assertRaises(ValueError) as the_err:
         Wind("no_file", "processed_no_file", [[0.1, 0.2], [0.4, 0.7]], self.g1)
      self.assertEqual(str(the_err.exception), pdf_error_string)


   def testBadRawFile(self):

      # Here we use the bad_inactive_active*.csv files to test for headerline problems
      for num in range(1, 12):
         num = str(num)
         w = Wind(os.path.join(findbin, "bad_inactive_active" + num + ".csv"), "no_file", [[1, 1]], self.g1)
         with self.assertRaises(ValueError) as the_err:
            w.parseRawFile()
         self.assertEqual(str(the_err.exception), "The header line in " + os.path.join(findbin, "bad_inactive_active" + num + ".csv") + " does not match #xmin=1.0,ymin=2.0,cell_size=3.0,nx=4,ny=3")

      num = "12"
      w = Wind(os.path.join(findbin, "bad_inactive_active" + num + ".csv"), "no_file", [[1, 1]], self.g1)
      with self.assertRaises(ValueError) as the_err:
         w.parseRawFile()
      self.assertEqual(str(the_err.exception), "Header line of the form #xmin=... not found in " + os.path.join(findbin, "bad_inactive_active" + num + ".csv"))

      # now use bad_wind*.csv to check for data problems
      data_string_append = " must be CSV formatted.  Each line in the file corresponds to a row (constant y) of cells, so must contain 8 entries (separated by commas).  Entries come (vel_x, vel_y) pairs, so a row is of the form vx0,vy0,vx1,vy1,vx2,vy2,....  There must be 3 such rows.  The first row corresponds to cells at y=ymin, the next row at y=ymin+cell_size, etc"
      for num in range(1, 5):
         num = str(num)
         w = Wind(os.path.join(findbin, "bad_wind" + num + ".csv"), "no_file", [[1, 1]], self.g1)
         with self.assertRaises(ValueError) as the_err:
            w.parseRawFile()
         self.assertEqual(str(the_err.exception), "Data in " + os.path.join(findbin, "bad_wind" + num + ".csv") + data_string_append)

   def testBadProcessedFile(self):

      # Here we use the bad_inactive_active*.csv files to test for xmin headerline problems
      for num in range(1, 12):
         num = str(num)
         w = Wind("no_file", os.path.join(findbin, "bad_inactive_active" + num + ".csv"), [[1, 1]], self.g1)
         with self.assertRaises(ValueError) as the_err:
            w.parseProcessedFile()
         self.assertEqual(str(the_err.exception), "The header line in " + os.path.join(findbin, "bad_inactive_active" + num + ".csv") + " does not match #xmin=1.0,ymin=2.0,cell_size=3.0,nx=4,ny=3")


      w = Wind("no_file", os.path.join(findbin, "bad_wind_processed1.csv"), [[1, 1]], self.g2)
      with self.assertRaises(ValueError) as the_err:
         w.parseProcessedFile()
      self.assertEqual(str(the_err.exception), "Active cell filename in " + os.path.join(findbin, "bad_wind_processed1.csv") + " incompatible with that specified in the grid, which is inactive_active.csv")

      w = Wind("no_file", os.path.join(findbin, "bad_wind_processed2.csv"), [[1, 1]], self.g2)
      with self.assertRaises(ValueError) as the_err:
         w.parseProcessedFile()
      self.assertEqual(str(the_err.exception), "Raw velocity filename in " + os.path.join(findbin, "bad_wind_processed2.csv") + " incompatible with that specified in the wind constructor, which is no_file")

      w = Wind("no_file", os.path.join(findbin, "bad_wind_processed3.csv"), [[1, 1]], self.g2)
      with self.assertRaises(ValueError) as the_err:
         w.parseProcessedFile()
      self.assertEqual(str(the_err.exception), "PDF specified in " + os.path.join(findbin, "bad_wind_processed3.csv") + " incompatible with that specified in the wind constructor.  Remember the file's version is processed to give dt values rather than t values")

      for num in ["4", "5", "6", "7"]:
         w = Wind("no_file", os.path.join(findbin, "bad_wind_processed" + num + ".csv"), [[1, 1]], self.g2)
         with self.assertRaises(ValueError) as the_err:
            w.parseProcessedFile()
         self.assertEqual(str(the_err.exception), "Not all header lines found in " + os.path.join(findbin, "bad_wind_processed" + num + ".csv"))

      for num in ["8", "9", "10"]:
         w = Wind("no_file", os.path.join(findbin, "bad_wind_processed" + num + ".csv"), [[1, 1]], self.g2)
         with self.assertRaises(ValueError) as the_err:
            w.parseProcessedFile()
         self.assertEqual(str(the_err.exception), "Data in " + os.path.join(findbin, "bad_wind_processed" + num + ".csv") + " must be CSV formatted.  Each line must contain integer1,integer2,float , where integer1 = active cell index from where mosquitoes are advecting; integer2 = active cell index to which mosquitoes are advecting; float = probability of this occuring")

      num = "11"
      w = Wind("no_file", os.path.join(findbin, "bad_wind_processed" + num + ".csv"), [[1, 1]], self.g2)
      with self.assertRaises(ValueError) as the_err:
         w.parseProcessedFile()
      self.assertEqual(str(the_err.exception), "Data in " + os.path.join(findbin, "bad_wind_processed" + num + ".csv") + " is incorrectly bounded.  Bad line = 1234,2,0.1")

      num = "12"
      w = Wind("no_file", os.path.join(findbin, "bad_wind_processed" + num + ".csv"), [[1, 1]], self.g2)
      with self.assertRaises(ValueError) as the_err:
         w.parseProcessedFile()
      self.assertEqual(str(the_err.exception), "Data in " + os.path.join(findbin, "bad_wind_processed" + num + ".csv") + " is incorrectly bounded.  Bad line = 0,12345,0.1")

      num = "13"
      w = Wind("no_file", os.path.join(findbin, "bad_wind_processed" + num + ".csv"), [[1, 1]], self.g2)
      with self.assertRaises(ValueError) as the_err:
         w.parseProcessedFile()
      self.assertEqual(str(the_err.exception), "Data in " + os.path.join(findbin, "bad_wind_processed" + num + ".csv") + " is incorrectly bounded.  Bad line = 0,0,-0.1")

      num = "14"
      w = Wind("no_file", os.path.join(findbin, "bad_wind_processed" + num + ".csv"), [[1, 1]], self.g2)
      with self.assertRaises(ValueError) as the_err:
         w.parseProcessedFile()
      self.assertEqual(str(the_err.exception), "Data in " + os.path.join(findbin, "bad_wind_processed" + num + ".csv") + " is incorrectly bounded.  Bad line = 0,0,1.1")



   def testGetRawWindFilename(self):
      self.assertEqual(self.w1.getRawWindFilename(), os.path.join(findbin, "wind1.csv"))

   def testGetProcessedWindFilename(self):
      self.assertEqual(self.w1.getProcessedWindFilename(), os.path.join(findbin, "wind1_processed.csv"))

   def testGetPDF(self):
      self.assertEqual(self.w1.getPDF(), [[1.0, 0.3], [1.0, 0.7]])
      
   def testGetProcessedDataComputed(self):
      self.assertEqual(self.w1.getProcessedDataComputed(), 0)
      self.w1.parseRawFile()
      self.assertEqual(self.w1.getProcessedDataComputed(), 1)

      self.assertEqual(self.w2.getProcessedDataComputed(), 0)
      self.w2.parseProcessedFile()
      self.assertEqual(self.w2.getProcessedDataComputed(), 1)

   def testAdvection(self):
      self.w1.parseRawFile()
      self.assertTrue(arrayequal(self.w1.getAdvectionFrom(), [0, 1, 2, 4, 5, 6, 6, 7, 8, 8, 9, 9, 10, 10, 11, 11]))
      self.assertTrue(arrayequal(self.w1.getAdvectionTo(), [0, 5, 0, 5, 5, 2, 0, 3, 1, 5, 10, 11, 11, 10, 10, 11]))
      self.assertTrue(arrayfuzzyequal(self.w1.getAdvectionP(), [1.0, 1.0, 1.0, 1.0, 1.0, 0.3, 0.7, 0.3, 0.3, 0.7, 0.3, 0.7, 0.3, 0.7, 0.3, 0.7], 1E-5))
      self.assertEqual(self.w1.getNumAdvection(), 16)
      
      self.w2.parseRawFile()
      self.assertTrue(arrayequal(self.w2.getAdvectionFrom(), [0, 1, 4, 4, 5, 5, 6, 6]))
      self.assertTrue(arrayequal(self.w2.getAdvectionTo(), [0, 0, 1, 0, 6, 5, 5, 6]))
      self.assertTrue(arrayfuzzyequal(self.w2.getAdvectionP(), [1.0, 1.0, 0.3, 0.7, 0.3, 0.7, 0.3, 0.7], 1E-5))
      self.assertEqual(self.w2.getNumAdvection(), 8)

   def testOutputProcessedCSV(self):
      for (windobj, fn) in [(self.w1, "None"), (self.w2, "inactive_active.csv")]:
         windobj.parseRawFile()
         windobj.outputProcessedCSV()
         with open(os.path.join(findbin, "wind1_processed.csv"), 'r') as f:
            data = f.readlines()
         self.assertEqual(data[1], "#Active cells defined by file " + fn + "\n")
         self.assertEqual(data[2], "#xmin=1.0,ymin=2.0,cell_size=3.0,nx=4,ny=3\n")
         self.assertEqual(data[3], "#raw_vel_filename=wind1.csv\n")
         self.assertEqual(data[4], "#processed_pdf=" + str(windobj.getPDF()) + "\n")
         afr = []
         ato = []
         app = []
         for line in data[7:]:
            line = line.strip().split(",")
            afr.append(int(line[0]))
            ato.append(int(line[1]))
            app.append(float(line[2]))
         self.assertTrue(arrayequal(windobj.getAdvectionFrom(), afr))
         self.assertTrue(arrayequal(windobj.getAdvectionTo(), ato))
         self.assertTrue(arrayfuzzyequal(windobj.getAdvectionP(), app, 1E-5))

   def testParseProcessedFile(self):
      self.w1.parseRawFile()
      self.w1.outputProcessedCSV()
      self.w1.parseProcessedFile()
      self.assertTrue(arrayequal(self.w1.getAdvectionFrom(), [0, 1, 2, 4, 5, 6, 6, 7, 8, 8, 9, 9, 10, 10, 11, 11]))
      self.assertTrue(arrayequal(self.w1.getAdvectionTo(), [0, 5, 0, 5, 5, 2, 0, 3, 1, 5, 10, 11, 11, 10, 10, 11]))
      self.assertTrue(arrayfuzzyequal(self.w1.getAdvectionP(), [1.0, 1.0, 1.0, 1.0, 1.0, 0.3, 0.7, 0.3, 0.3, 0.7, 0.3, 0.7, 0.3, 0.7, 0.3, 0.7], 1E-5))
      self.assertEqual(self.w1.getNumAdvection(), 16)
      self.assertEqual(self.w1.getProcessedDataComputed(), 1)

      self.w2.parseRawFile()
      self.w2.outputProcessedCSV()
      self.w2.parseProcessedFile()
      self.assertTrue(arrayequal(self.w2.getAdvectionFrom(), [0, 1, 4, 4, 5, 5, 6, 6]))
      self.assertTrue(arrayequal(self.w2.getAdvectionTo(), [0, 0, 1, 0, 6, 5, 5, 6]))
      self.assertTrue(arrayfuzzyequal(self.w2.getAdvectionP(), [1.0, 1.0, 0.3, 0.7, 0.3, 0.7, 0.3, 0.7], 1E-5))
      self.assertEqual(self.w2.getNumAdvection(), 8)
      self.assertEqual(self.w2.getProcessedDataComputed(), 1)


      



if __name__ == '__main__':
   unittest.main()

