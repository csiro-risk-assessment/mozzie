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
      self.w2 = Wind(os.path.join(findbin, "wind1.csv"), os.path.join(findbin, "wind2_processed.csv"), [[1.0, 0.3], [2.0, 0.7]], self.g2)

      self.w1b = Wind(os.path.join(findbin, "wind1.bin"), os.path.join(findbin, "wind1_processed.bin"), [[1.0, 0.3], [2.0, 0.7]], self.g1)
      self.w1b.setBinaryFileFormat(1)
      self.w2b = Wind(os.path.join(findbin, "wind1.bin"), os.path.join(findbin, "wind2_processed.bin"), [[1.0, 0.3], [2.0, 0.7]], self.g2)
      self.w2b.setBinaryFileFormat(1)
      

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
      for num in range(1, 13):
         num = str(num)
         w = Wind(os.path.join(findbin, "bad_inactive_active" + num + ".csv"), "no_file", [[1, 1]], self.g1)
         with self.assertRaises(ValueError) as the_err:
            w.parseRawFile()
         self.assertEqual(str(the_err.exception), "Header lines in " + os.path.join(findbin, "bad_inactive_active" + num + ".csv") + " must include #xmin=1.0,ymin=2.0,cell_size=3.0,nx=4,ny=3")
      
      for num in range(1, 3):
         num = str(num)
         w = Wind(os.path.join(findbin, "bad_wind" + num + ".csv"), "no_file", [[1, 1]], self.g1)
         with self.assertRaises(ValueError) as the_err:
            w.parseRawFile()
         self.assertEqual(str(the_err.exception), "There must be 8 entries per line in " + os.path.join(findbin, "bad_wind" + num + ".csv"))

      for num in range(3, 5):
         num = str(num)
         w = Wind(os.path.join(findbin, "bad_wind" + num + ".csv"), "no_file", [[1, 1]], self.g1)
         with self.assertRaises(ValueError) as the_err:
            w.parseRawFile()
         self.assertEqual(str(the_err.exception), "There must be 3 data lines in " + os.path.join(findbin, "bad_wind" + num + ".csv"))

   def testBadProcessedFile(self):

      # Here we use the bad_inactive_active*.csv files to test for xmin headerline problems
      for num in range(1, 12):
         num = str(num)
         w = Wind("no_file", os.path.join(findbin, "bad_inactive_active" + num + ".csv"), [[1, 1]], self.g1)
         with self.assertRaises(ValueError) as the_err:
            w.parseProcessedFile()
         self.assertEqual(str(the_err.exception), "Header lines in " + os.path.join(findbin, "bad_inactive_active" + num + ".csv") + " must include #Active cells defined by file None\n#raw_vel_filename=no_file\n#processed_pdf=[[1.0, 1.0]]\n#xmin=1.0,ymin=2.0,cell_size=3.0,nx=4,ny=3")

      for num in ["1", "2", "4", "5", "6", "7"]:
         w = Wind("no_file", os.path.join(findbin, "bad_wind_processed" + num + ".csv"), [[1, 1]], self.g2)
         with self.assertRaises(ValueError) as the_err:
            w.parseProcessedFile()
         self.assertEqual(str(the_err.exception), "Header lines in " + os.path.join(findbin, "bad_wind_processed" + num + ".csv") + " must include #Active cells defined by file inactive_active.csv\n#raw_vel_filename=no_file\n#processed_pdf=[[1.0, 1.0]]\n#xmin=1.0,ymin=2.0,cell_size=3.0,nx=4,ny=3")

      for num in ["1", "2", "3"]:
         wb = Wind("no_file", os.path.join(findbin, "bad_wind_processed" + num + ".bin"), [[1, 1]], self.g2)
         wb.setBinaryFileFormat(1)
         with self.assertRaises(ValueError) as the_err:
            wb.parseProcessedFile()
         self.assertEqual(str(the_err.exception), "Header lines in " + os.path.join(findbin, "bad_wind_processed" + num + ".bin") + " must include #Active cells defined by file inactive_active.csv\n#raw_vel_filename=no_file\n#processed_pdf=[[1.0, 1.0]]\n#xmin=1.0,ymin=2.0,cell_size=3.0,nx=4,ny=3")

      for num in ["8", "9"]:
         w = Wind("no_file", os.path.join(findbin, "bad_wind_processed" + num + ".csv"), [[1, 1]], self.g2)
         with self.assertRaises(ValueError) as the_err:
            w.parseProcessedFile()
         self.assertEqual(str(the_err.exception), "There must be 3 entries per line in " + os.path.join(findbin, "bad_wind_processed" + num + ".csv"))

      for num in ["11", "12", "13", "14"]:
         w = Wind("no_file", os.path.join(findbin, "bad_wind_processed" + num + ".csv"), [[1, 1]], self.g2)
         with self.assertRaises(ValueError) as the_err:
            w.parseProcessedFile()
         self.assertEqual(str(the_err.exception), "Data in " + os.path.join(findbin, "bad_wind_processed" + num + ".csv") + " is incorrectly bounded.  Bad data line number = 0")

      for num in ["11", "12", "13", "14"]:
         wb = Wind("no_file", os.path.join(findbin, "bad_wind_processed" + num + ".bin"), [[1, 1]], self.g2)
         wb.setBinaryFileFormat(1)
         with self.assertRaises(ValueError) as the_err:
            wb.parseProcessedFile()
         self.assertEqual(str(the_err.exception), "Data in " + os.path.join(findbin, "bad_wind_processed" + num + ".bin") + " is incorrectly bounded.  Bad data line number = 0")


   def testGetRawWindFilename(self):
      self.assertEqual(self.w1.getRawWindFilename(), os.path.join(findbin, "wind1.csv"))
      self.assertEqual(self.w2.getRawWindFilename(), os.path.join(findbin, "wind1.csv"))

   def testGetProcessedWindFilename(self):
      self.assertEqual(self.w1.getProcessedWindFilename(), os.path.join(findbin, "wind1_processed.csv"))
      self.assertEqual(self.w2.getProcessedWindFilename(), os.path.join(findbin, "wind2_processed.csv"))

   def testGetPDF(self):
      self.assertEqual(self.w1.getPDF(), [[1.0, 0.3], [1.0, 0.7]])
      self.assertEqual(self.w2.getPDF(), [[1.0, 0.3], [1.0, 0.7]])
      
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
      
      self.w1b.parseRawFile()
      self.assertTrue(arrayequal(self.w1b.getAdvectionFrom(), [0, 1, 2, 4, 5, 6, 6, 7, 8, 8, 9, 9, 10, 10, 11, 11]))
      self.assertTrue(arrayequal(self.w1b.getAdvectionTo(), [0, 5, 0, 5, 5, 2, 0, 3, 1, 5, 10, 11, 11, 10, 10, 11]))
      self.assertTrue(arrayfuzzyequal(self.w1b.getAdvectionP(), [1.0, 1.0, 1.0, 1.0, 1.0, 0.3, 0.7, 0.3, 0.3, 0.7, 0.3, 0.7, 0.3, 0.7, 0.3, 0.7], 1E-5))
      self.assertEqual(self.w1b.getNumAdvection(), 16)
      
      self.w2.parseRawFile()
      self.assertTrue(arrayequal(self.w2.getAdvectionFrom(), [0, 1, 4, 4, 5, 5, 6, 6]))
      self.assertTrue(arrayequal(self.w2.getAdvectionTo(), [0, 0, 1, 0, 6, 5, 5, 6]))
      self.assertTrue(arrayfuzzyequal(self.w2.getAdvectionP(), [1.0, 1.0, 0.3, 0.7, 0.3, 0.7, 0.3, 0.7], 1E-5))
      self.assertEqual(self.w2.getNumAdvection(), 8)

      self.w2b.parseRawFile()
      self.assertTrue(arrayequal(self.w2b.getAdvectionFrom(), [0, 1, 4, 4, 5, 5, 6, 6]))
      self.assertTrue(arrayequal(self.w2b.getAdvectionTo(), [0, 0, 1, 0, 6, 5, 5, 6]))
      self.assertTrue(arrayfuzzyequal(self.w2b.getAdvectionP(), [1.0, 1.0, 0.3, 0.7, 0.3, 0.7, 0.3, 0.7], 1E-5))
      self.assertEqual(self.w2b.getNumAdvection(), 8)

   def testOutputProcessedCSV(self):
      w3 = Wind(os.path.join(findbin, "wind1.csv"), os.path.join(findbin, "wind3_processed.csv"), [[1.0, 0.3], [2.0, 0.7]], self.g2)
      for (windobj, fn) in [(w3, "inactive_active.csv")]:
         windobj.parseRawFile()
         windobj.outputProcessedCSV()
         with open(os.path.join(findbin, "wind3_processed.csv"), 'r') as f:
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

      self.w1b.parseRawFile()
      self.w1b.outputProcessedCSV()
      self.w1b.parseProcessedFile()
      self.assertTrue(arrayequal(self.w1b.getAdvectionFrom(), [0, 1, 2, 4, 5, 6, 6, 7, 8, 8, 9, 9, 10, 10, 11, 11]))
      self.assertTrue(arrayequal(self.w1b.getAdvectionTo(), [0, 5, 0, 5, 5, 2, 0, 3, 1, 5, 10, 11, 11, 10, 10, 11]))
      self.assertTrue(arrayfuzzyequal(self.w1b.getAdvectionP(), [1.0, 1.0, 1.0, 1.0, 1.0, 0.3, 0.7, 0.3, 0.3, 0.7, 0.3, 0.7, 0.3, 0.7, 0.3, 0.7], 1E-5))
      self.assertEqual(self.w1b.getNumAdvection(), 16)
      self.assertEqual(self.w1b.getProcessedDataComputed(), 1)

      self.w2b.parseRawFile()
      self.w2b.outputProcessedCSV()
      self.w2b.parseProcessedFile()
      self.assertTrue(arrayequal(self.w2b.getAdvectionFrom(), [0, 1, 4, 4, 5, 5, 6, 6]))
      self.assertTrue(arrayequal(self.w2b.getAdvectionTo(), [0, 0, 1, 0, 6, 5, 5, 6]))
      self.assertTrue(arrayfuzzyequal(self.w2b.getAdvectionP(), [1.0, 1.0, 0.3, 0.7, 0.3, 0.7, 0.3, 0.7], 1E-5))
      self.assertEqual(self.w2b.getNumAdvection(), 8)
      self.assertEqual(self.w2b.getProcessedDataComputed(), 1)


   def testBinaryFormat(self):
      with self.assertRaises(ValueError) as the_err:
         self.w1.setBinaryFileFormat(3)
      self.assertEqual(str(the_err.exception), "FileFormat must be 0 or 1");
      self.assertEqual(self.w1.getBinaryFileFormat(), 0)
      self.w1.setBinaryFileFormat(1)
      self.assertEqual(self.w1.getBinaryFileFormat(), 1)

      



if __name__ == '__main__':
   unittest.main()

