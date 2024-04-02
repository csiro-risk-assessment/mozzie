import os
import sys
import unittest
import array

# so we can find our ../code no matter how we are called
findbin = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(findbin + "/../code")

from spatialDependence import SpatialDependence

delete_output = True

def arrayequal(a, b):
   return (len(a) == len(b)) and all([a[i] == b[i] for i in range(min(len(a), len(b)))])
def arrayfuzzyequal(a, b, eps):
   return (len(a) == len(b)) and all([(a[i] > b[i] - eps and a[i] < b[i] + eps) for i in range(min(len(a), len(b)))])

class TestSpatialDependence(unittest.TestCase):

   def setUp(self):
      self.s1 = SpatialDependence(1.0, 2.0, 3.0, 4, 3)

   def testBadFileType(self):
      with self.assertRaises(ValueError) as the_err:
         self.s1.parse(os.path.join(findbin, "inactive_active.csv"), "bad_filetype", [])
      self.assertEqual(str(the_err.exception), "filetype not recognized")
      with self.assertRaises(ValueError) as the_err:
         self.s1.parseWithPython(os.path.join(findbin, "inactive_active.csv"), "bad_filetype", [])
      self.assertEqual(str(the_err.exception), "filetype not recognized")

   def testBadFileName(self):
      with self.assertRaises(IOError) as the_err:
         self.s1.parse("badfilename", "active_inactive", [])
      self.assertEqual(str(the_err.exception), "Cannot open or read badfilename")
      with self.assertRaises(IOError) as the_err:
         self.s1.parseWithPython("badfilename", "active_inactive", [])
      self.assertEqual(str(the_err.exception), "Cannot open or read badfilename")

   def testBadHeader(self):
      expected_error = "Header lines in " + os.path.join(findbin, "inactive_active.csv") + " must include #missing header\n#xmin=1.0,ymin=2.0,cell_size=3.0,nx=4,ny=3"
      with self.assertRaises(Exception) as the_err:
         self.s1.parse(os.path.join(findbin, "inactive_active.csv"), "active_inactive", ["#missing header"])
      self.assertEqual(str(the_err.exception), expected_error)
      with self.assertRaises(Exception) as the_err:
         self.s1.parseWithPython(os.path.join(findbin, "inactive_active.csv"), "active_inactive", ["#missing header"])
      self.assertEqual(str(the_err.exception), expected_error)

      expected_error = "Header lines in " + os.path.join(findbin, "generic_float.csv") + " must include #missing header\n#xmin=1.0,ymin=2.0,cell_size=3.0,nx=4,ny=3"
      with self.assertRaises(Exception) as the_err:
         self.s1.parse(os.path.join(findbin, "generic_float.csv"), "generic_float", ["#missing header"])
      self.assertEqual(str(the_err.exception), expected_error)
      with self.assertRaises(Exception) as the_err:
         self.s1.parseWithPython(os.path.join(findbin, "generic_float.csv"), "generic_float", ["#missing header"])
      self.assertEqual(str(the_err.exception), expected_error)

      expected_error = "Header lines in " + os.path.join(findbin, "wind1.csv") + " must include #missing header\n#xmin=1.0,ymin=2.0,cell_size=3.0,nx=4,ny=3"
      with self.assertRaises(Exception) as the_err:
         self.s1.parse(os.path.join(findbin, "wind1.csv"), "wind_raw", ["#missing header"])
      self.assertEqual(str(the_err.exception), expected_error)
      with self.assertRaises(Exception) as the_err:
         self.s1.parseWithPython(os.path.join(findbin, "wind1.csv"), "wind_raw", ["#missing header"])
      self.assertEqual(str(the_err.exception), expected_error)

      expected_error = "Header lines in " + os.path.join(findbin, "wind2_processed_static.csv") + " must include #missing header\n#xmin=1.0,ymin=2.0,cell_size=3.0,nx=4,ny=3"
      with self.assertRaises(Exception) as the_err:
         self.s1.parse(os.path.join(findbin, "wind2_processed_static.csv"), "wind_processed", ["#missing header"])
      self.assertEqual(str(the_err.exception), expected_error)
      with self.assertRaises(Exception) as the_err:
         self.s1.parseWithPython(os.path.join(findbin, "wind2_processed_static.csv"), "wind_processed", ["#missing header"])
      self.assertEqual(str(the_err.exception), expected_error)

   def testBadActiveInactive(self):
      fn = os.path.join(findbin, "bad_inactive_active13.csv")
      with self.assertRaises(ValueError) as the_err:
         self.s1.parse(fn, "active_inactive", [])
      self.assertEqual(str(the_err.exception), "There must be 4 entries per line in " + fn)
      with self.assertRaises(ValueError) as the_err:
         self.s1.parseWithPython(fn, "active_inactive", [])
      self.assertEqual(str(the_err.exception), "There must be 4 entries per line in " + fn)

      fn = os.path.join(findbin, "bad_inactive_active14.csv")
      with self.assertRaises(ValueError) as the_err:
         self.s1.parse(fn, "active_inactive", [])
      self.assertEqual(str(the_err.exception), "There must be 4 entries per line in " + fn)
      with self.assertRaises(ValueError) as the_err:
         self.s1.parseWithPython(fn, "active_inactive", [])
      self.assertEqual(str(the_err.exception), "There must be 4 entries per line in " + fn)

      fn = os.path.join(findbin, "bad_inactive_active15.csv")
      with self.assertRaises(ValueError) as the_err:
         self.s1.parse(fn, "active_inactive", [])
      self.assertEqual(str(the_err.exception), "There must be 3 data lines in " + fn)
      with self.assertRaises(ValueError) as the_err:
         self.s1.parseWithPython(fn, "active_inactive", [])
      self.assertEqual(str(the_err.exception), "There must be 3 data lines in " + fn)

      fn = os.path.join(findbin, "bad_inactive_active16.csv")
      with self.assertRaises(ValueError) as the_err:
         self.s1.parse(fn, "active_inactive", [])
      self.assertEqual(str(the_err.exception), "The data entries in " + fn + " must be either 0 or 1")
      with self.assertRaises(ValueError) as the_err:
         self.s1.parseWithPython(fn, "active_inactive", [])
      self.assertEqual(str(the_err.exception), "The data entries in " + fn + " must be either 0 or 1")

   def testParseActiveInactive(self):
      self.s1.parse(os.path.join(findbin, "inactive_active.csv"), "active_inactive", [])
      self.assertTrue(arrayequal([1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1], self.s1.getData0()))
      self.s1.parseWithPython(os.path.join(findbin, "inactive_active.csv"), "active_inactive", [])
      self.assertTrue(arrayequal([1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1], self.s1.getData0()))

   def testBadGenericFloat(self):
      fn = os.path.join(findbin, "bad_inactive_active13.csv")
      with self.assertRaises(ValueError) as the_err:
         self.s1.parse(fn, "generic_float", [])
      self.assertEqual(str(the_err.exception), "There must be 4 entries per line in " + fn)
      with self.assertRaises(ValueError) as the_err:
         self.s1.parseWithPython(fn, "generic_float", [])
      self.assertEqual(str(the_err.exception), "There must be 4 entries per line in " + fn)

      fn = os.path.join(findbin, "bad_inactive_active14.csv")
      with self.assertRaises(ValueError) as the_err:
         self.s1.parse(fn, "generic_float", [])
      self.assertEqual(str(the_err.exception), "There must be 4 entries per line in " + fn)
      with self.assertRaises(ValueError) as the_err:
         self.s1.parseWithPython(fn, "generic_float", [])
      self.assertEqual(str(the_err.exception), "There must be 4 entries per line in " + fn)

      fn = os.path.join(findbin, "bad_inactive_active15.csv")
      with self.assertRaises(ValueError) as the_err:
         self.s1.parse(fn, "generic_float", [])
      self.assertEqual(str(the_err.exception), "There must be 3 data lines in " + fn)
      with self.assertRaises(ValueError) as the_err:
         self.s1.parseWithPython(fn, "generic_float", [])
      self.assertEqual(str(the_err.exception), "There must be 3 data lines in " + fn)

   def testParseGenericFloat(self):
      self.s1.parse(os.path.join(findbin, "generic_float.csv"), "generic_float", [])
      self.assertTrue(arrayequal([1.0, 0.0, -1.0, 1.25, 2.0, 3.0, 4.0, 5.0, -3.0, -4.0, -6.0, -8.5], self.s1.getData0()))
      self.s1.parseWithPython(os.path.join(findbin, "generic_float.csv"), "generic_float", [])
      self.assertTrue(arrayequal([1.0, 0.0, -1.0, 1.25, 2.0, 3.0, 4.0, 5.0, -3.0, -4.0, -6.0, -8.5], self.s1.getData0()))

   def testBadRawWind(self):
      fn = os.path.join(findbin, "bad_wind2.csv")
      with self.assertRaises(ValueError) as the_err:
         self.s1.parse(fn, "wind_raw", [])
      self.assertEqual(str(the_err.exception), "There must be 8 entries per line in " + fn)
      with self.assertRaises(ValueError) as the_err:
         self.s1.parseWithPython(fn, "wind_raw", [])
      self.assertEqual(str(the_err.exception), "There must be 8 entries per line in " + fn)

      fn = os.path.join(findbin, "bad_wind3.csv")
      with self.assertRaises(ValueError) as the_err:
         self.s1.parse(fn, "wind_raw", [])
      self.assertEqual(str(the_err.exception), "There must be 3 data lines in " + fn)
      with self.assertRaises(ValueError) as the_err:
         self.s1.parseWithPython(fn, "wind_raw", [])
      self.assertEqual(str(the_err.exception), "There must be 3 data lines in " + fn)

      fn = os.path.join(findbin, "bad_wind4.csv")
      with self.assertRaises(ValueError) as the_err:
         self.s1.parse(fn, "wind_raw", [])
      self.assertEqual(str(the_err.exception), "There must be 3 data lines in " + fn)
      with self.assertRaises(ValueError) as the_err:
         self.s1.parseWithPython(fn, "wind_raw", [])
      self.assertEqual(str(the_err.exception), "There must be 3 data lines in " + fn)

   def testParseWindRaw(self):
      self.s1.parse(os.path.join(findbin, "wind1.csv"), "wind_raw", [])
      self.assertTrue(arrayequal([0, 0, -6, 2, 3, -0.125, 0, 0, 3, 3, 3, -3], self.s1.getData0()))
      self.assertTrue(arrayequal([0, 3, 0, 0, 1, 0.125, -3, -3, -6, -0.5, -0.5, -0.5], self.s1.getData1()))
      self.s1.parseWithPython(os.path.join(findbin, "wind1.csv"), "wind_raw", [])
      self.assertTrue(arrayequal([0, 0, -6, 2, 3, -0.125, 0, 0, 3, 3, 3, -3], self.s1.getData0()))
      self.assertTrue(arrayequal([0, 3, 0, 0, 1, 0.125, -3, -3, -6, -0.5, -0.5, -0.5], self.s1.getData1()))

   def testBadProcessedWind(self):
      fn = os.path.join(findbin, "bad_wind_processed8.csv")
      with self.assertRaises(ValueError) as the_err:
         self.s1.parse(fn, "wind_processed", [])
      self.assertEqual(str(the_err.exception), "There must be 3 entries per line in " + fn)
      with self.assertRaises(ValueError) as the_err:
         self.s1.parseWithPython(fn, "wind_processed", [])
      self.assertEqual(str(the_err.exception), "There must be 3 entries per line in " + fn)

      fn = os.path.join(findbin, "bad_wind_processed9.csv")
      with self.assertRaises(ValueError) as the_err:
         self.s1.parse(fn, "wind_processed", [])
      self.assertEqual(str(the_err.exception), "There must be 3 entries per line in " + fn)
      with self.assertRaises(ValueError) as the_err:
         self.s1.parseWithPython(fn, "wind_processed", [])
      self.assertEqual(str(the_err.exception), "There must be 3 entries per line in " + fn)

   def testParseWindProcessed(self):
      self.s1.parse(os.path.join(findbin, "wind2_processed_static.csv"), "wind_processed", [])
      self.assertTrue(arrayequal([0, 1, 4, 4, 5, 5, 6, 6], self.s1.getData0()))
      self.assertTrue(arrayequal([0, 0, 1, 0, 6, 5, 5, 6], self.s1.getData1()))
      self.assertTrue(arrayfuzzyequal([1.0, 1.0, 0.3, 0.7, 0.3, 0.7, 0.3, 0.7], self.s1.getData2(), 1E-5))
      self.s1.parseWithPython(os.path.join(findbin, "wind2_processed_static.csv"), "wind_processed", [])
      self.assertTrue(arrayequal([0, 1, 4, 4, 5, 5, 6, 6], self.s1.getData0()))
      self.assertTrue(arrayequal([0, 0, 1, 0, 6, 5, 5, 6], self.s1.getData1()))
      self.assertTrue(arrayfuzzyequal([1.0, 1.0, 0.3, 0.7, 0.3, 0.7, 0.3, 0.7], self.s1.getData2(), 1E-5))

   def testRestrictToActive(self):
      global_index = array.array("I", [1, 2, 3, 4, 5, 6, 7, 9])

      self.s1.parse(os.path.join(findbin, "inactive_active.csv"), "active_inactive", [])
      self.s1.restrictToActive(global_index)
      self.assertTrue(arrayequal([0, 1, 1, 1, 0, 1, 0, 0], self.s1.getData0()))

      self.s1.parse(os.path.join(findbin, "generic_float.csv"), "generic_float", [])
      self.s1.restrictToActive(global_index)
      self.assertTrue(arrayequal([0.0, -1.0, 1.25, 2.0, 3.0, 4.0, 5.0, -4.0], self.s1.getData0()))

      self.s1.parse(os.path.join(findbin, "wind1.csv"), "wind_raw", [])
      self.s1.restrictToActive(global_index)
      self.assertTrue(arrayequal([0.0, -6.0, 2.0, 3.0, -0.125, 0.0, 0.0, 3.0], self.s1.getData0()))
      self.assertTrue(arrayequal([3.0, 0.0, 0.0, 1.0, 0.125, -3.0, -3.0, -0.5], self.s1.getData1()))

   def testOutputCSV(self):
      num_cells = 3 * 4
      active_list = [num_cells, 0, 1, 2,
                     num_cells, num_cells, 3, num_cells,
                     4, num_cells, 5, 6] # elements that = num_cells mean an inactive cell
      global_index = [1, 2, 3, 6, 8, 10, 11]
      self.s1.parse(os.path.join(findbin, "generic_float.csv"), "generic_float", [])
      self.s1.restrictToActive(array.array("I", global_index))
      out_fn = os.path.join(findbin, "testOutputCSV.csv")
      self.s1.outputCSV(out_fn, array.array("I", active_list), "999", "#additional header\n")
      with open(out_fn, "r") as f:
         data = f.readlines()
      self.assertEqual(len(data), 6)
      self.assertEqual(data[1].strip(), "#xmin=1.0,ymin=2.0,cell_size=3.0,nx=4,ny=3")
      self.assertEqual(data[2].strip(), "#additional header")
      self.assertTrue(arrayequal([999, 0.0, -1.0, 1.25], list(map(float, data[3].strip().split(",")))))
      self.assertTrue(arrayequal([999, 999, 4.0, 999], list(map(float, data[4].strip().split(",")))))
      self.assertTrue(arrayequal([-3.0, 999, -6.0, -8.5], list(map(float, data[5].strip().split(",")))))
      os.remove(out_fn)

if __name__ == '__main__':
   unittest.main()

