import os
import sys
import unittest

# so we can find our ../code no matter how we are called
findbin = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(findbin + "/../code")

from wind import Wind
from grid import Grid


grid = Grid(-4614384.0, -3967418.0, 5000.0, 1517, 1667)
wind = Wind("/scratch1/projects/mozzie/era5/2015.01.01.windvec.csv", "tmp.csv", [[1.0, 1.0]], grid)
wind.parseRawFile()
