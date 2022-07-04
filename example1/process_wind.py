# Processes wind vectors given by Nick into a single file, raw_wind.csv, and changing to m/day
import os
import sys

sys.stdout.write("Creating one wind file, raw_wind.csv, from Nick's two files.  Setting all NA to 0\n")

sys.stdout.write("Reading eastwind.csv\n")
with open("eastwind.csv", "r") as f:
   east = f.readlines()
sys.stdout.write("Reading northwind.csv\n")
with open("northwind.csv", "r") as f:
   north = f.readlines()

to_m_per_day = 24 * 3600
sys.stdout.write("Creating raw_wind.csv\n")
f = open("raw_wind.csv", "w")
f.write(east[0])
for i in range(1, len(east)):
   e = east[i].strip().split(",")
   n = north[i].strip().split(",")
   for j in range(len(e) - 1):
      try:
         f.write(str(float(e[j]) * to_m_per_day) + "," + str(float(n[j]) * to_m_per_day) + ",")
      except:
         f.write("0,0,")
   j = len(e) - 1
   try:
      f.write(str(float(e[j]) * to_m_per_day) + "," + str(float(n[j]) * to_m_per_day) + "\n")
   except:
      f.write("0,0\n")
f.close()
