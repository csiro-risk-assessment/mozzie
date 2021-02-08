#!/usr/bin/env python

######################################################################
#
# Converts two files, fn_base.eastwind.csv and fn_base.northwind.csv,
# in which the velocity is specified in m/s,
# into one file, fn_base.windvec_m_per_day.csv,
# in which the velocity is specified in m/day, rounded to the nearest integer
#
######################################################################

import sys

if len(sys.argv) != 2:
    sys.stdout.write("You must supply this python script with a filename base, for instance:\n")
    sys.stdout.write("   " + sys.argv[0] + " 1990.06.22\n")
    sys.exit(1)

# Read input data
fn_base = sys.argv[1]
sys.stdout.write("Reading " + fn_base + ".eastwind.csv and " + fn_base + ".northwind.csv\n")
convert_to_m_per_day = 3600 * 24
sys.stdout.write("These are assumed to specify wind velocity in m/s\n")
na = '0'
sys.stdout.write("All 'NA' strings will be replaced by " + na + "\n")
with open(fn_base + ".eastwind.csv", "r") as f:
    east = [v_x.replace("NA", na) for v_x in f.readlines()]
with open(fn_base + ".northwind.csv", "r") as f:
    north = [v_y.replace("NA", na) for v_y in f.readlines()]

# Check data
if east[0] != north[0]:
    sys.stdout.write("Header in the two files must be identical!\n")
    sys.exit(1)
ny = len(east) - 1
if len(north) - 1 != ny:
    sys.stdout.write("Number of lines (ny) in each file must be identical!\n")
    sys.exit(1)
nx = 0
if (ny > 0):
    nx = len(east[1].split(","))
    if len(north[1].split(",")) != nx:
        sys.stdout.write("nx must be the same in each file\n")
        sys.exit(1)

# Convert
out_fn = fn_base + ".windvec_m_per_day.csv"
sys.stdout.write("Writing result to " + out_fn + "\n")
with open(out_fn, "w") as f:
    f.write(east[0])
    for y in range(1, ny + 1):
        v_x = [int(v * convert_to_m_per_day + (0.5 if v > 0 else -0.5)) for v in map(float, east[y].strip().split(","))]
        v_y = [int(v * convert_to_m_per_day + (0.5 if v > 0 else -0.5)) for v in map(float, north[y].strip().split(","))]
        out_str = ""
        for x in range(nx):
            out_str += str(v_x[x]) + "," + str(v_y[x]) + ","
        f.write(out_str[:-1] + "\n")

sys.exit(0)
