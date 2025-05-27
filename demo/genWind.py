######################################################
# Creates daily wind velocities between
# 2022-01-01 and 2023-12-31
# It is assumed that:
# - The wind blows northeasterly in July and August at 10m/s
# - The wind blows southwesterly in January at 10m/s
# - The wind is set to 0m/s the other months
import os
import sys
import numpy as np

sys.stdout.write("Initialising paths and output directory\n")
######################################################
# Paths and filenames
#
# Set the working directory, which will contain the output directory
# Change working_dir to suit your needs.
working_dir = os.path.dirname(os.path.realpath(sys.argv[0]))
output_dir = "wind"

# Change working directory, and potentially create output directory
os.chdir(os.path.expanduser(working_dir))
if not os.path.exists(output_dir):
    os.makedirs(os.path.join(working_dir, output_dir))

sys.stdout.write("Computing wind velocity, and outputting\n")
# Define month days (not adjusting for leap years yet)
month_days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

years = range(2022, 2024)
months = range(1, 13)

# The wind blows northeasterly in July and August at 10m/s
# The wind blows southwesterly in January at 10m/s
# The wind is set to 0m/s the other months
# Iterate over years, months, and days
for y in years:
    sys.stdout.write("  year " + str(y) + "\n")
    for m in months:
        for d in range(1, month_days[m-1] + 1):
            wind_east = np.zeros((100, 100))
            wind_north = np.zeros((100, 100))

            # Wind direction adjustments
            # The magnitude of the speed is 10m/s in the indicated direction
            if m in [7, 8]:  # July and August: Northeasterly wind
                wind_east += 10 * np.cos(np.pi / 4)
                wind_north += 10 * np.sin(np.pi / 4)
            elif m == 1:  # January: Southwesterly wind
                wind_east += 10 * np.cos(5 * np.pi / 4)
                wind_north += 10 * np.sin(5 * np.pi / 4)

            # Convert wind speed from m/s to m/day
            wind_east = np.round(86400 * wind_east)
            wind_north = np.round(86400 * wind_north)

            # Flip vertically
            wind_east = np.flipud(wind_east)
            wind_north = np.flipud(wind_north)

            # Create filename
            newfile = os.path.join(output_dir, f"{y:02d}.{m:02d}.{d:02d}.windvec_m_per_day.csv")

            # Create header
            header = [
                "#Active cells defined by file active.csv",
                f"#raw_vel_filename={newfile}",
                "#processed_pdf=[[0.375, 1.0]]",
                "#xmin=-25000.0,ymin=-25000.0,cell_size=500.0,nx=100,ny=100"
            ]

            # Reshape matrix for saving
            wind_en = np.column_stack((wind_east, wind_north)).reshape(wind_east.shape[0], -1, order='F')

            # Write to file
            with open(newfile, 'w') as f:
                f.write("\n".join(header) + "\n")
                np.savetxt(f, wind_en, delimiter=",", fmt="%d")

sys.stdout.write("Done\n")
