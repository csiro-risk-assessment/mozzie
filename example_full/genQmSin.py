import os
import numpy as np
import pandas as pd
from scipy.ndimage import distance_transform_edt

# Set this to the path of your working directory
working_dir = "/scratch3/bee069/bee069/example_full" # "SET PATH HERE"
os.chdir(os.path.expanduser(working_dir))

# Create a directory for saving the CSV for qm daily between 01-01-2022 and 31-12-2023
writing_dir = "qm"
if not os.path.exists(writing_dir):
    os.makedirs(os.path.join(working_dir, writing_dir))

# Header of the CSV file
header = "#xmin=-1799134.0,ymin=-1299134.0,cell_size=5000.0,nx=100,ny=100"

# Read the active cells grid of 100 by 100 cells where active cells are represented by 1 and inactive cells are represented by 0
# Inactive cells form the shape of a "river"
active = pd.read_csv(os.path.join(working_dir, "activeDemo.csv"), skiprows=1, header=None)
active_mat = active.to_numpy()
active_mat = np.flipud(active_mat)  # Reverse the order of rows

# Compute the distance transform, it is going to be used in the vertical shift in the sinus formula
# The distance grid is a grid of of 100 by 100 cells where each value represented the distance to the closest inactive cells
# A distance of zero is equivalent to an inactive cell
distance_mat = distance_transform_edt(active_mat)
# Multiply the distances by 10 so they have some weight
distance_mat = distance_mat*10

# Define the years
years_simulated = list(range(2022, 2024))

# Define the species
species = ["Aa", "Ac", "Ag"]

# Define sinusoidal parameters, one for each species (Aa, Ac, Ag)
# Sinusoidal function is daily varying abundance at equilibrium of each species, a function of time in number of days
# f(t) = amplitude*sin(2*pi*frequency*t+phase)+shift
# The parameters are chosen so that:
# An. arabiensis abundance is varying between ~ 1300 and ~ 7300 females 
# An. coluzzii abundance is varying between ~ 3000 and ~ 83000 females 
# An. gambiae ss abundance is varying between ~ 1000 and ~ 45000 females
# The peaks of abundance of An. coluzzii and An. gambiae ss are shifted by 40 days approximately with respect to An. arabiensis peaks
# amplitude is the difference between the maximum and the minimum divided by 2
amplitude = [-3000, 40000, 22000] 
# frequency f = 1/T where T is the period, the period is approximately 365
frequency = [0.00274, 0.00274, 0.00274]
# horizontal shift or phase shift
phase = [0.9, -1/5, -1.5]
# vertical shift, the vertical shift must be sufficient to make the sinusoidal function always positive
# another vertical shift is going to be added later using the distance to the river
shift = [4300, 43000, 23000]

# Initialize parameters of the model
# Reproduction rate in number of females
lamb = 4.5
# Death rate, one for each species (Aa, Ac, Ag) 
z = np.repeat(0.125, 3)
# Competition matrix
A = np.array([[1.00000000, 0.34, 0.25],
                   [0.34, 1.00000000, 0.28],
                   [0.25, 0.28, 1.000000]])
         
# Create matrices
vec_X = np.zeros((100 * 100, 3))
mat_X = np.zeros((100, 100))
mat_d = np.zeros((100, 100))

# no of day
tt = 0
# Loop through days and years
for year in years_simulated:
    for doy in range(365):
        tt = tt + 1
        # 3 filenames, one for each species (Aa, Ac, Ag) 
        new_files = [f"{writing_dir}/qm.{species}.{year}.{doy+1}.csv" for species in ["Aa", "Ac", "Ag"]]

        # Check if the files already exist
        if not all(os.path.exists(file) for file in new_files):
            for sp in range(3):  # Loop through species
                # Matrix of distance is 100 by 100 cells
                mat_dist = np.array(distance_mat).reshape(100, 100)

                #A vertical shift is going to be added to the sinusoidal function depending on the distance grid and the species
                if sp == 0:
                    # The vertical shift for the cells with Aa far from the "river" is going to be increase by a quadratic factor of the distance
                    mat_d = (mat_dist / 20) ** 2
                elif sp == 1:
                    # The vertical shift for the cells with Aa close from the "river" is going to be increase by a linear factor of the distance
                    mat_d = -4 * mat_dist
                elif sp == 2:
                    # The vertical shift for the cells with Ag far from the "river" is going to be increase by a cubic factor of the distance
                    mat_d = (mat_dist / 35) ** 3

                # Calculate the sinusoidal for all grid to represent current species abundance in each cell
                mat_X = amplitude[sp] * np.sin(2 * np.pi * frequency[sp] * tt + phase[sp]) + (shift[sp] + mat_d)
                # Inactive cells from original grid have an abundance of zero
                mat_X = mat_X * active_mat
                # Represent the grid in a 1D vector, so that vec_X contains abundances for all three species in three different columns
                vec_X[:, sp] = mat_X.flatten()

                # Print a message if a value of abundance is inferior to zero (should not happen with the parameters values we chose)
                if np.any(vec_X < 0):
                    print(f"big problem for timestep: {tt}")
            
            # The values qm, species-specific, are obtained from the daily abundances X at equilibrium (X*)
            # X* = inv(A)*diag((lamb - z)/z)*qm
            # qm = A*diag(z/(lamb - z))
            vec_q = vec_X.dot(A)
            vec_q = vec_q.dot(np.diag(z / (lamb - z)))

            # Save each matrix representing qm for each species in a file, one for each day  
            for sp in range(3):  # Write CSV files
                mat_qm = vec_q[:, sp].reshape(100, 100)
                mat_qm = np.flipud(mat_qm)  # Flip rows
                file_path = new_files[sp]

                with open(file_path, 'w') as f:
                    f.write(header + "\n")
                    np.savetxt(f, mat_qm, delimiter=',', fmt='%f')
