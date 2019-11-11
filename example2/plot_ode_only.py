#!/usr/bin/env python

# plots results of ode_only.py
import sys
import numpy as np
import matplotlib.pyplot as plt

data1 = np.genfromtxt('ode_only_P_1.csv', delimiter = ',', names = True, dtype = float)
data20 = np.genfromtxt('ode_only_P_20.csv', delimiter = ',', names = True, dtype = float)

plt.figure()
plt.plot(data1['x'], data1['y'], marker=".", color='r', label='P=1')
plt.plot(data20['x'], data20['y'], marker='.', color='g', label='P=20')
plt.plot([0.02], [0.1], marker='o', color='k')
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.gca().set_xlim([-0.01, 0.1])
plt.gca().set_ylim([0, 0.25])
plt.title("Evolution dependent on carrying-capacity periodicity")
plt.savefig("ode_only.png")

sys.exit(0)
