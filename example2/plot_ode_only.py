#!/usr/bin/env python

# plots results of ode_only.py
import sys
import numpy as np
import matplotlib.pyplot as plt

data1 = np.genfromtxt('ode_only_P_1.csv', delimiter = ',', names = True, dtype = float)
data20 = np.genfromtxt('ode_only_P_20.csv', delimiter = ',', names = True, dtype = float)


plt.figure(0)
plt.plot(data1['x'], data1['y'], marker=".", color='r', label='P=1 day')
plt.plot(data20['x'], data20['y'], marker='.', color='g', label='P=20 day')
plt.plot([0.02], [0.1], marker='o', color='k')
plt.xlabel("x", fontsize=22)
plt.ylabel("y", fontsize=22)
plt.legend(fontsize=22)
plt.gca().set_xlim([-0.01, 0.1])
plt.gca().set_ylim([0, 0.25])
plt.title("Evolution dependent on carrying-capacity periodicity")
plt.tight_layout()
plt.savefig("ode_only.png")

plt.figure(1)
x = np.arange(0, 10, 0.01)
y = [0.2 + 0.8 * (1 - (int(i) % 2)) for i in x]
plt.plot(x, y)
plt.xlabel("t (day)", fontsize=22)
plt.ylabel("K", fontsize=22)
tics = [0, 2, 4, 6, 8, 10]
plt.gca().set_xticks(tics)
plt.gca().set_xticklabels([str(t) + "P" for t in tics], fontsize=22)
plt.title("Periodic carrying capacity", fontsize=24)
plt.tight_layout()
plt.savefig("ode_K.png")

sys.exit(0)
