import numpy as np
import math 
import matplotlib.pyplot as plt
import least_squares_approximation as lsa_lib

file = open("../data-sets/data1.txt")
lines = file.readlines()

xdata = np.zeros([len(lines)])
ydata = np.zeros([len(lines)])
yest = np.zeros([len(lines)])

for i in range(len(lines)):
   line = lines[i].strip().split(",")
   xdata[i] = line[0]
   ydata[i] = line[1]

def f(x):
   return np.array([1, x, math.cos(2.3*x), math.sin(2.3*x)])

coeff = lsa_lib.lsa(f, xdata, ydata)

for i in range(len(lines)):
   yest[i] = lsa_lib.lsa_estimate(f, coeff, xdata[i])

plt.plot(xdata, ydata, xdata, yest)
plt.show()