import math 
import numpy as np
import matplotlib.pyplot as plt

def f(x):
   return math.exp(x)

def df(x):
   return math.exp(x)

n = 30
#Create a list of n equally spaced points.
x = np.linspace(0,1,n)

h = x[1]-x[0]
y = [f(u) for u in x]

dy = np.zeros((n))

for i in range(1,n-1):
   #Central difference
   dy[i] = (y[i+1] - y[i-1])/(2*h)

#Forward difference
dy[0] = (-3*y[0] + 4*y[1] - y[2]) / (2*h)
dy[-1] = (3*y[-1] - 4*y[-2] + y[-3]) / (2*h)

plt.plot(x,dy)
plt.show()