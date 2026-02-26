import numpy as np

def lsa(farray, xdata, ydata):
   """
   param farray: Represents the basis functions. Function that for each x returns a vector with the values of 
   each basis function at point x. 
   param xdata: Vector with the x components of the data. 
   param ydata: Vector with the y components of the data. 
   """
   n = len(xdata)
   m = len(ydata)
   b = len(farray(1))

   assert n == m, "Dimensions of x and y data do not match!"

   M = np.zeros([b,b])
   f = np.zeros([n,n])

   for i in range(n):
      vec = farray(xdata[i])
      for j in range(b):
         f[i,j] = vec[j]

   for i in range(b):
      for j in range(i,b):
         M[i,j] = np.dot(f[:,i],f[:,j])
         M[j,i] = M[i,j]
   
   rhs = np.array([np.dot(f[:,i], ydata) for i in range(b)])
   sol = np.linalg.solve(M,rhs)

   return sol

def lsa_estimate(farray, coeff, x):
   return np.dot(farray(x), coeff)
