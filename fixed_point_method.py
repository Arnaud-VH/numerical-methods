import math
import numpy as np

def g(x):
   return math.cos(x)

#Implementation of the fixed point method
def fixedPM(g, x0, tol, maxiter):
   sol = x0
   i = 0
   oldsol = x0 + 2*tol
   while (i < maxiter and abs(sol - oldsol) > tol):
      oldsol = sol
      sol = g(sol)
      i += 1

   if i >= maxiter:
      print("HELP! We have passed the maximum number of iterations.") 
   print(f"The algorithm went through {i} iterations.")
   return sol

solution = fixedPM(g, 1, 1e-12, 100)
print(solution)

def fixedPointSystems(g, x0, nmax = 100, tol = 10**(-6)):
   """
   Fixed point method for systems of equations in the form g(x) = x

   :param g: Function that should return a vector (np.array) dimension compatible with x0
   :param x0: Initial approximation (np.array)
   :param nmax: Maximum numberof iterations (int)
   :param tol: Target difference between consecutive iterations (float)
   """
   sol = x0
   aux = 0
   iter = 0
   
   while ((iter < nmax) and np.linalg.norm(sol - aux, np.Inf) > tol) or iter == 0:
      aux = np.copy(sol)
      sol = g(sol)
      iter += 1
   
   if iter >= nmax:
      print("Warning: max number of iterations was reached.")
   return sol

def g1(x):
   out = np.array([-math.cos(x[0] + x[1]) / 3, math.sin(x[0] - x[1]) / 4])
   return out

print(fixedPointSystems(g1 , np.array([0,0])))

def g2(x):
   n = len(x)
   out = [0]*n
   h = 0.1
   for i in range(1, n-1):
      out[i] = (((-h**2))/(1+x[i]**2) + x[i + 1] + x[i - 1])/2
   return out

print(fixedPointSystems(g2,np.array([0.05 for i in range(50)])))