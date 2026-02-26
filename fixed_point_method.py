import math
import numpy as np

def g(x):
   return math.cos(x)

#Implementation of the fixed point method
def fixed_pm(g, x0, tol, maxiter):
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

def fixed_point_systems(g, x0, nmax = 100, tol = 10**(-6)):
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

def g2(x):
   n = len(x)
   out = [0]*n
   h = 0.1
   for i in range(1, n-1):
      out[i] = (((-h**2))/(1+x[i]**2) + x[i + 1] + x[i - 1])/2
   return out

###############################################################

def jacobi(A, b, x0, nmax = 100, tol = 10**(-6)):
   """
   Jacobi's iterative method for systems of equations. 

   param A: Coefficient matrix (n x n). 
   param b: Right-hand side vector of length n. 
   param x0: Initial estimate vector of length n.
   param nmax: Maximum number of iterations.
   param tol: Tolerance for the difference between two iterations before convergence. 
   """
   n = len(A)
   assert len(b) == n and len(x0) == n , "Dimensions do not match"

   sol = x0
   aux = np.add(sol, 2*tol)
   iterations = 0

   while ((iterations < nmax) and np.linalg.norm(sol - aux) > tol):
      aux = np.copy(sol)
      for i in range(n):
         sum = np.dot(A[i],aux) - A[i,i] * aux[i]
         sol[i] = (b[i] - sum) / A[i,i]
      iterations += 1
   
   if iterations >= nmax:
      print("Warning: Maximum number of iterations has been reached.")

   return sol


def guass_seidel(A, b, x0, nmax = 100, tol = 10**(-6)):
   """
   Solves Ax = b using Gauss-Seidel iterative method. 

   param A: Coefficient matrix (n x n). 
   param b: Right-hand side vector of length n. 
   param x0: Initial estimate vector of length n.
   param nmax: Maximum number of iterations.
   param tol: Tolerance for the difference between two iterations before convergence. 
   """
   n = len(A)
   assert len(b) == n and len(x0) == n, "Dimensions do not match"

   sol = x0.copy()
   aux = np.add(sol, 2*tol)

   iterations = 0

   while ((iterations < nmax) and np.linalg.norm(sol - aux) > tol):
      aux = np.copy(sol)
      for i in range(n):
         sum = 0
         for j in range(i):
            sum += A[i,j] * sol[j]
         for j in range(i + 1, n):
            sum += A[i,j] * aux[j]
         sol[i] = (b[i] - sum) / A[i,i]
      iterations += 1
   
   if iterations >= nmax:
      print("Warning: Maximum number of iterations has been reached.")

   return sol

def successive_over_relaxation(A, b, x0, w = 1.01, nmax = 1000, tol = 10**(-10)):
   """
   Solves the linear system Ax = b using Successive Over-Relaxation. 
   Method that improves on Gauss-Seidel with the relaxation parameter.
   
   param A: Coefficient matrix (n x n). 
   param b: Right-hand side vector of length n. 
   param x0: Initial estimate vector of length n.
   param w: Relaxation parameter.
   param nmax: Maximum number of iterations.
   param tol: Tolerance for the difference between two iterations before convergence. 
   """
   n = len(A)
   assert len(b) == n and len(x0) == n, "Dimensions do not match"

   sol = np.copy(x0)
   aux = np.add(sol, 2*tol)
   iterations = 0

   while (iterations < nmax and np.linalg.norm(sol - aux) < tol):
      aux = np.copy(sol)
      for i in range(n):
         sum = 0
         for j in range(i):
            sum += A[i,j] * sol[j]
         for j in range(i + 1, n):
            sum += A[i,j] * aux[j]
         sol[i] = (b[i] - sum) / A[i,i]
      iterations += 1

      sol = (1 - w) * aux + w * sol

   if iterations >= nmax:
      print("Warning: Maximum number of iterations has been reached.")

   return sol

n = 100
A = np.zeros((n,n))

#Create the diagonally dominant matrix
for i in range(n):
   A[i,i] = -2

for i in range(n-1):
   A[i,i + 1] = 1

for i in range(1,n):
   A[i,i - 1] = 1

b = np.array([1 for i in range(n)])
x0 = np.zeros(n)

print(successive_over_relaxation(A,b,x0))