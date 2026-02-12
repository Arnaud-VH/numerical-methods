import math

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