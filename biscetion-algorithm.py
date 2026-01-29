import math
from scipy import optimize as opt

def f(r):
   M = 30000; n = 10; v= 900
   return M-v*(1+r)*((1+r)**n-1)/r

def f2(x):
   return x**2 - 2

#Should verify that b > a, the tolerance is positive
#Verify that f(a)*f(b) is negative.
def bisection(f, a, b, tol = 0.5*10**(-6)):
   left = a; right = b
   k = math.ceil(math.log((b-a)/tol)/math.log(2)-1)
   for _ in range(k):
      mid = (left + right) / 2
      if (f(left)*f(mid) < 0):
         #Solution on the left 
         right = mid
      else:
         #Solution on the right
         left = mid

   return mid

ir = bisection(f, 0.001, 0.5, 10**-12)
print(ir)

#Built in bisection method from scipy
ir= opt.bisect(f, 0.001, 1)
print(ir)