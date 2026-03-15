#include "bisection.h"
#include <math.h>

double bisection(function f, double a, double b, double tol) {
   double left = a; 
   double right = b;
   double mid;
   int k = ceil(log((b - a) / tol) / log(2) - 1);
   for (int i = 0; i <= k; i++) {
      mid = (left + right) / 2;
      if (f(left) * f(mid) < 0) {
         right = mid;
      } else {
         left = mid;
      }
   } 
   return mid;
}