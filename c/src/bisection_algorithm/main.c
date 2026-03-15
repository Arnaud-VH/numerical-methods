#include <stdio.h>
#include <math.h>
#include "bisection.h"

double f(double x) {
   return (x * x) - 2;
}

int main() {
   double root;
   root = bisection(f, 1.0, 2.0, 1e-6);
   printf("Root found %f\n", root);
   
   return 0;
}