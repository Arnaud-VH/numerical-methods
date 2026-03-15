#ifndef BISECTION_H
#define BISECTION_H

typedef double (*function)(double);

double bisection(function f, double a, double b, double tol);


#endif