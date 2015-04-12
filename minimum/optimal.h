#ifndef __OPTIMAL_H
#define __OPTIMAL_H

#include "methods.h"
#include "print.h"
#define M_E1 2.7182818284590452353602874713527

void sample1();
void sample2();

double _f0(double t, double *x, int n, double u);
double _f1(double t, double *x, int n, double u);
double _f1(double t, double *x, int n, double u);
double _p1(double t, double *x, int n, double *p);
double _p2(double t, double *x, int n, double *p);
double _H1(double t, double *x, int n, double u, double *p);

#endif // __OPTIMAL_H
