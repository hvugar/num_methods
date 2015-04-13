#ifndef __OPTIMAL_H
#define __OPTIMAL_H

#include "methods.h"
#include "print.h"
#define M_E1 2.7182818284590452353602874713527

void control1();
void sample2();

double f0(double t, double *x, int n, double u);
double f1(double t, double *x, int n, double u);
double f2(double t, double *x, int n, double u);
double Hamilton(double t, double *x, int n, double u, double *p);

double dPhi1(double t, double *x, int n, double *p, double u);
double dPhi2(double t, double *x, int n, double *p, double u);
double dU(double t, double *x, int n, double u, double *p);

double JSum(double *t, double *x1, double *x2, double n, double *u, int N);
double Integral(double *u, int N);

void _RungaKuttaSystem1(double x0, double *y0, double x, double *y, int n, double h, double u);
void _RungaKuttaSystem2(double x0, double *y0, double x1, double *y1, int n, double h, double *x, double u);

double _p1(double t, double *x, int n, double *p, double u);
double _p2(double t, double *x, int n, double *p, double u);


#endif // __OPTIMAL_H
