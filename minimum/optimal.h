#ifndef __OPTIMAL_H
#define __OPTIMAL_H

#include "methods.h"
#include "print.h"
#define M_E1 2.7182818284590452353602874713527

void smp2_control();
void smp2_RungaKuttaSystem1(double t0, double *x0, double t, double *x, int n, double h, double u);
void smp2_RungaKuttaSystem2(double t0, double *p0, double t, double *p, int n, double h, double *x, double u);
double smp2_f0(double t, double *x, int n, double u);
double smp2_f1(double t, double *x, int n, double u);
double smp2_f2(double t, double *x, int n, double u);
double smp2_Hamilton(double t, double *x, double *psi, int n, double u);
double smp2_dPsi1(double t, double *x, double *psi, int n, double u);
double smp2_dPsi2(double t, double *x, double *psi, int n, double u);
double smp2_dU(double t, double *x, double *psi, int n, double u);
double smp2_JSum(double *t, double *x1, double *x2, int n, double *u, int N);

#endif // __OPTIMAL_H
