#ifndef __OPTIMAL_H
#define __OPTIMAL_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "methods.h"
#include "print.h"

#define M_E1 2.7182818284590452353602874713527
#define dx 0.000000001

typedef struct 
{
    double *t;
    double *u;
    double *x1;
    double *x2;
    double *psi1;
    double *psi2;
    double *gradJ;

    double t1;
    double t2;
    double h;
    double n;

    double x01;
    double x02;

    double *s;
} Process;

double fx0(double t, double x1, double x2, double u);
double T(double t, double x1, double x2, double u);
double fx1(double t, double x1, double x2, double u);
double fx2(double t, double x1, double x2, double u);
double H(double t, double x1, double x2, double u, double psi1, double psi2);
double fp1(double t, double x1, double x2, double psi1, double psi2, double u);
double fp2(double t, double x1, double x2, double psi1, double psi2, double u);
double gradJ(double t, double x1, double x2, double psi1, double psi2, double u);
double JSum(Process *p);
void init_process(Process *p);
void free_process(Process *p);
void calculate_x(Process *p);
void calculate_psi(Process *p);
void calculate_gradient(Process *p);

#endif // __OPTIMAL_H
