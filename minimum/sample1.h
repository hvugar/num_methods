#ifndef _SAMPLE1_H_
#define _SAMPLE1_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "methods.h"
#include "print.h"

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
} Process1;

double fx0(double t, double x1, double x2, double u);
double T(double t, double x1, double x2, double u);
double fx1(double t, double x1, double x2, double u);
double fx2(double t, double x1, double x2, double u);
double H(double t, double x1, double x2, double u, double psi1, double psi2);
double fp1(double t, double x1, double x2, double psi1, double psi2, double u);
double fp2(double t, double x1, double x2, double psi1, double psi2, double u);
double gradJ(double t, double x1, double x2, double psi1, double psi2, double u);
double JSum(Process1 *p);
void init_process(Process1 *p);
void free_process(Process1 *p);
void calculate_x(Process1 *p);
void calculate_psi(Process1 *p);
void calculate_gradient(Process1 *p);

#endif // _SAMPLE1_H_
