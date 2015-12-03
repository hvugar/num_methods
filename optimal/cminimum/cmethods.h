#ifndef METHODS
#define METHODS

#include "cglobal.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

typedef double (*CR1Function)(double);
typedef double (*CRnFunction)(double*, unsigned int);

double derivative1(CR1Function fx, double x, double h);
double derivative2(CR1Function fx, double x, double h);
double derivative3(CR1Function fx, double x, double h);

void gradient1(CRnFunction fx, double *x, double *g, unsigned int n, double h);
void gradient2(CRnFunction fx, double *x, double *g, unsigned int n, double h);
void gradient3(CRnFunction fx, double *x, double *g, unsigned int n, double h);

double trapesium1(CR1Function fx, unsigned int n, double a, double b);
double trapesium2(CR1Function fx, double h, double a, double b);

void tomasAlgorithm(const double *a, const double *b, const double *c, const double *d, double *x, unsigned int n);

#endif // METHODS

