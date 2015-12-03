#ifndef METHODS
#define METHODS

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

typedef double (*R1Function)(double);
typedef double (*RnFunction)(double*, unsigned int);

double derivative1(R1Function fx, double x, double h);
double derivative2(R1Function fx, double x, double h);
double derivative3(R1Function fx, double x, double h);

void gradient1(RnFunction fx, double *x, double *g, unsigned int n, double h);
void gradient2(RnFunction fx, double *x, double *g, unsigned int n, double h);
void gradient3(RnFunction fx, double *x, double *g, unsigned int n, double h);

double trapesium1(R1Function fx, unsigned int n, double a, double b);
double trapesium2(R1Function fx, double h, double a, double b);

void tomasAlgorithm(const double *a, const double *b, const double *c, const double *d, double *x, unsigned int n);

#endif // METHODS

