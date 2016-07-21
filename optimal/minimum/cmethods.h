#ifndef METHODS
#define METHODS

#include <global.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef double (*CR1Function)(double);
typedef double (*CRnFunction)(double*, unsigned int);

typedef double (*ODE1stOrderEquation) (double x, double y);
typedef double (*ODE1stOrderEquationN) (double x, double *y, size_t n);

double derivative(CR1Function fx);
double derivative1(CR1Function fx, double x, double h);
double derivative2(CR1Function fx, double x, double h);
double derivative3(CR1Function fx, double x, double h);

void gradient1(CRnFunction fx, double *x, double *g, unsigned int n, double h);
void gradient2(CRnFunction fx, double *x, double *g, unsigned int n, double h);
void gradient3(CRnFunction fx, double *x, double *g, unsigned int n, double h);

double trapesium1(CR1Function fx, unsigned int n, double a, double b);
double trapesium2(CR1Function fx, double h, double a, double b);

void MINIMUMSHARED_EXPORT tomasAlgorithm(const double *a, const double *b, const double *c, const double *d, double *x, unsigned int n);
void gaussianElimination(double **a, double *b, double *x, unsigned int n);
void gaussJordanElimination(double **a, double *b, double *x, unsigned int n);


void MINIMUMSHARED_EXPORT euler(double x0, double y0, double xN, double yN, unsigned int N, double *x, double *y, ODE1stOrderEquation eq);
void eulerMod(double x0, double y0, double xN, double yN, unsigned int N, double *x, double *y, ODE1stOrderEquation eq);
void runge_kutta_rk3(double x0, double y0, double xN, double yN, unsigned int N, double *x, double *y, ODE1stOrderEquation eq);
void runge_kutta_rk4(double x0, double y0, double xN, double yN, unsigned int N, double *x, double *y, ODE1stOrderEquation eq);

void MINIMUMSHARED_EXPORT runge_kutta_rk4_system(double x0, double x1, double *y0, double **y1, size_t n, unsigned int N, double h, ODE1stOrderEquationN *eq);
//void RungaKuttaSystem(RmFunction *f, double x0, const double *y0, double x, double *y, const int n, double h);

#ifdef __cplusplus
}
#endif


#endif // METHODS

