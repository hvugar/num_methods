#ifndef _PRINT_H_
#define _PRINT_H_

#include <stdio.h>
#include <stdarg.h>

#include "methods.h"

typedef double (*R1Function)(double);
typedef double (*RnFunction)(double*, int);

typedef void (*Printer)(RnFunction f, double *x, int n, ...);

void printer1(RnFunction f, double *x, int n, ...);
void printer2(RnFunction f, double *x, int n, ...);
void printer3(RnFunction f, double *x, int n, ...);
void printer4(RnFunction f, double *x, int n, ...);

void printX(char *label, double *x, int n);
void _print1(char *s, double *a, int n);
void _seperator();

void _printM(double **x, int m, int n);
void _printV(double *x, int m, int n);

#endif //_PRINT_H_