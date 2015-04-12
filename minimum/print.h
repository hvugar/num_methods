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

#endif //_PRINT_H_