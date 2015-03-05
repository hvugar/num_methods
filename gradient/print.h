#ifndef _PRINT_H_
#define _PRINT_H_

#include <stdio.h>
#include <stdarg.h>
#include "minimum.h"

typedef void (*Printer)(RnFunction f, double *x, int n, ...);

void printer1(RnFunction f, double *x, int n, ...);
void printer2(RnFunction f, double *x, int n, ...);


#endif //_PRINT_H_