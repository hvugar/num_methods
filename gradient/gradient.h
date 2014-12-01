#ifndef __GRADIENT_H
#define __GRADIENT_H

#include "minimum.h"

#ifdef __cplusplus
extern "C" {
#endif

void gradient(RnFunction f, double *x, int n, double dx, double *gradients);
double grad_module(double *grads, int n);

#ifdef __cplusplus
}
#endif


#endif // __GRADIENT_H