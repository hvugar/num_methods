#ifndef TOMASMETHOD_H
#define TOMASMETHOD_H

#include "global.h"
#include "doublevector.h"

#ifdef __cplusplus
extern "C" {
#endif

void MINIMUMSHARED_EXPORT TomasAlgorithm(const DoubleVector &a, const DoubleVector &b, const DoubleVector &c, const DoubleVector &d, DoubleVector &x);
void MINIMUMSHARED_EXPORT tomasAlgorithm(const double *a, const double *b, const double *c, const double *d, double *x, unsigned int n);

#ifdef __cplusplus
}
#endif


#endif // TOMASMETHOD_H
