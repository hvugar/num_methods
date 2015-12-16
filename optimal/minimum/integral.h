#ifndef INTEGRAL_H
#define INTEGRAL_H

#include "global.h"
#include "function.h"

double MINIMUMSHARED_EXPORT Trapesium(R1Function* f, unsigned int N, double a, double b);
double MINIMUMSHARED_EXPORT Trapesium(R2Function *f, unsigned int N1, unsigned int N2, double a1, double b1, double a2, double b2);
double MINIMUMSHARED_EXPORT Trapesium(R3Function *f, unsigned int N1, unsigned int N2, unsigned int N3, double a1, double b1, double a2, double b2, double a3, double b3);

#ifdef __cplusplus
extern "C" {
#endif

double trapesium(R1Function* f, unsigned int N, double a, double b);

#ifdef __cplusplus
}
#endif



#endif // INTEGRAL_H
