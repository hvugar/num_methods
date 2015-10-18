#ifndef INTEGRAL_H
#define INTEGRAL_H

#include "global.h"
#include "function.h"

//#ifdef __cplusplus
//extern "C" {
//#endif
double MINIMUMSHARED_EXPORT Trapesium(R1Function* f, unsigned int N, double a, double b);
double MINIMUMSHARED_EXPORT Trapesium(R2Function *f, unsigned int N1, unsigned int N2, double a1, double b1, double a2, double b2);
double MINIMUMSHARED_EXPORT Trapesium(R3Function *f, unsigned int N1, unsigned int N2, unsigned int N3, double a1, double b1, double a2, double b2, double a3, double b3);
//#ifdef __cplusplus
//}
//#endif

class MINIMUMSHARED_EXPORT Integral
{
public:
    Integral();
    ~Integral();
};

#endif // INTEGRAL_H
