#ifndef TOMASMETHOD_H
#define TOMASMETHOD_H

#include "global.h"
#include "doublevector.h"

class TomasMethod
{
public:
    TomasMethod();
};

void MINIMUMSHARED_EXPORT TomasAlgorithm(const DoubleVector &a, const DoubleVector &b, const DoubleVector &c, const DoubleVector &d, DoubleVector &x);


#endif // TOMASMETHOD_H
