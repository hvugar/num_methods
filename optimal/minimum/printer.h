#ifndef PRINTER_H
#define PRINTER_H

#include "doublevector.h"
#include "function.h"
#include "global.h"

void MINIMUMSHARED_EXPORT printMatrix1(const DoubleMatrix& x);

struct MINIMUMSHARED_EXPORT Printer
{
    virtual void print(unsigned int iteration, const DoubleVector& x, const DoubleVector &gradient, double alpha, RnFunction* fn) const = 0;
};

#endif // PRINTER_H
