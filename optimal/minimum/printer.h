#ifndef PRINTER_H
#define PRINTER_H

#include "doublevector.h"
#include "function.h"

struct Printer
{
    virtual void print(unsigned int iteration, const DoubleVector& x, const DoubleVector &gradient, double alpha, RnFunction* fn) const = 0;
};

#endif // PRINTER_H
