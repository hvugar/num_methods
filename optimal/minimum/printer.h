#ifndef PRINTER_H
#define PRINTER_H

#include "doublevector.h"
#include "function.h"
#include "global.h"

struct MINIMUMSHARED_EXPORT Printer
{
    virtual void print(unsigned int iteration, const DoubleVector& x, const DoubleVector &gradient, double alpha, RnFunction* fn) const = 0;

    static void printMatrix(const DoubleMatrix& x, int m, int n, FILE* f = stdout, const char* s=NULL);
    static void printVector(const DoubleVector& x, int n, FILE* f = stdout, const char* s=NULL);
};

#endif // PRINTER_H
