#ifndef PRINTER_H
#define PRINTER_H

#include "doublevector.h"
#include "function.h"
#include "global.h"

struct MINIMUMSHARED_EXPORT Printer
{
    virtual void print(unsigned int iteration, const DoubleVector& x, const DoubleVector &gradient, double alpha, RnFunction* fn) const = 0;

    static void printMatrix(const DoubleMatrix &x, unsigned int m = 10, unsigned int n = 10, const char* s = NULL, FILE* f = stdout);
    static void printVector(const DoubleVector &x, unsigned int n = 10, const char* s=NULL, FILE* f = stdout);
    static void printVector(const DoubleVector &x, const char *s = NULL, unsigned int n = 10, unsigned int start = 0, unsigned int end = 0, FILE *file = stdout);
};

#endif // PRINTER_H
