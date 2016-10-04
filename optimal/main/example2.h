#ifndef EXAMPLE2_H
#define EXAMPLE2_H

#include <vector2d.h>
#include <matrix2d.h>
#include <printer.h>

class Example2
{
public:
    Example2();

    void calculate();
    void calculate(unsigned int N, unsigned int k, const DoubleMatrix &a, DoubleMatrix &beta, DoubleVector &qamma, DoubleVector &x);
};

#endif // EXAMPLE2_H
