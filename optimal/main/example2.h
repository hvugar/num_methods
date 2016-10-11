#ifndef EXAMPLE2_H
#define EXAMPLE2_H

#include <vector2d.h>
#include <matrix2d.h>
#include <printer.h>

class Example2
{
public:
    Example2();

    void init();
    void calculate();
    void calculate1(unsigned int N, unsigned int k, const DoubleMatrix &a, DoubleMatrix &beta, DoubleVector &qamma, DoubleVector &x);
    void calculate2(unsigned int N, unsigned int k, const DoubleMatrix &a, DoubleMatrix &beta, DoubleVector &qamma, DoubleVector &x);

    static void Main(int argc, char* argv[]);
};

#endif // EXAMPLE2_H
