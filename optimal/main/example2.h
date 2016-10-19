#ifndef EXAMPLE2_H
#define EXAMPLE2_H

#include <vector2d.h>
#include <matrix2d.h>
#include <printer.h>

class Example2
{
public:
    Example2();

    void init1();
    void init2();
    void init3();
    void init4();
    void init5();
    void calculateLeft2Right(unsigned int N, unsigned int k, const DoubleMatrix &a, DoubleMatrix &beta, DoubleVector &qamma, DoubleVector &x, const DoubleVector &x0);
    void calculateRight2Left(unsigned int N, unsigned int k, const DoubleMatrix &a, DoubleMatrix &beta, DoubleVector &qamma, DoubleVector &x, const DoubleVector &x0);

    static void Main(int argc, char* argv[]);
};

#endif // EXAMPLE2_H
