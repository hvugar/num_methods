#ifndef EXAMPLE2_H
#define EXAMPLE2_H

#include <vector2d.h>
#include <matrix2d.h>
#include <printer.h>
#include <math.h>

//#define SAMPLE_1
//#define SAMPLE_2
#define SAMPLE_3

class Example2
{
public:
    Example2();

    void calculateLeft2RightSample();
    void calculateRight2LeftSample();
    void calculateLeft2Right(unsigned int N, unsigned int k, const DoubleMatrix &A, DoubleMatrix &beta, DoubleVector &qamma, DoubleVector &x);
    void calculateLeft2Right4(unsigned int N, unsigned int k, const DoubleMatrix &A, DoubleMatrix &beta, DoubleVector &qamma, DoubleVector &x, const DoubleVector &x0);
    void calculateRight2Left(unsigned int N, unsigned int k, const DoubleMatrix &A, DoubleMatrix &beta, DoubleVector &qamma, DoubleVector &x, const DoubleVector &x0);

    //
    double A(double t) const;
    double B(double t) const;
    double X(double t) const;

    void sample1();
    void sample2();

    void sample_n3();
    void sample_n4();
    void sample_n6();

    void sample_n4_test();

    static void Main(int argc, char* argv[]);
};

#endif // EXAMPLE2_H
