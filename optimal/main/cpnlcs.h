#ifndef CAUCHYPROBLEMNONLOCALCONTIONS_H
#define CAUCHYPROBLEMNONLOCALCONTIONS_H

#include <vector2d.h>
#include <matrix2d.h>

class CauchyProblemNonLocalContions
{
public:
    CauchyProblemNonLocalContions();

    DoubleVector times;
    unsigned int n = 2;
    double h = 0.1;
    unsigned int N = 10;
    unsigned int L = 3;

    double A(unsigned int k, unsigned int i, unsigned int j) const;
    double B(unsigned int k, unsigned int i) const;

    DoubleMatrix alpa0;
    DoubleMatrix alpa1;
    DoubleMatrix alpa2;
    DoubleVector betta;

    double x1(unsigned int k) const;
    double x2(unsigned int k) const;
};

#endif // CAUCHYPROBLEMNONLOCALCONTIONS_H
