#ifndef PROBLEM1_H
#define PROBLEM1_H

#include <parabolicequation.h>
#include <cmethods.h>

class Problem1 : IParabolicEquation
{
public:
    Problem1();

    double initial(unsigned int i) const;
    double boundary(Boundary type, unsigned int j) const;
    double f(unsigned int i, unsigned int j) const;
    double v(unsigned int j) const;

    void calculate();

private:
    double t0;
    double t1;
    double x0;
    double x1;
    double a;
    double lambda;
    double Te;
    DoubleMatrix u;
    double hx;
    double ht;
    unsigned int M;
    unsigned int N;
};

#endif // PROBLEM1_H
