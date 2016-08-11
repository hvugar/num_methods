#ifndef PROBLEM1_H
#define PROBLEM1_H

#include <parabolicequation.h>
#include <cmethods.h>

class Problem1
{
public:
    Problem1();
    double v(unsigned int j) const;

    void calculate1();
    void calculate2();
    void calculate3();

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
