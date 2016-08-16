#ifndef PROBLEM2_H
#define PROBLEM2_H

#include "newtonheatprocess.h"
#include <cmethods.h>

class Problem2 : public NewtonHeatProcess
{
public:
    Problem2();
    virtual ~Problem2();

    virtual double vm(unsigned int j) const;
    virtual double vl(unsigned int j) const;
    virtual double vr(unsigned int j) const;

    virtual double initial(unsigned int i) const;

    void calculate1(DoubleMatrix &m, double ht, double hx, unsigned int M, unsigned int N, double lambdaM, double lambdaL, double lambdaR, double a);

public:
    double t0;
    double t1;
    double x0;
    double x1;
    double hx;
    double ht;
    unsigned int M;
    unsigned int N;
    double a;

    double lambdaM;
    double lambdaL;
    double lambdaR;

    unsigned int L;
    DoubleVector k;
    DoubleVector xi;
    DoubleVector z;

    DoubleMatrix *pm;
};

#endif // PROBLEM2_H
