#ifndef NEWTONHEATPROCESS_H
#define NEWTONHEATPROCESS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <global.h>
#include <matrix2d.h>

class NewtonHeatProcess
{
public:
    virtual double vm(unsigned int j) const = 0;
    virtual double vl(unsigned int j) const = 0;
    virtual double vr(unsigned int j) const = 0;

    virtual double initial(unsigned int i) const = 0;

    void calculate(DoubleMatrix &m, double ht, double hx, unsigned int M, unsigned int N, double lambdaM, double lambdaL, double lambdaR, double a = 1.0);
};



#endif // NEWTONHEATPROCESS_H
