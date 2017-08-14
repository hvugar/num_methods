#ifndef NUMINTEGRALEXP1_H
#define NUMINTEGRALEXP1_H

#include <grid/integral1.h>
#include <stdio.h>

class NumIntegralExp1 : public NumericalIntegral
{
public:
    static void Main(int agrc, char *argv[]);

    NumIntegralExp1(const ODEGrid &grid);

    virtual double f(double x, int n) const;
};

#endif // NUMINTEGRALEXP1_H
