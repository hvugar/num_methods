#ifndef SYSTEM_LINEAR_ODE_NONLOCALCONTIONSM_H
#define SYSTEM_LINEAR_ODE_NONLOCALCONTIONSM_H

#include <load_sys/islodenlcsm.h>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <printer.h>

class SystemLinearODENonLocalContionsM : public ISystemLinearODENonLocalContionsM
{
public:
    static void Main(int agrc, char *argv[]);

    void initialize();
public:
    virtual double A(double t, unsigned int k, unsigned int row, unsigned int col) const;
    virtual double B(double t, unsigned int k, unsigned int row, unsigned int col) const;

    double X(double t, unsigned int k, unsigned int row, unsigned int col) const;
    double dX(double t, unsigned int k, unsigned int row, unsigned int col) const;
};

#endif // SYSTEM_LINEAR_ODE_NONLOCALCONTIONSM_H
