#ifndef SYSTEM_LINEAR_ODE_NONLOCALCONTIONSV_H
#define SYSTEM_LINEAR_ODE_NONLOCALCONTIONSV_H

#define SAMPLE_3

#include <load_sys/islodenlcsv.h>

class SystemLinearODENonLocalContionsV : public ISystemLinearODENonLocalContionsV
{
public:
    SystemLinearODENonLocalContionsV(const ODEGrid &grid);
    static void Main(int agrc, char *argv[]);

    void initialize();
    double x(double t, int i) const;

public:
    virtual double A(double t UNUSED_PARAM, unsigned int k, unsigned int row, unsigned int col) const;
    virtual double B(double t UNUSED_PARAM, unsigned int k, unsigned int row) const;
};

#endif // SYSTEM_LINEAR_ODE_NONLOCALCONTIONSV_H
