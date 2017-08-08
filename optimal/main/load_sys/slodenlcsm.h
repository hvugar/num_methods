#ifndef SYSTEMLINEARODENONLOCALCONTIONSM_H
#define SYSTEMLINEARODENONLOCALCONTIONSM_H

#include <load_sys/islodenlcsm.h>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


class SystemLinearODENonLocalContionsM : public ISystemLinearODENonLocalContionsM
{
public:
    static void Main(int agrc, char *argv[]);

    SystemLinearODENonLocalContionsM(const ODEGrid& grid);

protected:
    virtual double A(double t UNUSED_PARAM, unsigned int k, unsigned int row, unsigned int col) const = 0;
    virtual double B(double t UNUSED_PARAM, unsigned int k, unsigned int row, unsigned int col) const = 0;

    double X(double t, unsigned int k, unsigned int row, unsigned int col) const;
    double dX(double t, unsigned int k, unsigned int row, unsigned int col) const;
};

#endif // SYSTEMLINEARODENONLOCALCONTIONSM_H
