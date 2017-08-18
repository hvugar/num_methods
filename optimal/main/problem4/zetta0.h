#ifndef ZETTA0_H
#define ZETTA0_H

#include <load_sys/islodenlcsv.h>

class Problem4Ex1;

class Zetta0 : public ISystemLinearODENonLocalContionsV
{
public:
    Zetta0(const ODEGrid &grid, const Problem4Ex1 *p5);

    virtual double A(double t, unsigned int k, unsigned int row, unsigned int col) const;
    virtual double B(double t, unsigned int k, unsigned int row) const;
private:
    const Problem4Ex1 *p5;
};

#endif // ZETTA0_H
