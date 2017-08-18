#ifndef ZETTA1_H
#define ZETTA1_H

#include <load_sys/islodenlcsm.h>

class Problem4Ex1;

class Zetta1 : public ISystemLinearODENonLocalContionsM
{
public:
    Zetta1(const ODEGrid &grid, const Problem4Ex1 *p4);

    virtual double A(double t, unsigned int k, unsigned int row, unsigned int col) const;
    virtual double B(double t, unsigned int k, unsigned int row, unsigned int col) const;
private:
    const Problem4Ex1 *p4;
};

#endif // ZETTA1_H
