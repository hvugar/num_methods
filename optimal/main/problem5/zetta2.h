#ifndef ZETTA2_H
#define ZETTA2_H

#include <load_sys/islodenlcsm.h>

class Problem5Ex1;

class Zetta2 : public ISystemLinearODENonLocalContionsM
{
public:
    Zetta2(const ODEGrid &grid, const Problem5Ex1 *p5);

    virtual double A(double t UNUSED_PARAM, unsigned int k, unsigned int row, unsigned int col) const;
    virtual double B(double t UNUSED_PARAM, unsigned int k, unsigned int row) const;
private:
    const Problem5Ex1 *p5;
};

#endif // ZETTA2_H
