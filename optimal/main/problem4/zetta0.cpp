#include "zetta0.h"
#include "problem4ex1.h"

Zetta0::Zetta0(const ODEGrid &grid, const Problem4Ex1 *p5)
    : ISystemLinearODENonLocalContionsV(grid), p5(p5) {}

double Zetta0::A(double t, unsigned int k, unsigned int row, unsigned int col) const
{
    return p5->A(t,k,row,col);
}

double Zetta0::B(double t, unsigned int k, unsigned int row) const
{
    return p5->C(t, k, row);
}

