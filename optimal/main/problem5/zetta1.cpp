#include "zetta1.h"
#include "problem5ex1.h"

Zetta1::Zetta1(const ODEGrid &grid, const Problem5Ex1 *p5) : ISystemLinearODENonLocalContionsM(grid), p5(p5) {}

double Zetta1::A(double t, unsigned int k, unsigned int row, unsigned int col) const
{
    return p5->A(t, k, row, col);
}

double Zetta1::B(double t, unsigned int k, unsigned int row, unsigned int col) const
{
    return p5->C(t, k, 0, row, col);
}
