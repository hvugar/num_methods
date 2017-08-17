#include "zetta2.h"
#include "problem5ex1.h"

Zetta2::Zetta2(const ODEGrid &grid, const Problem5Ex1 *p5) : ISystemLinearODENonLocalContionsM(grid), p5(p5) {}

double Zetta2::A(double t, unsigned int k, unsigned int row, unsigned int col) const
{
    return p5->A(t, k, row, col);
}

double Zetta2::B(double t, unsigned int k, unsigned int row, unsigned int col) const
{
    return p5->C(t, k, 1, row, col);
}
