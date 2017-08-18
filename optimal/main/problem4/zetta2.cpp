#include "zetta2.h"
#include "problem4ex1.h"

Zetta2::Zetta2(const ODEGrid &grid, const Problem4Ex1 *p5) : ISystemLinearODENonLocalContionsM(grid), p4(p5) {}

double Zetta2::A(double t, unsigned int k, unsigned int row, unsigned int col) const
{
    return p4->A(t, k, row, col);
}

double Zetta2::B(double t, unsigned int k, unsigned int row, unsigned int col) const
{
    return p4->B(t, k, 1, row, col);
}
