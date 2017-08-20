#include "zettai.h"
#include "problem4ex1.h"

Zettai::Zettai(const Problem4Ex1 &p4, unsigned int i)
    : ISystemLinearODENonLocalContionsM(), p4(p4), i(i)
{}

double Zettai::A(double t, unsigned int k, unsigned int row, unsigned int col) const
{
    return p4.A(t, k, row, col);
}

double Zettai::B(double t, unsigned int k, unsigned int row, unsigned int col) const
{
    return p4.B(t, k, i, row, col);
}
