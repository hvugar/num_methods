#include "zetta0.h"
#include "problem4ex1.h"

//Zetta0::Zetta0(const Problem4Ex1 &p4) : ISystemLinearODENonLocalContionsV(), p4(p4) {}

//double Zetta0::A(double t, unsigned int k, unsigned int row, unsigned int col) const
//{
//    return p4.A(t,k,row,col);
//}

//double Zetta0::B(double t, unsigned int k, unsigned int row) const
//{
//    return p4.B(t, k, row);
//}

Zetta01::Zetta01(const Problem4Ex1 &p) : p(p)
{}

unsigned int Zetta01::equationsNumber() const
{
    return 3;
}

double Zetta01::A(const GridNodeODE &node, unsigned int row, unsigned int col) const
{
    return p.A(node,row,col);
}

double Zetta01::B(const GridNodeODE &node, unsigned int row) const
{
    return p.B(node,row);
}
