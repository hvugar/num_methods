#include "zettai.h"
#include "problem4ex1.h"

Zettai::Zettai(const Problem4Ex1 &p4, unsigned int i)
    : ISystemLinearODENonLocalContionsM(), p4(p4), i(i)
{}

double Zettai::A(const GridNodeODE &node, unsigned int row, unsigned int col) const
{
    return p4.A(node, row, col);
}

double Zettai::B(const GridNodeODE &node, unsigned int row, unsigned int col) const
{
    return p4.C(node, i, row, col);
}

////////////////////////////////////////////////////////////////////////////////////////////////

Zettai1::Zettai1(const Problem4Ex1 &p, unsigned int i) : p(p), i(i)
{}

double Zettai1::A(const GridNodeODE &node, unsigned int row, unsigned int col) const
{
    return p.A(node, row, col);
}

double Zettai1::B(const GridNodeODE &node, unsigned int row) const
{
    return p.C(node, i, row, cur_col);
}

double Zettai1::C(const GridNodeODE &node, unsigned int row, unsigned int col) const
{
    return p.C(node, i, row, col);
}

unsigned int Zettai1::equationsNumber() const
{
    return 3;
}

void Zettai1::calculateM(const std::vector<LinearODE1stOrder::Condition> &cs, const DoubleMatrix &betta, std::vector< std::vector<DoubleVector> > &zmi)
{
    unsigned int n = equationsNumber();

    zmi.resize(n);
    for (unsigned int row=0; row<n; row++)
        zmi[row].resize(n);

    for (unsigned int col=0; col<n; col++)
    {
        cur_col = col;

        std::vector<DoubleVector> z;
        calculate(cs, betta.col(cur_col), z);

        for (unsigned int row=0; row<n; row++)
        {
            zmi[row][cur_col] = z[row];
            z[row].clear();
        }
        z.clear();
    }
}
