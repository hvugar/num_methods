#ifndef CAUCHYPROBLEMNONLOCALCONTIONS_H
#define CAUCHYPROBLEMNONLOCALCONTIONS_H

#include <vector2d.h>
#include <matrix2d.h>
#include <vector>
#include <grid/cauchyp.h>

using namespace std;

class CauchyProblemNonLocalContions : public CauchyProblemM
{
public:
    static void Main(int agrc, char *argv[]);

    CauchyProblemNonLocalContions(const Dimension &grid);

    DoubleVector times;
    unsigned int n = 2;
    unsigned int L = 3;

    double A(unsigned int k, unsigned int i, unsigned int j) const;
    double B(unsigned int k, unsigned int i) const;

    std::vector<DoubleMatrix> alpha;
    DoubleVector betta;

    double x1(unsigned int k) const;
    double x2(unsigned int k) const;

    void calculate1();

    double S0(unsigned int i, double t, const DoubleVector &x, unsigned int k) const;

protected:
    virtual double f(double x, const DoubleVector &y, unsigned int k, unsigned int i) const;
};

#endif // CAUCHYPROBLEMNONLOCALCONTIONS_H
