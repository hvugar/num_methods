#ifndef CAUCHYPROBLEMMEX_H
#define CAUCHYPROBLEMMEX_H

#include <grid/cauchyp.h>

#define SAMPLE_3

class CauchyProblemMEx : public CauchyProblemM1stOrder
{
public:
    static void Main(int agrc, char *argv[]);

    CauchyProblemMEx(const Dimension &grid);

protected:
    virtual double f(double x, const DoubleVector &y, unsigned int k, unsigned int i) const;

    double y(double x, unsigned int k, unsigned int i) const;


};

#endif // CAUCHYPROBLEMMEX_H
