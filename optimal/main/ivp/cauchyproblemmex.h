#ifndef CAUCHYPROBLEMMEX_H
#define CAUCHYPROBLEMMEX_H

#include <grid/cauchyp.h>

class CauchyProblemMEx : public CauchyProblemM
{
public:
    static void Main(int agrc, char *argv[]);

    CauchyProblemMEx(const Dimension &grid);

protected:
    virtual double f(double x, const DoubleVector &y, unsigned int k, unsigned int i) const;


};

#endif // CAUCHYPROBLEMMEX_H
