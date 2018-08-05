#ifndef PROBLEM2HNM_H
#define PROBLEM2HNM_H

#include "common.h"

class PROBLEM2HSHARED_EXPORT Problem2HNM
{
protected:
    virtual double initial1(const SpaceNodePDE &sn) const;
    virtual double initial2(const SpaceNodePDE &sn) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
public:
    double a;
    double lambda;

    void layerInfo(DoubleMatrix &u);
    void layerInfo(DoubleVector &u);

    void solveEquation2D1(const Dimension &time, const Dimension &dimx, const Dimension &dimy, DoubleMatrix &u);
    void solveEquation2D2(const Dimension &time, const Dimension &dimx, const Dimension &dimy, DoubleMatrix &u);
    void solveEquation2D3(const Dimension &time, const Dimension &dimx, const Dimension &dimy, DoubleMatrix &u);

    void solveEquation1D1(const Dimension &time, const Dimension &dimx, DoubleVector &u);

    static void Main(int argc, char* agr[]);
};

#endif // PROBLEM2HNM_H
