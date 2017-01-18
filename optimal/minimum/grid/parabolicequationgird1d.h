#ifndef PARABOLICEQUATIONGIRD1D_H
#define PARABOLICEQUATIONGIRD1D_H

#include "gridmethod.h"

class MINIMUMSHARED_EXPORT ParabolicEquationGird1D : public GridMethod
{
public:
    ParabolicEquationGird1D();
    virtual ~ParabolicEquationGird1D() {}

protected:
    virtual double initial(unsigned int n) const = 0;
    virtual double boundary(unsigned int m, Boundary boundary) const = 0;
    virtual double f(unsigned int n, unsigned int m) const = 0;
    virtual double a(unsigned int n, unsigned int m) const = 0;

public:
    void gridMethod(DoubleMatrix &u);

    /* dirichlet conditions */
    void calculateN4L2RD(DoubleMatrix &u);
    void calculateN4R2LD(DoubleMatrix &u);

    void calculateN6L2RD(DoubleMatrix &u);
    void calculateN6R2LD(DoubleMatrix &u);
};

#endif // PARABOLICEQUATIONGIRD1D_H
