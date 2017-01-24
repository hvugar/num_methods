#ifndef PARABOLICIBVP1_H
#define PARABOLICIBVP1_H

#include <grid/pibvp.h>

class ParabolicIBVP1 : public ParabolicIBVP
{
public:
    inline double U(unsigned int n, unsigned int m) const;

protected:
    virtual double initial(unsigned int n) const;
    virtual double boundary(unsigned int m, BoundaryType boundary) const;
    virtual double f(unsigned int n, unsigned int m) const;
    virtual double a(unsigned int n, unsigned int m) const;
};

#endif // PARABOLICIBVP1_H
