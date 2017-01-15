#ifndef PARABOLICEQUATIOND_H
#define PARABOLICEQUATIOND_H

#include <vector2d.h>
#include <matrix2d.h>
#include <cmethods.h>

class ParabolicEquationD
{
public:
    virtual double initial(unsigned int n, double hx) const = 0;
    virtual double boundary(Boundary type, unsigned int m, double ht) const = 0;
    virtual double f(unsigned int n, unsigned int m, double hx, double ht) const = 0;

    void calculateU(DoubleMatrix &u, double hx, double ht, unsigned int N, unsigned int M, double a=1.0) const;
};

#endif // PARABOLICEQUATIOND_H
