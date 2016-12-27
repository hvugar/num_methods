#ifndef BORDERPARABOLIC2D_H
#define BORDERPARABOLIC2D_H

#include <global.h>
#include <parabolicequation.h>

class MINIMUMSHARED_EXPORT BorderParabolic2D : public IParabolicEquation2D
{
public:
    virtual double initial(unsigned int i, unsigned int j) const;
    virtual double boundary(unsigned int i, unsigned int j, unsigned int k) const;
    virtual double f(unsigned int i, unsigned int j, unsigned int k) const;

    void calculate(DoubleMatrix &u, double h1, double h2, double ht, unsigned int N1, unsigned int N2, unsigned int M, double a1, double a2) const;

    double u(unsigned int i, unsigned int j, unsigned int k) const;

    double h1 = 0.01;
    double h2 = 0.01;
    double ht = 0.01;
    unsigned int N1 = 100;
    unsigned int N2 = 100;
    unsigned int M  = 100;
    double a1 = 1.0;
    double a2 = 1.0;

    static void Main(int argc, char **argv);
};

#endif // BORDERPARABOLIC2D_H
