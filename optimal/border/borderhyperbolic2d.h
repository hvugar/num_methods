#ifndef BORDERHYPERBOLIC2D_H
#define BORDERHYPERBOLIC2D_H

#include <pde_old/hyperbolicequation.h>

class MINIMUMSHARED_EXPORT BorderHyperbolic2D : public IHyperbolicEquation2D
{
public:
    virtual double initial1(unsigned int i, unsigned int j) const;
    virtual double initial2(unsigned int i, unsigned int j) const;
    virtual double boundary(unsigned int i, unsigned int j, unsigned int k) const;
    virtual double f(unsigned int i, unsigned int j, unsigned int k) const;

    double u(unsigned int i, unsigned int j, unsigned int k) const;

    double h1;
    double h2;
    double ht;
    unsigned int N1;
    unsigned int N2;
    unsigned int M;
    double a1;
    double a2;
    double lambda;

    void calculate(DoubleMatrix &u, double h1, double h2, double ht, unsigned int N1, unsigned int N2, unsigned int M, double a1, double a2) const;

    static void Main(int argc, char **argv);
};

class MINIMUMSHARED_EXPORT BorderHyperbolic2DN : public IHyperbolicEquation2D
{
public:
    virtual double initial1(unsigned int i, unsigned int j) const;
    virtual double initial2(unsigned int i, unsigned int j) const;
    virtual double boundary(unsigned int i, unsigned int j, unsigned int k) const;
    virtual double f(unsigned int i, unsigned int j, unsigned int k) const;

    double u(unsigned int i, unsigned int j, unsigned int k) const;
    double mu(unsigned int k, double t) const;
    double h(unsigned int i, unsigned int j, unsigned int k, double t, unsigned int n) const;

    double hx;
    double hy;
    double ht;
    unsigned int Nx;
    unsigned int Ny;
    unsigned int M;
    double a1;
    double a2;
    double lambda0;
    double lambda1;

    void calculateMVD3(DoubleMatrix &u, double hx, double hy, double ht, unsigned int Nx, unsigned int Ny, unsigned int M, double a1, double a2) const;

    static void Main(int argc, char **argv);
};

#endif // BORDERHYPERBOLIC2D_H
