#ifndef PROBLEM1KZX_H
#define PROBLEM1KZX_H

#include <cmethods.h>
#include <function.h>
#include <printer.h>
#include <gradient_cjt.h>
#include <vector>
#include <rungekutta.h>
#include <example3.h>

class Problem1KZX : public RnFunction, public IGradient, public IPrinter, public ConjugateGradient
{
public:
    Problem1KZX();
    virtual ~Problem1KZX() {}

    virtual double fx(const DoubleVector &x);
    virtual void gradient(const DoubleVector &x, DoubleVector &g);
    virtual void print(unsigned int iteration, const DoubleVector &x, const DoubleVector &gradient, double alpha, RnFunction *fn) const;
    virtual double mu(unsigned int i UNUSED_PARAM) const { return 1.0; }

    static void Main(int argc, char* argv[]);

    void calculateU(DoubleMatrix &u, double ht, double hx, unsigned int M, unsigned int N, double alpha, double lambda0, double lambdal, double a);
    // Duz xetler
    void calculateU1(DoubleMatrix &u, double ht, double hx, unsigned int M, unsigned int N, double alpha, double lambda0, double lambdal, double a);
    // Qaus
    void calculateU2(DoubleMatrix &u, double ht, double hx, unsigned int M, unsigned int N, double alpha, double lambda0, double lambdal, double a);
    // Qovma
    void calculateU3(DoubleMatrix &u, double ht, double hx, unsigned int M, unsigned int N, double alpha, double lambda0, double lambdal, double a);

    // Qovma
    void calculateP3(DoubleMatrix &m, double ht, double hx, unsigned int M, unsigned int N, double alpha, double lambda0, double lambdal, double a);

    virtual double initial(unsigned int i) const;
    virtual double vm(unsigned int j) const;
    virtual double vr(unsigned int j) const;
    virtual double vl(unsigned int j) const;
    virtual double vs(double j) const;

private:
    double hx = 0.01;
    double ht = 0.01;
    unsigned int M = 100;
    unsigned int N = 100;
    double a = 1.0;

    double alpha = 1.0;
    double lambda0 = 1.0;
    double lambdal = 1.0;

    unsigned int L = 2;
    double alpha1 = 1.0;
    double alpha2 = 0.0;
    DoubleVector V;

    double Te = 3.0;
    double Ti = 2.0;

    DoubleVector *pk;
    DoubleVector *pz;
    DoubleVector *pe;

    DoubleMatrix *pu;
    DoubleMatrix *pp;

    DoubleVector ks;
    DoubleVector zs;
    DoubleVector es;

};

#endif // PROBLEM1KZX_H
