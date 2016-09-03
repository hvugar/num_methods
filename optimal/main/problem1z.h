#ifndef PROBLEM1Z_H
#define PROBLEM1Z_H

#include <cmethods.h>
#include <function.h>
#include <printer.h>
#include <gradient_cjt.h>
#include <vector>
#include <rungekutta.h>

class Problem1Z : public RnFunction, public IGradient, public IPrinter, public ConjugateGradient
{
public:
    Problem1Z();
    virtual ~Problem1Z();

    virtual double fx(const DoubleVector &x);
    virtual void gradient(const DoubleVector &x, DoubleVector &g);
    virtual void print(unsigned int iteration, const DoubleVector &k, const DoubleVector &g, double alpha, RnFunction *fn) const;

    virtual double vm(unsigned int j) const;
    virtual double vl(unsigned int j) const;
    virtual double vr(unsigned int j) const;

    virtual double initial(unsigned int i) const;

    virtual double mu(unsigned int i) const;

    void calculateU(DoubleMatrix &m, double ht, double hx, unsigned int M, unsigned int N, double alpha, double lambda0, double lambdal, double a);
    void calculateU2(DoubleMatrix &m, double ht, double hx, unsigned int M, unsigned int N, double alpha, double lambda0, double lambdal, double a);

    void calculateP(DoubleMatrix &m, double ht, double hx, unsigned int M, unsigned int N, double alpha, double lambda0, double lambdal, double a);
    void calculateP2(DoubleMatrix &m, double ht, double hx, unsigned int M, unsigned int N, double alpha, double lambda0, double lambdal, double a);

    void calculateV(const DoubleVector &k);

    double vs(double j) const;

    static void Main(int argc, char* argv[]);

public:
    double t0;
    double t1;
    double x0;
    double x1;
    double hx;
    double ht;
    unsigned int M;
    unsigned int N;
    double a;

    double alpha;
    double lambda0;
    double lambdal;
    double alpha1;
    double alpha2;

    double Te;

    unsigned int L;
    DoubleVector k;
    DoubleVector xi;
    std::vector<unsigned int> Xi;

    DoubleVector V;
    const DoubleVector *pz;
    const DoubleMatrix *pu;
    const DoubleMatrix *pp;

    DoubleVector zs;
};

#endif // PROBLEM1Z_H
