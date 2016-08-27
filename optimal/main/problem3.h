#ifndef PROBLEM3_H
#define PROBLEM3_H

#include "newtonheatprocess.h"
#include <cmethods.h>
#include <function.h>
#include <printer.h>
#include <gradient_cjt.h>
#include <vector>

class Problem3 : public RnFunction, public IGradient, public IPrinter, public ConjugateGradient
{
public:
    Problem3();
    virtual ~Problem3();

    virtual double fx(const DoubleVector &x);
    virtual void gradient(const DoubleVector &v, DoubleVector &g);
    virtual void print(unsigned int iteration, const DoubleVector &k, const DoubleVector &g, double alpha, RnFunction *fn) const;

    virtual double v(unsigned int j) const;
    virtual double vs(unsigned int j) const;
    //virtual double vm(unsigned int j) const;
    //virtual double vr(unsigned int j) const;

    virtual double initial(unsigned int i) const;

    virtual double mu(unsigned int i) const;

    void calculateU(DoubleMatrix &m, double ht, double hx, unsigned int M, unsigned int N, double lambdaM, double lambdaL, double lambdaR, double a);
    void calculateP(DoubleMatrix &m, double ht, double hx, unsigned int M, unsigned int N, double lambdaM, double lambdaL, double lambdaR, double a);
    void calculateV(const DoubleVector &k);


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

    double lambdaM;
    double lambdaL;
    double lambdaR;

    unsigned int L;
    double Te;
    double Ti;

    double alpha1;
    double alpha2;

    DoubleVector V;
    const DoubleVector *pv;
    const DoubleMatrix *pu;
    const DoubleMatrix *pp;

    DoubleVector ks;
};

#endif // PROBLEM3_H
