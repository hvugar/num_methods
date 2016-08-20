#ifndef PROBLEM2_H
#define PROBLEM2_H

#include "newtonheatprocess.h"
#include <cmethods.h>
#include <function.h>
#include <printer.h>
#include <gradient_cjt.h>
#include <vector>

class Problem2 : public RnFunction, public IGradient
{
public:
    Problem2();
    virtual ~Problem2();

    virtual double fx(const DoubleVector &x);
    virtual void gradient(const DoubleVector &x, DoubleVector &g);

    virtual double vm(unsigned int j) const;
    virtual double vl(unsigned int j) const;
    virtual double vr(unsigned int j) const;

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
    //DoubleVector k;
    DoubleVector xi;
    DoubleVector z;
    std::vector<unsigned int> Xi;
    double Te;
    double Ti;


    double alpha1;
    double alpha2;

    DoubleVector V;
    const DoubleVector *pk;
    const DoubleMatrix *pu;
    const DoubleMatrix *pp;
};

#endif // PROBLEM2_H
