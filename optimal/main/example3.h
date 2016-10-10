#ifndef EXAMPLE3_H
#define EXAMPLE3_H

#include <gradient_cjt.h>
#include <function.h>
#include <gradient.h>
#include <vector2d.h>
#include <matrix2d.h>
#include <printer.h>
#include <projection.h>

class Example3 : public RnFunction, public IGradient, public ConjugateGradient, public IPrinter, public IProjection
{
public:
    Example3();
    virtual ~Example3() {}

    virtual double fx(const DoubleVector &x);
    virtual void gradient(const DoubleVector &x, DoubleVector &g);
    virtual void print(unsigned int iteration, const DoubleVector &x, const DoubleVector &gradient, double alpha, RnFunction *fn) const;
    virtual void project(DoubleVector &x, int index);

    double mu(unsigned int i UNUSED_PARAM) const { return 1.0; }
    // qovma 1
    void calculateU(DoubleMatrix &u, const DoubleVector &x);
    void calculateP(DoubleMatrix& p, const DoubleMatrix &u, const DoubleVector &x);

    double initial(unsigned int i) const;

    //GaussianElimination
    void calculateU1(DoubleMatrix &u, const DoubleVector &x);
    // qovma E
    void calculateU2(DoubleMatrix &u, const DoubleVector &x);
    // teze qovma
    void calculateU3(DoubleMatrix &u, const DoubleVector &x);
    void calculate1(unsigned int N, unsigned int K, const DoubleMatrix &a, DoubleMatrix &beta, DoubleVector &qamma, DoubleVector &x);


private:
    unsigned int N = 1000;
    unsigned int M = 1000;
    double hx = 0.001;
    double ht = 0.001;
    unsigned int L = 2;
    double a = 1.0;

    double Ti = 2.0;
    double Te = 3.0;
    double alpha = 1.0;
    double lambda0 = 1.0;
    double lambdal = 1.0;

    double alpha0 = 1.0;
    double alpha1 = 1.0;
    double alpha2 = 1.0;
    double alpha3 = 1.0;

    const DoubleVector *px;
    const DoubleVector *sx;
    DoubleVector V;

//    DoubleVector k;
//    DoubleVector z;
//    DoubleVector e;

public:
    static void Main(int argc, char* argv[]);
};

#endif // EXAMPLE3_H
