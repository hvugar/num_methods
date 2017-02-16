#ifndef HYPERBOLICCONTROL2DM_H
#define HYPERBOLICCONTROL2DM_H

#include <function.h>
#include <hyperbolicequation.h>
#include <printer.h>
#include <projection.h>
#include <gradient_cjt.h>


class MINIMUMSHARED_EXPORT HyperbolicControl2DM : public R1Function, public RnFunction,
        public IGradient, public IHyperbolicEquation2D, public IBackwardHyperbolicEquation2D,
        public IPrinter, public IProjection
{
public:
    HyperbolicControl2DM() {}
    virtual ~HyperbolicControl2DM() {}

    virtual double fx(double x) const;
    virtual double fx(const DoubleVector &x) const;
    virtual void gradient(const DoubleVector &v, DoubleVector &e);

    virtual double initial1(unsigned int i, unsigned int j) const;
    virtual double initial2(unsigned int i, unsigned int j) const;
    virtual double boundary(unsigned int i, unsigned int j, unsigned int k) const;
    virtual double f(unsigned int i, unsigned int j, unsigned int k) const;

    virtual double binitial1(unsigned int i, unsigned int j) const;
    virtual double binitial2(unsigned int i, unsigned int j) const;
    virtual double bboundary(unsigned int i, unsigned int j, unsigned int k) const;
    virtual double bf(unsigned int i, unsigned int j, unsigned int k) const;

    virtual void print(unsigned int i, const DoubleVector &v, const DoubleVector &g, double fx) const;
    virtual void project(DoubleVector &x, int index);
    double fxt(unsigned int i, unsigned int j, unsigned int k) const;
    double norm(const DoubleVector& v) const;

    void calculateGX(const DoubleVector& x, const DoubleMatrix& psi, DoubleVector& g, unsigned int k);
    void psiDerivative(double &psiX1, double &psiX2, double x1, double x2, const DoubleMatrix &psi);

    static void Main(int argc, char *argv[]);

private:
    double h1;
    double h2;
    double ht;
    unsigned int L;
    unsigned int N1;
    unsigned int N2;
    unsigned int M;

    double alpha0;
    double alpha1;
    double alpha2;
    double alpha3;
    double qamma;
    double a;

    double U0;
    double U1;

    //DoubleVector d;
    DoubleVector e;
    const DoubleVector *px;
    const DoubleCube *pu;
    double h;

    double vd;
    double vu;
};

#endif // HYPERBOLICCONTROL2DM_H
