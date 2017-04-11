#ifndef ART_PROBLEM1L2_H
#define ART_PROBLEM1L2_H

#include <function.h>
#include <gradient.h>
#include <vector2d.h>
#include <matrix2d.h>
#include <printer.h>
#include <projection.h>
#include <gradient_cjt.h>
#include <gradient_sd.h>
#include <cmethods.h>
#include <gradient/igradient.h>

#include "../iproblem1.h"

class ArtProblem1L2 : protected IProblem1/*protected RnFunction, protected IGradient,*/, public IPrinter, public IProjection
{
public:
    ArtProblem1L2();
    virtual ~ArtProblem1L2() {}

    void initialize();
    void startOptimize();

    void optimize(DoubleVector &x0) const;

//    virtual double fx(const DoubleVector &y) const;
//    virtual void gradient(const DoubleVector &y, DoubleVector &g);

    virtual void print(unsigned int i, const DoubleVector &y, const DoubleVector &g, double fx, GradientMethod::MethodResult result) const;
    virtual void project(DoubleVector &x, int index);

//    void calculateU(DoubleMatrix &u, const DoubleVector &k, const DoubleVector &z, const DoubleVector &e) const;
//    void calculateP(DoubleMatrix &p, const DoubleMatrix &u, const DoubleVector &k, const DoubleVector &z, const DoubleVector &e);

//    double initial(unsigned int n) const;
//    double mu(unsigned int i UNUSED_PARAM) const { return 1.0; }

//    void getParameters(DoubleVector &k, DoubleVector &z, DoubleVector &e, const DoubleVector &x) const;
//    void qovmaFirstColM(double *a, double *b, double *c, double *d, double *x, unsigned int n, double *E) const;
//    void qovmaFirstRowM(double *a, double *b, double *c, double *d, double *x, unsigned int n, double *E) const;

//    double vf(unsigned int m, const DoubleVector &k, const DoubleVector &z, const DoubleVector &e, const DoubleMatrix &u) const;

//    double vd0(unsigned int m, const DoubleVector &k, const DoubleVector &z, const DoubleVector &e, const DoubleMatrix &u) const;
//    double gf(unsigned int m, const DoubleVector &k, const DoubleVector &z, const DoubleVector &e, const DoubleMatrix &u) const;
//    double sgn_min(unsigned int m, const DoubleVector &k, const DoubleVector &z, const DoubleVector &e, const DoubleMatrix &u) const;
//    double delta(unsigned int n, const DoubleVector &e, unsigned int i) const;
//    double sign(double x) const;

    void table1Generate();
    void table2Generate();
    void table3Generate();

    void image1Generate();
    void image2Generate();
    void image3Generate();

    void imager2L();
    void imager3L();

private:
//    unsigned int L;
//    unsigned int N;
//    unsigned int M;
//    double hx;
//    double ht;

    double hk;
    double hz;
    double he;

//    double fi;
//    double tt;

//    DoubleVector vfi;
//    DoubleVector vtt;

//    double a;
//    double lambda0;
//    double lambda1;
//    double lambda2;

//    double alpha0;
//    double alpha1;
//    double alpha2;
//    double alpha3;
//    double R;

//    double vmin;
//    double vmax;
//    double d0;
//    double d1;

//    DoubleVector V;
//    const DoubleVector *py;

//    DoubleVector K;
//    DoubleVector Z;
//    DoubleVector E;

//    DoubleVector k0;
//    DoubleVector z0;
//    DoubleVector e0;

//    bool optimizeK;
//    bool optimizeZ;
//    bool optimizeE;

    bool withError = false;
    bool hello = false;
    double persent = 0.01;

//    unsigned int DD;

//    double u_xi(unsigned int m, double xi, const DoubleMatrix &u) const;
//    double u_xi_d(unsigned int m, double xi, const DoubleMatrix &u) const;

public:
    static void Main(int argc, char* argv[]);
};

#endif // ART_PROBLEM1L2_H
