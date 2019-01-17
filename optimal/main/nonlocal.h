#ifndef NONLOCAL_H
#define NONLOCAL_H

#include <vector2d.h>
#include <limits>
#include <cmath>
#include <cfloat>

//class NonLocal
//{
//public:
//    NonLocal();
//    virtual ~NonLocal();

//    void solve(DoubleVector &c, double d, double h);

//    virtual double a(double t) const;
//    virtual double b(double t) const;

//    bool equalZero(double) const;
//};

//NonLocal::~NonLocal() {}

//bool NonLocal::equalZero(double x) const
//{
//    return fabs(x) <= DBL_EPSILON;
//}

//void NonLocal::solve(DoubleVector &c, double d, double h)
//{
//    unsigned int N = c.length();
//    if (N < 2) throw std::exception();
//    if (equalZero(c[0])) throw std::exception();
//    if (equalZero(c[N])) throw std::exception();

//    unsigned int M = 0;
//    for (unsigned int i=2; i<=N; i++) if (!equalZero(c[i])) { M++; }

//    double *p = new double[N];
//    double *q = new double[N];
//    double *k = new double[N];
//    double **C = new double*[M];

//    p[0] = c[0]; q[0] = c[1]; k[0] = d;

//    for (unsigned int i=1; i<=N; i++)
//    {
//        double t = (i-1)*h;
//        double m = 3.0+2.0*h*a(t);
//        double a1 = +4.0/m;
//        double a2 = -1.0/m;
//        double a3 = -b(t)/m;

//        p[i] = p[i-1]*a1 + q[i-1];
//        q[i] = p[i-1]*a2;
//        k[i] = k[i-1] - p[i-1]*a2;


//    }
//}

#endif // NONLOCAL_H
