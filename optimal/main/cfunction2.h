#ifndef CFUNCTION2_H
#define CFUNCTION2_H

#include <math.h>
#include <function.h>
#include "utils.h"

class CFunction2 : public RnFunction
{
public:
    CFunction2(double t0, double t1, double h);
    virtual ~CFunction2();

    virtual double fx(const std::vector<double>& u);
    void gradientJ(double grad_step, std::vector<double>& g, const std::vector<double>& u);

    double t0;
    double t1;
    double h;
    int n;
    std::vector<double> t;
    std::vector<double> x1;
    std::vector<double> x2;
    std::vector<double> psi1;
    std::vector<double> psi2;

private:
    double fx0(double t, double x1, double x2, double u);
    double fx1(double t, double x1, double x2, double u);
    double fx2(double t, double x1, double x2, double u);
    double fp1(double t, double x1, double x2, double p1, double p2, double u);
    double fp2(double t, double x1, double x2, double p1, double p2, double u);
    double H(double t, double x1, double x2, double u, double p1, double p2);
    void test(const std::vector<double>& u);
};

#endif // CFUNCTION2_H
