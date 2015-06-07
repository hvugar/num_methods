#ifndef CFUNCTION_H
#define CFUNCTION_H

#include "function.h"

void printX(char* s, std::vector<double> x);

struct CFunction : public RnFunction
{
    double t0;
    double t1;
    double h;
    int N;
    std::vector<double> t;
    std::vector<double> x1;
    std::vector<double> x2;
    std::vector<double> psi1;
    std::vector<double> psi2;
    double fx0(double t, std::vector<double> x, double u);
    double fx1(double t, std::vector<double> x, double u);
    double fx2(double t, std::vector<double> x, double u);
    double fp1(double t, std::vector<double> x, std::vector<double> psi, double u);
    double fp2(double t, std::vector<double> x, std::vector<double> psi, double u);
    double gradJ(double t, std::vector<double> x, std::vector<double> psi, double u);
    void gradientJ(double grad_step, std::vector<double>& g, const std::vector<double>& u);
    double H(double t, std::vector<double> x, double u, int r, std::vector<double> psi);
    virtual double fx(const std::vector<double>& u);

    CFunction(double t0, double t1, double h);
    virtual ~CFunction();

    void test(const std::vector<double>& u);
};

#endif // CFUNCTION_H
