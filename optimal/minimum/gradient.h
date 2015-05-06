#ifndef GRADIENT_H
#define GRADIENT_H

#include "global.h"
#include "methods.h"

class MINIMUMSHARED_EXPORT Gradient
{
public:
    Gradient();

    virtual double fx(double *x, int n) = 0;

    void fast_proximal_gradient_method(RnFunction *f, double *x, int n, double line_step, double gold_step, double grad_step, double epsilon);
    void conjugate_gradient_method(RnFunction *f, double *x, int n, double line_step, double gold_step, double grad_step, double epsilon);

    void fast_proximal_gradient_method();
    double argmin(double);
    double R1Minimize();
protected:
    virtual void show();

private:
    int count;
    std::vector<double> x0;
    std::vector<double> x1;
    double *x;
    double *grads;
    double epsilon;
    int n;
    double alpha;
};

#endif // GRADIENT_H
