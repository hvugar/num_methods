#ifndef PFUNCTIONAL_H
#define PFUNCTIONAL_H

#include "ifunctional.h"
#include <gradient_cjt.h>
#include <gradient_sd.h>
#include <gradient_cs.h>

class PFunctional
{
public:
    PFunctional();

    void calculate(DoubleVector &x, double R);
    void calculate(DoubleVector &x, const DoubleVector &r, const DoubleVector &e);

    double rfactor;

    IFunctional* func;
    GradientMethod *grad;
};

#endif // PFUNCTIONAL_H
