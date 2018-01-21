#ifndef PFUNCTIONAL_H
#define PFUNCTIONAL_H

#include "jfunctional.h"
#include <gradient_cjt.h>
#include <gradient_sd.h>

class PFunctional
{
public:
    PFunctional();

    void calculate(DoubleVector &x, double R);
    void calculate(DoubleVector &x);

    double rfactor;
    JFunctional* jfunc;
};

#endif // PFUNCTIONAL_H
