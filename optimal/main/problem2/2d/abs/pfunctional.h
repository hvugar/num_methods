#ifndef PFUNCTIONAL_H
#define PFUNCTIONAL_H

#include "jfunctional.h"
#include <gradient_cjt.h>

class PFunctional
{
public:
    PFunctional();

    void calculate(DoubleVector &x, double R);

    double rfactor;
    JFunctional* jfunc;
};

#endif // PFUNCTIONAL_H
