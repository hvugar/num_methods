#ifndef GRADIENT_H
#define GRADIENT_H

#include "methods.h"

class Gradient
{
public:
    Gradient();

    void fast_proximal_gradient_method(RnFunction *f, double *x, int n, double line_step, double gold_step, double grad_step, double epsilon);
    void conjugate_gradient_method(RnFunction *f, double *x, int n, double line_step, double gold_step, double grad_step, double epsilon);
};

#endif // GRADIENT_H
