#ifndef C_GRADIENT_H
#define C_GRADIENT_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <global.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef double (*c_function)(double *x, unsigned int n);
typedef void (*c_gradient)(double *x, double *g, unsigned int n);

struct c_gradient_t {
    unsigned int count;
    double grad_norm_esp;
    double distance_eps;
    double func_diff_eps;
    double min_step;
    double min_epsilon;
    c_function func;
    c_gradient grad;
    double *x;
    unsigned int n;
    int normalize;
    int showMsg;
};

double c_vector_L2Norm(double *x, unsigned int n);
void c_vector_L2Normalize(double *x, unsigned int n);

void conjucate_gradient(struct c_gradient_t *gr);

void steepest_descent_gradient(struct c_gradient_t *gr);

//void stranghLineSearch(double x, double step, double &a, double &b, R1Function *fn);
//void goldenSectionSearch(double &a, double &b, double &x, R1Function *f, double epsilon);


#ifdef __cplusplus
}
#endif


#endif // C_GRADIENT_H
