#ifndef _DIFF_METHOD_H
#define _DIFF_METHOD_H

typedef double (*FxtFunction)(double x, double t);
typedef double (*FxytFunction)(double x, double y, double t);
typedef double (*FxyztFunction)(double x, double y, double z, double t);

//Initial condition function type
typedef double (*fiFunction)(double x);

//Boundary condition function type
typedef double (*m1Function)(double t);
typedef double (*m2Function)(double t);

void implicit_difference_scheme(FxtFunction f, fiFunction fi, m1Function m1, m2Function m2, double a, double delta_x, double delta_t, double X, double T);

#endif //_DIFF_METHOD_H