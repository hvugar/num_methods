#ifndef FUNCTION_H
#define FUNCTION_H

#ifdef __cplusplus
extern "C" {
#endif

#define M_E1 2.7182818284590452353602874713527

typedef double (*R1Function)(double);
typedef double (*R2Function)(double x, double y);
typedef double (*R3Function)(double x, double y, double z);
typedef double (*RnFunction)(double*, int);

typedef void (*GFunction)(RnFunction f, double grad_step, double* x, int n, double *g);

typedef double (*RmFunction)(double x, double *y, int n);
typedef void   (*Printer)(RnFunction f, double *x, int n, ...);

#ifdef __cplusplus
}
#endif

#endif // FUNCTION_H
