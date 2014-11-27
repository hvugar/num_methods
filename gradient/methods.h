#ifndef __METHODS_H
#define __METHODS_H

#include "gradient.h"

#ifdef __cplusplus
extern "C" {
#endif

//Метод наискорейшего спуска
//Fast proximal gradient method
void fast_proximal_gradient_method(RnFunction f, R1Function g, double* x, int N, double dx, double epsilon);

#ifdef __cplusplus
}
#endif


#endif // __METHODS_H