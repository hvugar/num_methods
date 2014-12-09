#ifndef __METHODS_H
#define __METHODS_H

#include "gradient.h"

#ifdef __cplusplus
extern "C" {
#endif

//Метод наискорейшего спуска
void fast_proximal_gradient_method(RnFunction f, R1Function g, double *x, int n, double dx, double epsilon);

// Метод сопряженных градиентов Флетчера — Ривса 
void conjugate_gradient_method(RnFunction f, R1Function g, double *x, int n, double dx, double epsilon);

//Метод сопряженных градиентов
void conjugate_gradient_method1(RnFunction f, R1Function g, double *x, int n, double dx, double epsilon);


#ifdef __cplusplus
}
#endif


#endif // __METHODS_H