#ifndef __METHODS_H
#define __METHODS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "minimum.h"
#include "gradient.h"

#ifdef __cplusplus
extern "C" {
#endif

//Метод наискорейшего спуска
void fast_proximal_gradient_method(RnFunction f, double *x, int n, double dx, double epsilon);

// Метод сопряженных градиентов Флетчера — Ривса 
void conjugate_gradient_method(RnFunction f, double *x, int n, double dx, double epsilon);


#ifdef __cplusplus
}
#endif


#endif // __METHODS_H