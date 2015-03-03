#ifndef __METHODS_H
#define __METHODS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "minimum.h"
#include "gradient.h"
#include "print.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Метод наискорейшего спуска
 * @param f Целевая функция
 * @param x Независимые переменные n - измерение
 * @param n
 * @param line_eps
 * @param gold_eps
 * @param grad_eps
 * @param epsilon
 */
void fast_proximal_gradient_method(RnFunction f, double *x, int n, double line_step, double gold_step, double grad_step, double epsilon);

/**
 * @brief Метод сопряженных градиентов Флетчера — Ривса
 * @param f Целевая функция
 * @param x Независимые переменные n - измерение
 * @param n
 * @param line_eps
 * @param gold_eps
 * @param grad_eps
 * @param epsilon
 */
void conjugate_gradient_method(RnFunction f, double *x, int n, double line_step, double gold_step, double grad_step, double epsilon);

/**
 * @brief Метод сопряженных градиентов Флетчера — Ривса
 * @param f Целевая функция
 * @param x Независимые переменные n - измерение
 * @param n
 * @param line_step
 * @param gold_step
 * @param grad_step
 * @param epsilon
 */
void conjugate_gradient_method1(RnFunction f, double *x, int n, double line_step, double gold_step, double grad_step, double epsilon);

#ifdef __cplusplus
}
#endif


#endif // __METHODS_H
