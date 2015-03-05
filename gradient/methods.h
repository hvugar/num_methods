#ifndef __METHODS_H
#define __METHODS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef double (*R1Function)(double);
typedef double (*RnFunction)(double*, int);
typedef void (*Printer)(RnFunction f, double *x, int n, ...);

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief
 * @param
 */
void gradient(RnFunction f, double *x, int n, double dx, double *gradients);

/**
 * @brief
 * @param
 */
double vertor_norm(double *x, int n);

/**
 * @brief
 * @param
 */
double grad_module(double *grads, int n);

/**
 * @brief
 * @param
 */
double distance(double *x1, double *x2, int n);

/**
 * @brief Метод золотого сечения
 * @param f
 * @param a
 * @param b
 * @param epsilon
 * @return
 */
double straight_line_search_metod(R1Function fx, double x0, double dx, double *a, double *b);

/**
 * @brief
 * @param f
 * @param x0
 * @param dx
 * @param a
 * @param b
 * @return
 */
double golden_section_search_min(R1Function fx, double a, double b, double epsilon);

/**
 * @brief
 * @param f
 * @param x0
 * @param dx
 * @param a
 * @param b
 * @return
 */
double search_method_dck(R1Function f, double x0, double dx, double *a, double *b);

/**
 * @brief
 * @param f
 * @param x0
 * @param dx
 * @param a
 * @param b
 * @return
 */
double search_method_pauella(R1Function f, double x0, double dx, double epsilon, double *a, double *b);

/**
 * @brief Этап установления границ интервала. Метод Свенна
 * @param f
 * @param x0
 * @param dx
 * @param a
 * @param b
 * @return
 */
void search_interval_svenn(R1Function f, double x0, double dx, double *a, double *b);

/**
 * @brief Метод деления интервала пополам
 * @param f
 * @param epsilon
 * @param a
 * @param b
 * @return
 */
void halph_interval_method(R1Function f, double epsilon, double *a, double *b);

/**
 * @brief Метод Ньютона - Рафсона
 * @param f
 * @param x0
 * @param epsilon
 * @return
 */
double newton_raphson(R1Function f, double x0, double epsilon);

/**
 * @brief
 * @param
 */
double derivative_1(R1Function f, double x, double h);

/**
 * @brief
 * @param
 */
double derivative_2(R1Function f, double x, double h);

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
void conjugate_gradient_method(RnFunction f, double *x, int n, double line_step, double gold_step, double grad_step, double epsilon, Printer printer);

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
void conjugate_gradient_method1(RnFunction f, double *x, int n, double line_step, double gold_step, double grad_step, double epsilon, Printer printer);

/**
 * @brief
 * @param
 * @param
 * @param
 * @param
 * @return
 */
void penalty_method(RnFunction f, double *x, int n, RnFunction* h, int m, RnFunction* g, int p, double r1, double r2, double epsilon);

#ifdef __cplusplus
}
#endif


#endif // __METHODS_H
