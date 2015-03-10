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
 * @brief Gradient of function
 * @param f
 * @param x
 * @param n
 * @param dx
 * @param gradients
 */
void gradient(RnFunction f, double *x, int n, double dx, double *gradients);


/**
 * @brief Gradient of function with derivative formula f(x+dx)-f(x) / dx
 * @param f
 * @param x
 * @param n
 * @param dx
 * @param gradients
 */
void gradient1(RnFunction f, double *x, int n, double dx, double *gradients);

/**
 * @brief Gradient of function with derivative formula f(x+dx)-f(x-dx) / 2dx
 * @param f
 * @param x
 * @param n
 * @param dx
 * @param gradients
 */
void gradient2(RnFunction f, double *x, int n, double dx, double *gradients);

/**
 * @brief vertor_norm
 * @param x
 * @param n
 * @return
 */
double vertor_norm(double *x, int n);

/**
 * @brief grad_module
 * @param grads
 * @param n
 * @return
 */
double grad_module(double *grads, int n);

/**
 * @brief distance
 * @param x1
 * @param x2
 * @param n
 * @return
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
 * @brief golden_section_search_min
 * @param fx
 * @param a
 * @param b
 * @param epsilon Число эпсилон для останова метода золотого сечение
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
 * @brief derivative_1
 * @param f
 * @param x
 * @param h
 * @return
 */
double derivative_1(R1Function f, double x, double h);

/**
 * @brief derivative_2
 * @param f
 * @param x
 * @param h
 * @return
 */
double derivative_2(R1Function f, double x, double h);

/**
 * @brief Метод наискорейшего спуска
 * @param f Целевая функция
 * @param x Независимые переменные
 * @param n Число переменных
 * @param line_step Длина шагов метода прямого поиска
 * @param gold_eps  Число эпсилон для останова метода золотого сечение
 * @param grad_step Длина шагов для нахождение градиента
 * @param epsilon   Число эпсилон для останова метода наискорейшего спуска
 */
void fast_proximal_gradient_method(RnFunction f, double *x, int n, double line_step, double gold_step, double grad_step, double epsilon, Printer printer);

/**
 * @brief Метод сопряженных градиентов Флетчера — Ривса
 * @param f         Целевая функция
 * @param x         Независимые переменные
 * @param n         Число переменных
 * @param line_step Длина шагов метода прямого поиска
 * @param gold_eps  Число эпсилон для останова метода золотого сечение
 * @param grad_step Длина шагов для нахождение градиента
 * @param epsilon   Число эпсилон для останова метода сопряженных градиентов
 */
void conjugate_gradient_method(RnFunction f, double *x, int n, double line_step, double gold_step, double grad_step, double epsilon, Printer printer);

/**
 * @brief Метод сопряженных градиентов Флетчера — Ривса
 * @param f         Целевая функция
 * @param x         Независимые переменные
 * @param n         Число переменных
 * @param line_step Длина шагов метода прямого поиска
 * @param gold_eps  Число эпсилон для останова метода золотого сечение
 * @param grad_step Длина шагов для нахождение градиента
 * @param epsilon   Число эпсилон для останова метода сопряженных градиентов
 */
void conjugate_gradient_method1(RnFunction f, double *x, int n, double line_step, double gold_step, double grad_step, double epsilon, Printer printer);

/**
 * @brief Методы штрафных функций
 * @param f       Целевая функция
 * @param x       Независимые переменные
 * @param n       Число переменных
 * @param h       Функции ограничений в виде равенств
 * @param m       Число ограничений в виде равенств
 * @param g       Функции ограничений в виде неравенств
 * @param p       Число ограничений в виде неравенств
 * @param r1      Штрафной коэффициент для ограничений в виде равенств
 * @param r2      Штрафной коэффициент для ограничений в виде неравенств
 * @param epsilon Число эпсилон для останова метода
 */
void penalty_method(RnFunction f, double *x, int n, RnFunction* h, int m, RnFunction* g, int p, double r1, double r2, double epsilon);

/**
 * @brief Методы штрафных функций
 * @param f       Целевая функция
 * @param x       Независимые переменные
 * @param n       Число переменных
 * @param h       Функции ограничений в виде равенств
 * @param m       Число ограничений в виде равенств
 * @param g       Функции ограничений в виде неравенств
 * @param p       Число ограничений в виде неравенств
 * @param r1      Штрафной коэффициент для ограничений в виде равенств
 * @param r2      Штрафной коэффициент для ограничений в виде неравенств
 * @param epsilon Число эпсилон для останова метода
 */
void penalty_method1(RnFunction f, double *x, int n, RnFunction* h, int m, RnFunction* g, int p, double r1, double r2, double epsilon);

/**
 * @brief Методы штрафных функций
 * @param f       Целевая функция
 * @param x       Независимые переменные
 * @param n       Число переменных
 * @param h       Функции ограничений в виде равенств
 * @param m       Число ограничений в виде равенств
 * @param g       Функции ограничений в виде неравенств
 * @param p       Число ограничений в виде неравенств
 * @param r1      Штрафной коэффициент для ограничений в виде равенств
 * @param r2      Штрафной коэффициент для ограничений в виде неравенств
 * @param epsilon Число эпсилон для останова метода
 */
void penalty_method2(RnFunction f, double *x, int n, RnFunction* h, int m, RnFunction* g, int p, double r1, double r2, double epsilon);

#ifdef __cplusplus
}
#endif


#endif // __METHODS_H
