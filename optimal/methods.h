#ifndef CXX_METHODS_H
#define CXX_METHODS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef double (*R1Function)(double);
typedef double (*R2Function)(double x, double y);
typedef double (*RnFunction)(double *x, int n);

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief    Градиент функции
 * @param f  Целевая функция
 * @param x  Независимые переменные
 * @param n  Число переменных
 * @param dx Длина шагов для нахождение градиента
 * @param gr Вектор градиента функции
 */
void gradient(RnFunction f, double *x, int n, double dx, double *gr);

/**
 * @brief    Градиент функции по формуле f(x+dx)-f(x) / dx
 * @param f  Целевая функция
 * @param x  Независимые переменные
 * @param n  Число переменных
 * @param dx Длина шагов для нахождение градиента
 * @param gr Вектор градиента функции
 */
void gradient1(RnFunction f, double *x, int n, double dx, double *gradients);

/**
 * @brief    Градиент функции по формуле f(x+dx)-f(x-dx) / 2dx
 * @param f  Целевая функция
 * @param x  Независимые переменные
 * @param n  Число переменных
 * @param dx Длина шагов для нахождение градиента
 * @param gr Вектор градиента функции
 */
void gradient2(RnFunction f, double *x, int n, double dx, double *gradients);

/**
 * @brief derivative1
 * @param f
 * @param x
 * @param h
 * @return
 */
double derivative1(R1Function f, double x, double h);

/**
 * @brief derivative2
 * @param f
 * @param x
 * @param h
 * @return
 */
double derivative2(R1Function f, double x, double h);

/**
 * @brief    Норма вектора
 * @param vc Вектор
 * @param n  Число елементов вектора
 * @return
 */
double vertor_norm(double *vctr, int n);

/**
 * @brief grad_module
 * @param gr
 * @param n
 * @return
 */
double grad_module(double *gr, int n);

/**
 * @brief distance
 * @param x1
 * @param x2
 * @param n
 * @return
 */
double distance(double *x1, double *x2, int n);

/**
 * @brief         Методы прямого поиска
 * @param f       Целевая функция
 * @param a       Начальная точка отрезка
 * @param b       Конечнная точка отрезка
 * @param epsilon Число эпсилон для останова метода
 * @return
 */
double straight_line_search_metod(R1Function fx, double x0, double dx, double &a, double &b);

/**
 * @brief         Метод золотого сечения
 * @param fx      Целевая функция
 * @param a       Начальная точка отрезка
 * @param b       Конечнная точка отрезка
 * @param epsilon Число эпсилон для останова метода золотого сечение
 * @return
 */
double golden_section_search_min(R1Function fx, double a, double b, double epsilon);

#ifdef __cplusplus
}
#endif

#endif // CXX_METHODS_H

