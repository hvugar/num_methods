#ifndef CXX_METHODS_H
#define CXX_METHODS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef double (*R1Function)(double);
typedef double (*R2Function)(double x, double y);
typedef double (*RnFunction)(double *x, int n);


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

#endif // CXX_METHODS_H

