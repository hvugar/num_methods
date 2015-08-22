#ifndef GRADIENT_H
#define GRADIENT_H

#include "function.h"
#include "methods.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief      Градиент функции
 * @param f    Целевая функция
 * @param dx   Длина шагов для нахождение градиента
 * @param x    Независимые переменные
 * @param n    Число переменных
 * @param gr   Вектор градиента функции
 */
void Gradient(RnFunction f, double dx, double *x, int n, double *gr);

/**
 * @brief Метод наискорейшего градиентного спуска.
 * @param f Целевая функция
 * @param g
 * @param x Независимые переменные
 * @param n Число переменных
 * @param line_step Длина шагов метода прямого поиска
 * @param gold_eps  Число эпсилон для останова метода золотого сечение
 * @param grad_step Длина шагов для нахождение градиента
 * @param epsilon   Число эпсилон для останова метода наискорейшего спуска
 */
void SteepestDescentMethod(RnFunction f, GFunction g, double *x, int n, double line_step, double gold_eps, double grad_eps, double epsilon1, double espsilon2);


/**
 * @brief Метод Флетчера-Ривса (Метод сопряженных градиентов).
 * @param f Целевая функция
 * @param g
 * @param x Независимые переменные
 * @param n Число переменных
 * @param line_step Длина шагов метода прямого поиска
 * @param gold_eps  Число эпсилон для останова метода золотого сечение
 * @param grad_step Длина шагов для нахождение градиента
 * @param epsilon   Число эпсилон для останова метода наискорейшего спуска
 */
void ConjugateGradientMethod(RnFunction f, GFunction g, double *x, int n, double line_step, double gold_eps, double grad_eps, double epsilon1, double espsilon2);

#ifdef __cplusplus
}
#endif

#endif // GRADIENT_H
