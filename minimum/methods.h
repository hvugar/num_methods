#ifndef __METHODS_H
#define __METHODS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef double (*R1Function)(double);
typedef double (*R2Function)(double x, double y);
typedef double (*RnFunction)(double*, int);
typedef double (*RmFunction)(double, double*, int);
typedef void   (*Printer)(RnFunction f, double *x, int n, ...);
typedef void   (*GetInfo)(RnFunction f, double *x, int n, int iteration, double *grad, double *s, R1Function min, double alpha, double a, double b);

typedef struct {
	R1Function func;
	double x;
} FunctionR1;

typedef struct {
	RnFunction func;
	double *x;
	int n;
} FunctionRn;


/**
 * @brief Градиент функции
 * @param f         Целевая функция
 * @param x         Независимые переменные
 * @param n         Число переменных
 * @param dx        Длина шагов для нахождение градиента
 * @param gradients Вектор градиента функции
 */
void gradient(RnFunction f, double *x, int n, double dx, double *gradients);


/**
 * @brief Градиент функции по формуле f(x+dx)-f(x) / dx
 * @param f         Целевая функция
 * @param x         Независимые переменные
 * @param n         Число переменных
 * @param dx        Длина шагов для нахождение градиента
 * @param gradients Вектор градиента функции
 */
void gradient1(RnFunction f, double *x, int n, double dx, double *gradients);

/**
 * @brief Градиент функции по формуле f(x+dx)-f(x-dx) / 2dx
 * @param f         Целевая функция
 * @param x         Независимые переменные
 * @param n         Число переменных
 * @param dx        Длина шагов для нахождение градиента
 * @param gradients Вектор градиента функции
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
 * @brief Норма вектора
 * @param vctr Вектор
 * @param n    Число елементов вектора
 * @return
 */
double vertor_norm(double *vctr, int n);

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
 * @brief Методы прямого поиска
 * @param f       Целевая функция
 * @param a       Начальная точка отрезка
 * @param b       Конечнная точка отрезка
 * @param epsilon Число эпсилон для останова метода
 * @return
 */
double straight_line_search_metod(R1Function fx, double x0, double dx, double *a, double *b);

/**
 * @brief Метод золотого сечения
 * @param fx      Целевая функция
 * @param a       Начальная точка отрезка
 * @param b       Конечнная точка отрезка
 * @param epsilon Число эпсилон для останова метода золотого сечение
 * @return
 */
double golden_section_search_min(R1Function fx, double a, double b, double epsilon);

/**
 * @brief Метод наискорейшего спуска
 * @param f         Целевая функция
 * @param x         Независимые переменные
 * @param n         Число переменных
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
void conjugate_gradient_method(RnFunction f, double *x, int n, double line_step, double gold_step, double grad_step, double epsilon, Printer printer, GetInfo info);

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
 * @brief         Метод штрафных функций
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
 * @brief         Метод барьерных функций [обратная штрафная функция]. Смешанной вспомогательная функция
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
 * @brief         Метод барьерных функций [логарифмическая штрафная функция]. Смешанной вспомогательная функция
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

double search_method_dck(R1Function f, double x0, double dx, double *a, double *b);

double search_method_pauella(R1Function f, double x0, double dx, double epsilon, double *a, double *b);

void search_interval_svenn(R1Function f, double x0, double dx, double *a, double *b);

void halph_interval_method(R1Function f, double epsilon, double *a, double *b);

double newton_raphson(R1Function f, double x0, double epsilon);

double runga_kutta1(R2Function y, double y0, double x0, double x, double h);

void runga_kutta2(RmFunction *f, double *y0, double x0, double *y1, double x1, int n, double h);

#ifdef __cplusplus
}
#endif


#endif // __METHODS_H
