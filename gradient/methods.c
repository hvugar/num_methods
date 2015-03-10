#include "methods.h"

/**
 * @brief vertor_norm
 * @param x
 * @param n
 * @return
 */
double vertor_norm(double *vctr, int n)
{
    double sum = 0.0;
	int i;
    for (i=0; i<n; i++)
    {
        sum += vctr[i] * vctr[i];
    }
    return sqrt(sum);
}

/**
 * @brief distance
 * @param x1
 * @param x2
 * @param n
 * @return
 */
double distance(double *x1, double *x2, int n)
{
	double dist = 0.0;
	int i;
    for (i=0; i<n; i++)
    {
		double dx = x1[i] - x2[i];
        dist += dx*dx;
    }
    return sqrt(dist);
}

/**
 * @brief Градиент функции
 * @param f         Целевая функция
 * @param x         Независимые переменные
 * @param n         Число переменных
 * @param dx        Длина шагов для нахождение градиента
 * @param gradients Вектор градиента функции
 */
void gradient(RnFunction f, double *x, int n, double dx, double *gradients)
{
	gradient2(f, x, n, dx, gradients);
}

/**
 * @brief Градиент функции по формуле f(x+dx)-f(x) / dx
 * @param f         Целевая функция
 * @param x         Независимые переменные
 * @param n         Число переменных
 * @param dx        Длина шагов для нахождение градиента
 * @param gradients Вектор градиента функции
 */
void gradient1(RnFunction f, double *x, int n, double dx, double *gradients)
{
    double f0 = f(x, n);
	int i;
    for (i=0; i<n; i++)
    {
        x[i] = x[i] + dx;
		double f1 = f(x, n);
        gradients[i] = (f1 - f0) / dx;
        x[i] = x[i] - dx;
    }
}

/**
 * @brief Градиент функции по формуле f(x+dx)-f(x-dx) / 2dx
 * @param f         Целевая функция
 * @param x         Независимые переменные
 * @param n         Число переменных
 * @param dx        Длина шагов для нахождение градиента
 * @param gradients Вектор градиента функции
 */
void gradient2(RnFunction f, double *x, int n, double dx, double *gradients)
{
	int i = 0;
	for (i=0; i<n; i++)
	{
		x[i] = x[i] - dx;
		double f1 = f(x, n);
		x[i] = x[i] + 2*dx;
		double f2 = f(x, n);
		x[i] = x[i] - dx;
		gradients[i] = (f2 - f1) / (2 * dx);
	}
}