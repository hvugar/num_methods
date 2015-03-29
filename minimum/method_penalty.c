#include "methods.h"
#include "print.h"

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
void penalty_method(RnFunction f, double *x, int n, RnFunction* h, int m, RnFunction* g, int p, double r1, double r2, double epsilon)
{
	double G(RnFunction g, double *x, int n)
	{
		double y = g(x,n);
		double s = (y >= 0.0) ? 0.0 : y;
		return s*s;
	}
	
	double H(RnFunction h, double *x, int n)
	{
		double y = h(x, n);
		return y*y;
	}
	
	double R(double *x, int n)
	{
		int i;
		double sum_g = 0.0;
		double sum_h = 0.0;
		for (i=0; i<p; i++) sum_g = sum_g + G(g[i], x, n);
		for (i=0; i<m; i++) sum_h = sum_h + H(h[i], x, n);
		return sum_g*r2 + sum_h*r2;
	}
	
	double P(double *x, int n)
	{	
		return f(x,n) + R(x,n);
	}

	// Qoshma qradient usulu ucun parametrler
	double min_epsilon = 0.001;       //dovrun sona catma meyari
	double grad_step   = 0.005;       //gradient
	double line_step   = 0.1;         //parcani bolme
	double gold_step   = 0.0001;      //qizil qayda ucun
	
	double* x1 = (double*) malloc( sizeof(double) * n );
	do
	{
		memcpy( x1, x, sizeof(double) * n );
		printf("Minimization...\n");
		conjugate_gradient_method(P, x, n, line_step, gold_step, grad_step, min_epsilon, printer2, NULL);
		printf("Minimized...\n");

		//printf("r1 = %.10f\n", r1);
		printf("R  = %.1f\n", r2);
		printf("x1 = %.10f\n", x[0]);
		printf("x2 = %.10f\n", x[1]);
		printf("fx = %.10f\n", f(x,n));
		
		//r1 = r1 * 0.10;
		r2 = r2 * 10.0;
		
	} while ( distance(x1, x, n) > epsilon*0.001 );
	free(x1);
}

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
void penalty_method1(RnFunction f, double *x, int n, RnFunction* h, int m, RnFunction* g, int p, double r1, double r2, double epsilon)
{
	double G(RnFunction g, double *x, int n)
	{
		return - exp (g(x,n));
	}
	
	double H(RnFunction h, double *x, int n)
	{
		return h(x,n) * h(x,n);
	}
	
	double R(double *x, int n)
	{
		int i;
		double sum = 0.0;
		for (i=0; i<p; i++) sum = sum + r1 * G(g[i], x, n);
		for (i=0; i<m; i++) sum = sum + r2 * H(h[i], x, n);
		return sum;
	}
	
	double P(double *x, int n)
	{	
		return f(x,n) + R(x,n);
	}

	// Qoshma qradient usulu ucun parametrler
	double min_epsilon = 0.001;       //dovrun sona catma meyari
	double grad_step   = 0.005;       //gradient
	double line_step   = 0.1;         //parcani bolme
	double gold_step   = 0.0001;      //qizil qayda ucun
	
	double* x1 = (double*) malloc( sizeof(double) * n );
	do
	{
		memcpy( x1, x, sizeof(double) * n );
		printf("\nr1 = %.10f\nr2 = %.10f\n", r1, r2);
		printf("Minimization...\n");
		conjugate_gradient_method(P, x, n, line_step, gold_step, grad_step, min_epsilon, printer2, NULL);
		printf("Minimized...\n");
		printf("********************************************************\n");
		
		r1 = r1 * 0.10;
		r2 = r2 * 10.0;
	} while ( distance(x1, x, n) > epsilon );
	free(x1);
	printf("x1 = %.10f\nx2 = %.10f\nf  = %.10f\n", x[0], x[1], f(x,n));
}

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
void penalty_method2(RnFunction f, double *x, int n, RnFunction* h, int m, RnFunction* g, int p, double r1, double r2, double epsilon)
{
	double G(RnFunction g, double *x, int n)
	{
		return 1.0 / g(x,n);
	}
	
	double H(RnFunction h, double *x, int n)
	{
		return h(x,n) * h(x,n);
	}
	
	double R(double *x, int n)
	{
		int i;
		double sum = 0.0;
		for (i=0; i<p; i++) sum = sum + r1 * G(g[i], x, n);
		for (i=0; i<m; i++) sum = sum + r2 * H(h[i], x, n);
		return sum;
	}
	
	double P(double *x, int n)
	{	
		return f(x,n) + R(x,n);
	}

	/*
	// Qoshma qradient usulu ucun parametrler
	double min_epsilon = 0.001;       //dovrun sona catma meyari
	double grad_step   = 0.005;       //gradient
	double line_step   = 0.1;         //parcani bolme
	double gold_step   = 0.0001;      //qizil qayda ucun
	
	while ( r1 * R(x,n) > epsilon_p ) 
	{
		printf("\nr1 = %.10f\nr2 = %.10f\n", r1, r2);
		printf("Minimization...\n");
		conjugate_gradient_method(P, x, n, line_step, gold_step, grad_step, min_epsilon, printer2);
		printf("Minimized...\nx1 = %.10f\nx2 = %.10f\n", x[0], x[1]);
		printf("********************************************************\n");
		
		r1 = r1 * 0.10;
		r2 = r2 * 10.0;
	}
	*/
}
