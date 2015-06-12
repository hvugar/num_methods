#include "methods.h"

double minimize1(RnFunction f, double *x, double *grad, int n, double line_eps, double gold_eps)
{
	double *x1 = (double*) malloc(sizeof(double) * n);
    double argmin(double alpha)
    {
        int j;
        for (j=0; j<n; j++) x1[j] = x[j] - alpha * grad[j];
        double result = f(x1, n);
        return result;
    }
	double alpha0 = 0.0;
    double a = 0.0;
	double b = 0.0;
    straight_line_search_metod(argmin, alpha0, line_eps, &a, &b);
    double alpha = golden_section_search_min(argmin, a, b, gold_eps);
	if ( argmin(alpha) > argmin(alpha0) ) alpha = alpha0;
	free(x1);
    return alpha;
}

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
void steepest_descent_gradient_method(RnFunction f, double *x, int n, double step, double gold_eps, double grad_eps, double epsilon, Printer printer)
{
    int i = 0;
	int k = 0;
    double grad_norm = 0.0;
	double dstnc = 0.0;
	double alpha = 0.0;

    double* grad = (double*) malloc(sizeof(double) * n);
    
	do
    {   
		/* calculating function gradient at current point */
        gradient(f, x, n, grad_eps, grad);
		
		/* if gradinet norm at current point is less than epsilon then break. no minimize */
		grad_norm = vertor_norm(grad, n);
		//if (grad_norm < epsilon) break;
		
		k++;
		
		/* calculating unit vectors */
        for (i=0; i<n; i++) grad[i] = grad[i] / grad_norm;
		
		/* R1 minimization in direct of antigradient */
		alpha = minimize1(f, x, grad, n, step, gold_eps);

        if (printer != NULL) printer(f, x, n, k-1, grad, grad_norm, alpha);

		dstnc = 0.0;
        for (i=0; i<n; i++)
        {
            x[i] = x[i] - alpha * grad[i];
			dstnc = dstnc + (alpha * grad[i])*(alpha * grad[i]);
        }
		dstnc = sqrt(dstnc);
    }
    while ( dstnc > epsilon );
	if (printer != NULL) printer(f, x, n, k-1, grad, grad_norm, alpha);

    free(grad);
}
