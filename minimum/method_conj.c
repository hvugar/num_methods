#include "methods.h"
#include "print.h"

double minimize2(RnFunction f, double *x, double *s, int n, double step, double gold_step)
{
	double *x1 = (double*) malloc(sizeof(double) * n);
    double argmin(double alpha)
    {
        int j;
        for (j=0; j<n; j++) x1[j] = x[j] + alpha * s[j];
        double result = f(x1, n);
        return result;
    }
    
	double alpha0 = 0.0;
    double a = 0.0;
	double b = 0.0;
    straight_line_search_metod(argmin, alpha0, step, &a, &b);
    double alpha = golden_section_search_min(argmin, a, b, gold_step);
	if (argmin(alpha)>argmin(alpha0)) alpha = alpha0;
	free(x1);
    return alpha;
}

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
void conjugate_gradient_method(RnFunction f, double *x, int n, double step, double gold_step, double grad_step, double epsilon, Printer printer)
{
    int i = 0;
    int k = 0;
    int c = 0;

    // Gradient of x
    double* gr = (double*) malloc(sizeof(double) * n);
    // Direction
    double *s  = (double*) malloc(sizeof(double) * n);

	double gr0_mod = 0.0;
    double gr1_mod = 0.0;
    double gr2_mod = 0.0;
	double grad_norm = 0.0;
	double dstnc = 0.0;
	double alpha = 0.0;
    do
    {
        // Gradient of objectiv function in current point
        gradient(f, x, n, grad_step, gr);
		
		grad_norm = vertor_norm(gr, n);
        if (grad_norm < epsilon) break;

		c++;
		
        // Module of gradient
        gr0_mod = 0.0;
        for (i=0; i<n; i++) gr0_mod = gr0_mod + gr[i]*gr[i];
		//gr0_mod = sqrt(gr0_mod);
		
        // First iteration
        if (k == 0)
        {
			gr1_mod = gr0_mod;
            // First direction is antigradient
            for (i=0; i<n; i++) s[i] = -gr[i];
        }
        else
        {
            gr2_mod = gr0_mod;
            double w = gr2_mod / gr1_mod;
            gr1_mod = gr2_mod;
            // Direction in next (k+1) iteration
            for (i=0; i<n; i++) s[i] = -gr[i] + s[i] * w;
        }
		
        double sn = 0.0;
        for (i=0; i<n; i++) sn = sn + s[i]*s[i];
        sn = 1.0;//sqrt(sn);
        for (i=0; i<n; i++) s[i] = s[i] / sn;
		
		alpha = minimize2(f, x, s, n, step, gold_step);
		step /= 1.2;

        if (printer != NULL) printer(f, x, n, c-1, 0, gr, s, sn, alpha);

		dstnc = 0.0;
        for (i=0; i<n; i++)
        {
            x[i] = x[i] + alpha * s[i];
			dstnc = dstnc + (alpha * s[i])*(alpha * s[i]);
        }
		dstnc = sqrt(dstnc);

        if ( k == n ) { k = 0; } else { k++; }

    } while ( dstnc > epsilon );

	if (printer != NULL) printer(f, x, n, c-1, 0, gr, s, 0.0, alpha);

    free(gr);
    free(s);
}
