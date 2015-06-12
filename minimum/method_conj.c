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

        if (printer != NULL) printer(f, x, n, c-1, 0, gr, s, gr, gr);

		dstnc = 0.0;
        for (i=0; i<n; i++)
        {
            x[i] = x[i] + alpha * s[i];
			dstnc = dstnc + (alpha * s[i])*(alpha * s[i]);
        }
		dstnc = sqrt(dstnc);

        if ( k == n ) { k = 0; } else { k++; }

    } while ( dstnc > epsilon );

	if (printer != NULL) printer(f, x, n, c-1, 0, gr, s, gr, gr);

    free(gr);
    free(s);
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
void conjugate_gradient_method1(RnFunction f, double *x, int n, double line_step, double gold_step, double grad_step, double epsilon, Printer printer)
{
    int i = 0;
    int k = 0;
    int iter = 0;
    double mod_s = 0.0;
    double dist = 0.0;
    double *s  = (double*) malloc(sizeof(double) * n);
    double *s1 = (double*) malloc(sizeof(double) * n);
    double *x1 = (double*) malloc(sizeof(double) * n);
    double *x2 = (double*) malloc(sizeof(double) * n);
    int count = 0;
    
    for (i=0; i<n; i++)
    {
        s[i]  = 0.0;
        s1[i] = 0.0;
    }

    double* gr1 = (double*) malloc(sizeof(double) * n);
    double* gr2 = (double*) malloc(sizeof(double) * n);

    double sn = 0.0;
    double gr1_mod = 0.0;
    double gr2_mod = 0.0;
    do
    {
        // First iteration
        if (k == 0)
        {
            // First direction is gradient direction
            gradient(f, x, n, grad_step, gr1);

            for (i=0; i<n; i++) s[i] = -gr1[i];

            // Norm of direction
            sn = vertor_norm(s, n);

            // Divide direction to its norm
            for (i=0; i<n; i++) s1[i] = s[i] / sn;

            // Module of gradient
            gr1_mod = 0.0;
            for (i=0; i<n; i++) gr1_mod = gr1_mod + gr1[i]*gr1[i];
        }
        else
        {
            // Calculating gradient in next coordinates
            gradient(f, x, n, grad_step, gr2);

            // Module of next gradient
            gr2_mod = 0.0;
            for (i=0; i<n; i++) gr2_mod = gr2_mod + gr2[i]*gr2[i];

            double w = gr2_mod / gr1_mod;
            gr1_mod = gr2_mod;

            // Calculating direction for next (k+1) iteration
            for (i=0; i<n; i++) s[i] = -gr2[i] + s[i] * w;

            // Norm of direction
            sn = vertor_norm(s, n);

            // Divide direction to its module
            for (i=0; i<n; i++) s1[i] = s[i] / sn;
        }

        if (printer != NULL) printer(f, x, n, iter, count, s, s1, gr1, gr2);
        iter++;

        memcpy(x1, x, sizeof(double) * n);

        // Minimization in one dimensional direction
        double argmin(double alpha)
        {
            int j;
            for (j=0; j<n; j++) x2[j] = x[j] + alpha * s1[j];
            double result = f(x2, n);
            return result;
        }

        double a,b;
        double alpha0 = 0.0;
        straight_line_search_metod(argmin, alpha0, line_step, &a, &b);
        double alpha = golden_section_search_min(argmin, a, b, gold_step);
        //double alpha = minimize(f, x, s1, n, alpha0, line_step, gold_step);
        //line_step /= 1.2;

        if (argmin(alpha)>argmin(alpha0)) alpha = alpha0;

        // Calculating next coordinates
        for (i=0; i<n; i++)
        {
            x[i] = x[i] + alpha * s1[i];
        }

        mod_s = 0.0;
        //for (i=0; i<n; i++) mod_s = mod_s + s[i]*s[i];
        mod_s = vertor_norm(s, n);
        dist = distance(x1, x, n);

        if ( k == n ) { k = 0; } else { k++; }

    } while ( mod_s > epsilon && dist > epsilon );

    if (printer != NULL) printer(f, x, n, iter, count, s, s1, gr1, gr2);

    free(gr1);
    free(gr2);
    free(s1);
    free(s);
    free(x1);
    free(x2);

    gr1 = gr2 = s1 = s = NULL;
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
void conjugate_gradient_method2(RnFunction f, double *x, int n, double line_step, double gold_step, double grad_step, double epsilon, Printer printer)
{
    int iter = 0;
    int i = 0;
    double *s = (double*) malloc(sizeof(double) * n);
    double *s1 = (double*) malloc(sizeof(double) * n);
    int count = 0;
    
    for (i=0; i<n; i++)
    {
        s[i] = 0.0;
        s1[i] = 0.0;
    }

    do
    {
        // First iteration
        double* gr1 = (double*) malloc(sizeof(double) * n);
        gradient(f, x, n, grad_step, gr1);

        for (i=0; i<n; i++)
            s[i] = -gr1[i];

        double ss = vertor_norm(s, n);

        for (i=0; i<n; i++)
            s1[i] = s[i] / ss;

        double gr1_mod = 0.0;
        for (i=0; i<n; i++)
            gr1_mod += gr1[i]*gr1[i];

        int k = 0;
        do
        {
            if (printer != NULL) printer(f, x, n, iter, count, s, s1);

            double alpha = minimize2(f, x, s1, n, line_step, gold_step);
            line_step /= 1.2;

            for (i=0; i<n; i++)
            {
                x[i] = x[i] + alpha * s1[i];
            }

            double* gr2 = (double*) malloc(sizeof(double) * n);
            gradient(f, x, n, grad_step, gr2);
            count += 2*n+1;

            double gr2_mod = 0.0;
            for (i=0; i<n; i++)
                gr2_mod += gr2[i]*gr2[i];

            double w = gr2_mod / gr1_mod;
            gr1_mod = gr2_mod;

            for (i=0; i<n; i++) s[i] = -gr2[i] + s[i] * w;

            ss = vertor_norm(s, n);
            for (i=0; i<n; i++) s1[i] = s[i] / ss;

            free(gr2);
            gr2 = NULL;

            iter++;
            k++;
        } while ( k < n );
        
        free(gr1);
        gr1 = NULL;
        
        double mod_s = s[0]*s[0] + s[1]*s[1];

        if ( mod_s < epsilon ) break;
    } while ( 1 );
    
    free(s1);
    s1 = NULL;
    free(s);
    s = NULL;
}
