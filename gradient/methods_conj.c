#include "methods.h"
#include "print.h"

double minimize(RnFunction f, double *x, double *grad, int n, double alpha0, double line_step, double gold_step);

/**
 * @brief Метод сопряженных градиентов Флетчера — Ривса
 * @param f Целевая функция
 * @param x Независимые переменные n - измерение
 * @param n
 * @param line_eps
 * @param gold_eps
 * @param grad_eps
 * @param epsilon
 */
void conjugate_gradient_method(RnFunction f, double *x, int n, double line_step, double gold_step, double grad_step, double epsilon, Printer printer)
{
    int i = 0;
    int k = 0;
	
	int iter = 0;
    int count = 0;

	// Direction
	double *s  = (double*) malloc(sizeof(double) * n);
	// Saves last point coordinates
    double *x1 = (double*) malloc(sizeof(double) * n);
	// Used for one dimention minimization for stopring next point coordinates
    double *x2 = (double*) malloc(sizeof(double) * n);

	// Gradient of x(k) point
    double* gr1 = (double*) malloc(sizeof(double) * n);
	// Gradinet of x(k+1) point
    double* gr2 = (double*) malloc(sizeof(double) * n);

    double gr1_mod = 0.0;
    double gr2_mod = 0.0;
    do
    {
        // First iteration
        if (k == 0)
        {
            // Gradient of objectiv function in current point
            gradient(f, x, n, grad_step, gr1);

			// First direction is antigradient
            for (i=0; i<n; i++) s[i] = -gr1[i];

            // Module of gradient
            gr1_mod = 0.0;
            for (i=0; i<n; i++) gr1_mod += gr1[i]*gr1[i];
        }
        else
        {
            /// Gradient of objectiv function in next point
            gradient(f, x, n, grad_step, gr2);

            // Module of next gradient
            gr2_mod = 0.0;
            for (i=0; i<n; i++) gr2_mod = gr2_mod + gr2[i]*gr2[i];

            double w = gr2_mod / gr1_mod;
            gr1_mod = gr2_mod;

            // Direction in next (k+1) iteration
            for (i=0; i<n; i++) s[i] = -gr2[i] + s[i] * w;
        }

		if (printer != NULL) printer(f, x, n, iter, count, s, s, gr1, gr2);
        iter++;

        // Minimization in one dimensional direction
		double argmin(double alpha)
		{
			int j;
			for (j=0; j<n; j++) x2[j] = x[j] + alpha * s[j];
			return f(x2, n);
		}
    
		double a,b;
		double alpha0 = 0.0;
		straight_line_search_metod(argmin, alpha0, line_step, &a, &b);
		double alpha = golden_section_search_min(argmin, a, b, gold_step);
        //double alpha = minimize(f, x, s1, n, alpha0, line_step, gold_step);
        //line_step /= 1.2;
		
		//
		if (argmin(alpha)>argmin(alpha0)) alpha = alpha0;
	
		// Saving last point coordinates
		memcpy(x1, x, sizeof(double) * n);
		
        // Calculating next point coordinates
        for (i=0; i<n; i++)
        {
            x[i] = x[i] + alpha * s[i];
        }

        //mod_s = 0.0;
        //for (i=0; i<n; i++) mod_s = mod_s + s[i]*s[i];
		//mod_s = vertor_norm(s, n);
		//dist = distance(x1, x, n);
		
        if ( k == n ) { k = 0; } else { k++; }
		
    } while ( vertor_norm(s, n) > epsilon && distance(x1, x, n) > epsilon );
	
	if (printer != NULL) printer(f, x, n, iter, count, s, s, gr1, gr2);

    free(gr1);
    free(gr2);
    free(s);
	free(x1);
	free(x2);

    gr1 = gr2 = s = NULL;
}
 
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

            double alpha0 = 0.0;
            double alpha = minimize(f, x, s1, n, alpha0, line_step, gold_step);
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

double minimize(RnFunction f, double *x, double *s, int n, double alpha0, double line_step, double gold_step)
{
    double argmin(double alpha)
    {
        int j;
        for (j=0; j<n; j++) x[j] = x[j] + alpha * s[j];
        double result = f(x, n);
	    for (j=0; j<n; j++) x[j] = x[j] - alpha * s[j];
        return result;
    }
    
    double a,b;
    straight_line_search_metod(argmin, alpha0, line_step, &a, &b);
    double min = golden_section_search_min(argmin, a, b, gold_step);
	
	if (argmin(min)>argmin(alpha0)) 
	{
		fprintf(stderr, "Value of function in min point is greater than value of initial point...\n");
		fprintf(stderr, "min  := %.10f f(x) := %.10f\ninit := %.10f f(x) := %.10f\n", min, argmin(min), alpha0, argmin(alpha0));
		system("pause");
		min = alpha0;
	}
	
    return min;
}
