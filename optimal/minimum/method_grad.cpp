#include "methods.h"

void fast_proximal_gradient_method(RnFunction f, double *x, int n, double line_eps, double gold_eps, double grad_eps, double epsilon, Printer printer)
{
    int i = 0;
    double grad_norm = 0.0;

    double* grads = (double*) malloc(sizeof(double) * n);
    // Saves last point coordinates
    double *x1 = (double*) malloc(sizeof(double) * n);
    // Used for one dimention minimization for stopring next point coordinates
    double *x2 = (double*) malloc(sizeof(double) * n);
    do
    {
        gradient(f, x, n, grad_eps, grads);

        // Funksiyanin minimumuniu tapmaq ucun qizil bolgu qaydasinda istifade edib alphani tapiriq
        //double alpha = minimize1(f, x, grads, n, alpha0, line_eps, gold_eps);
        double argmin(double alpha)
        {
            int j;
            for (j=0; j<n; j++) x2[j] = x[j] - alpha * grads[j];
            return f(x2, n);
        }

        double a,b;
        double alpha0 = 0.0;
        straight_line_search_metod(argmin, alpha0, line_eps, &a, &b);
        double alpha = golden_section_search_min(argmin, a, b, gold_eps);

        if ( argmin(alpha) > argmin(alpha0) ) alpha = alpha0;

        grad_norm = vertor_norm(grads, n);
        if (printer != NULL) printer(f, x, n, i, grads, grad_norm, alpha);

        memcpy(x1, x, sizeof(double) * n);
        int j;
        for (j=0; j<n; j++)
        {
            x[j] = x[j] - alpha * grads[j];
        }
        i++;
    }
    while ( grad_norm > epsilon && distance(x1, x, n) > epsilon /*&& fabs(f(x1,n) - f(x,n)) > epsilon*/ );
    if (printer != NULL) printer(f, x, n, i, grads, grad_norm, 0);
    free(x1);
    free(x2);
    free(grads);
}

double minimize1(RnFunction f, double *x, double *grad, int n, double alpha0, double line_eps, double gold_eps)
{
    double argmin(double alpha)
    {
        int j;
        for (j=0; j<n; j++) x[j] = x[j] - alpha * grad[j];
        double result = f(x, n);
        for (j=0; j<n; j++) x[j] = x[j] + alpha * grad[j];
        return result;
    }
    double a,b;
    straight_line_search_metod(argmin, alpha0, line_eps, &a, &b);
    double min = golden_section_search_min(argmin, a, b, gold_eps);
    return min;
}
