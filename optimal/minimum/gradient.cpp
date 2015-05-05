#include "gradient.h"

Gradient::Gradient()
{
}

void Gradient::conjugate_gradient_method(RnFunction *f, double *x, int n, double line_step, double gold_step, double grad_step, double epsilon)
{
}

void Gradient::fast_proximal_gradient_method(RnFunction *f, double *x, int n, double line_step, double gold_step, double grad_step, double epsilon)
{
}

double Gradient::argmin(double alpha)
{
    double *x2;
    for (int i=0; i<n; i++)
        x2[i] = x[i] - alpha * grads[i];
    return fx(x2, n);
}

double Gradient::R1Minimize()
{
    //    double a,b;
    //    double alpha0 = 0.0;
    //    straight_line_search_metod(argmin, alpha0, line_eps, &a, &b);
    //    double alpha = golden_section_search_min(argmin, a, b, gold_eps);
    //    if ( argmin(alpha) > argmin(alpha0) ) alpha = alpha0;
}

void Gradient::fast_proximal_gradient_method()
{
    double grad_norm = 0.0;
    grads = (double*) malloc(sizeof(double) * n);
    // Saves last point coordinates
    double *x1 = (double*) malloc(sizeof(double) * n);
    // Used for one dimention minimization for stopring next point coordinates
    double *x2 = (double*) malloc(sizeof(double) * n);
    do
    {
        //gradient(f, x, n, grad_eps, grads);
        //grad_norm = vertor_norm(grads, n);
        //memcpy(x1, x, sizeof(double) * n);

        alpha = R1Minimize();

        show();

        for (int i=0; i<n; i++)
        {
            x[i] = x[i] - alpha * grads[i];
        }
        count++;
    }
    while ( true/*grad_norm > epsilon && distance(x1, x, n) > epsilon*/ );

    free(x1);
    free(x2);
    free(grads);
}

void Gradient::show()
{
    printf("%d\n", count);
}
