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

double Gradient::minimize()
{
    class ArgMin : public R1Minimize
    {
    public:
        int n;
        std::vector<double> x0;
        std::vector<double> x1;
        double *grads;
        Gradient *gradient;

        virtual double fx(double alpha)
        {
            for (int i=0; i<n; i++) x1[i] = x0[i] - alpha * grads[i];
            return gradient->fx(x1);
        }
    };

    ArgMin argmin;
    argmin.n = n;
    argmin.x0 = x0;
    argmin.grads = grads;
    argmin.gradient = this;

    double alpha0 = 0.0;
    argmin.setX0(alpha0);
    argmin.setStep(0.1);
    argmin.setEpsilon(0.0001);

    argmin.straightLineSearch();
    double alpha = argmin.goldenSectionSearch();
    if ( argmin.fx(alpha) > argmin.fx(alpha0) ) alpha = alpha0;
    return alpha;
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

        alpha = minimize();

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
