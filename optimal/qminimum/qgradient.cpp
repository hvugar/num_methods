#include "qgradient.h"


QGradient::QGradient(RnFunction f, double *x0, int n)
{
    x = (double*) malloc( sizeof(double) * n );
    memcpy(x, x0, n);
    this->fn = f;
    this->n = n;
}

QGradient::~QGradient()
{
    free(x);
}

double QGradient::minimize(double alpha)
{
    Q_UNUSED(alpha);
    double *x1 = (double*) malloc( sizeof(double) * n );
    //int j;
    //for (j=0; j<n; j++) x1[j] = x[j] - alpha * grads[j];
    double min = fn(x1, n);
    free(x1);
    return min;
}

void QGradient::calculate()
{
    int i = 0;
    //double grad_norm = 0.0;

    double* grads = (double*) malloc(sizeof(double) * n);
    // Saves last point coordinates
    double *x1 = (double*) malloc(sizeof(double) * n);
    // Used for one dimention minimization for stopring next point coordinates
    double *x2 = (double*) malloc(sizeof(double) * n);
    //do
    {
        //gradient(f, x, n, grad_eps, grads);

//        double argmin(double alpha)
//        {
//            int j;
//            for (j=0; j<n; j++) x2[j] = x[j] - alpha * grads[j];
//            return f(x2, n);
//        }

//        double a,b;
//        double alpha0 = 0.0;
        //double alpha = 0.0;
//        straight_line_search_metod(argmin, alpha0, line_eps, &a, &b);
//        double alpha = golden_section_search_min(argmin, a, b, gold_eps);
//        if ( argmin(alpha) > argmin(alpha0) ) alpha = alpha0;

        //grad_norm = vertor_norm(grads, n);

        memcpy(x1, x, sizeof(double) * n);
        int j;
        for (j=0; j<n; j++)
        {
        //    x[j] = x[j] - alpha * grads[j];
        }
        i++;
    }
    //while ( grad_norm > epsilon && distance(x1, x, n) > epsilon /*&& fabs(f(x1,n) - f(x,n)) > epsilon*/ );
    free(x1);
    free(x2);
    free(grads);
}

double* QGradient::getX() const
{
    return x;
}
