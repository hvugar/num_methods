#include "cmethods.h"

double derivative(CR1Function fx)
{
    return fx(5.0);
}

double derivative1(CR1Function fx, double x, double h)
{
    return (fx(x+h) - fx(x)) / h;
}

double derivative2(CR1Function fx, double x, double h)
{
    return (fx(x) - fx(x-h)) / h;
}

double derivative3(CR1Function fx, double x, double h)
{
    return (fx(x+h) - fx(x-h)) / (2.0*h);
}

void gradient1(CRnFunction fx, double *x, double *g, unsigned int n, double h)
{
    unsigned int i;
    double tx;
    double f0 = fx(x, n);
    for (i=0; i<n; i++)
    {
        tx = x[i];
        x[i] += h;
        double f = fx(x, n);
        g[i] = (f-f0)/h;
        x[i] = tx;
    }
}

void gradient2(CRnFunction fx, double *x, double *g, unsigned int n, double h)
{
    unsigned int i;
    double tx;
    double f0 = fx(x, n);
    for (i=0; i<n; i++)
    {
        tx = x[i];
        x[i] -= h;
        double f = fx(x, n);
        g[i] = (f-f0)/h;
        x[i] = tx;
    }
}

void gradient3(CRnFunction fx, double *x, double *g, unsigned int n, double h)
{
    unsigned int i;
    double tx;
    for (i=0; i<n; i++)
    {
        tx = x[i];
        x[i] = tx + h;
        double f1 = fx(x, n);
        x[i] = tx - h;
        double f2 = fx(x, n);
        g[i] = (f2-f1)/(2.0*h);
        x[i] = tx;
    }
}

double trapesium1(CR1Function fx, unsigned int n, double a, double b)
{
    double sum = 0.0;
    double h = (a-b)/n;
    unsigned int i;
    for (i=0; i<=n-1; i++)
    {
        double x = a + i*h;
        double f1 = fx(x);
        double f2 = fx(x+h);
        sum = sum + (f1+f2);
    }
    sum = (h/2.0)*sum;
    return sum;
}

double trapesium2(CR1Function fx, double h, double a, double b)
{
    unsigned int n = (unsigned int)round((b - a)/h);
    double sum = 0.0;
    unsigned int i;
    for (i=0; i<=(n-1); i++)
    {
        double x = a + i*h;
        double f1 = fx(x);
        double f2 = fx(x+h);
        sum = sum + (f1+f2);
    }
    sum = (h/2.0)*sum;
    return sum;
}

void tomasAlgorithm(const double *a, const double *b, const double *c, const double *d, double *x, unsigned int n)
{
    double *p = (double*)malloc(sizeof(double)*n);
    double *q = (double*)malloc(sizeof(double)*n);

    const unsigned int j = (unsigned)0-1;
    unsigned int i;
    for (i=0; i<n; i++)
    {
        if (i==0)
        {
            p[0] = +d[0]/b[0];
            q[0] = -c[0]/b[0];
        }
        else
        {
            if(i==(n-1))
            {
                p[n-1] = (d[i]-a[i]*p[i-1])/(b[i]+a[i]*q[i-1]);
                q[n-1] = 0.0;
            }
            else
            {
                p[i] = +(d[i]-a[i]*p[i-1])/(b[i]+a[i]*q[i-1]);
                q[i] = -c[i]/(b[i]+a[i]*q[i-1]);
            }
        }
    }

    for (i=n-1; i != j; i--)
    {
        if (i==(n-1))
        {
            x[i] = p[i];
        }
        else
        {
            x[i] = p[i] + q[i]*x[i+1];
        }
    }

    free(p);
    free(q);
}

void gaussianElimination(double **A, double *b, double *x, unsigned int n)
{
    unsigned int k;
    for (k=0; k<n-1; k++)
    {
        unsigned int j;
        unsigned int i;

        for (j=(k+1); j<n; j++)
        {
            if (A[k][k] != 0.0)
            {
                double f = A[j][k]/A[k][k];
                for (i=k; i<n; i++) A[j][i] += A[k][i] * f;
                b[j] += b[k] *f;
            }
        }
    }

    //    int k;
    //    for (k=0; k<n-1; k++)
    //    {
    //        int i,j;
    //        for (i=(k+1); i<n; i++)
    //        {
    //            if (a[k][k] != 0.0)
    //            {
    //                double f = a[i][k] / a[k][k];
    //                for (j=k; j<n; j++)
    //                {
    //                    a[i][j] = a[i][j] - a[k][j] * f;
    //                }
    //                b[i] = b[i] - b[k] * f;
    //            }
    //            else
    //            {
    //                // printf("%d %f\n", k, a[k][k]);
    //            }
    //            check_matrix(a, b, n);
    //        }
    //    }

    //    int i;
    //    for (i=(n-1); i>=0; i--)
    //    {
    //        int j;
    //        for (j=(n-1); j>i; j--) b[i] -= (a[i][j] * x[j]);
    //        x[i] = b[i] / a[i][i];
    //    }
}

void gaussJordanElimination(double **a, double *b, double *x, unsigned int n)
{
}

void euler(double x0, double y0, double xN, double yN, unsigned int N, double *x, double *y, ODE1stOrderEquation eq)
{

}

void eulerMod(double x0, double y0, double xN, double yN, unsigned int N, double *x, double *y, ODE1stOrderEquation eq)
{

}

void runge_kutta_rk3(double x0, double y0, double xN, double yN, unsigned int N, double *x, double *y, ODE1stOrderEquation eq)
{

}

void runge_kutta_rk4(double x0, double y0, double xN, double yN, unsigned int N, double *x, double *y, ODE1stOrderEquation eq)
{
    if (N == 0) return;

    double k1 = 0.0;
    double k2 = 0.0;
    double k3 = 0.0;
    double k4 = 0.0;

    double h = (xN - x0) / N;
    //    x = (double*)malloc(sizeof(double)*(N+1));
    //    y = (double*)malloc(sizeof(double)*(N+1));

    if (h > 0)
    {
        x[0] = x0;
        y[0] = y0;
        int i = 1;
        for (i=1; i<=N; i++)
        {
            k1 = eq(x0, y0);
            k2 = eq(x0+h/2.0, y0+(h/2.0)*k1);
            k3 = eq(x0+h/2.0, y0+(h/2.0)*k2);
            k4 = eq(x0+h, y0+h*k3);
            y0 = y0 + (h/6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);
            x0 = x0 + h;
            x[i] = x0;
            y[i] = y0;
        }
        yN = y[N];
    }

    if (h < 0)
    {
        x[N] = x0;
        y[N] = y0;
        int i = N-1;
        for (i=N-1; i>=0; i--)
        {
            k1 = eq(x0,       y0);
            k2 = eq(x0+h/2.0, y0+(h/2.0)*k1);
            k3 = eq(x0+h/2.0, y0+(h/2.0)*k2);
            k4 = eq(x0+h,     y0+h*k3);
            y0 = y0 + (h/6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);
            x0 = x0 + h;
            x[i] = x0;
            y[i] = y0;
        }
        yN = y[0];
    }
}

void runge_kutta_rk4_system(double x0, double x1, double *y0, double **y, size_t n, unsigned int N, double h, ODE1stOrderEquationN *eq)
{
    double *k1 = (double *)malloc(sizeof(double) * n);
    double *k2 = (double *)malloc(sizeof(double) * n);
    double *k3 = (double *)malloc(sizeof(double) * n);
    double *k4 = (double *)malloc(sizeof(double) * n);

    unsigned int i, j;
    for (i=0; i<n; i++) y[i][0] = y0[i];

    for (i=0; i<N; i++)
    {
        double *arg = (double*) malloc(sizeof(double) * n);

        // Calculating k1 vector
        for (j = 0; j<n; j++) arg[j] = y[j][i];
        for (j = 0; j<n; j++) k1[j] = eq[j](x0, arg, n);

        // Calculating k2 vector
        for (j = 0; j<n; j++) arg[j] = y[j][i] + (h/2.0) * k1[j];
        for (j = 0; j<n; j++) k2[j] = eq[j](x0+h/2.0, arg, n);

        // Calculating k3 vector
        for (j = 0; j<n; j++) arg[j] = y[j][i] + (h/2.0) * k2[j];
        for (j = 0; j<n; j++) k3[j] = eq[j](x0+h/2.0, arg, n);

        // Calculating k4 vector
        for (j = 0; j<n; j++) arg[j] = y[j][i] + h * k3[j];
        for (j = 0; j<n; j++) k4[j] = eq[j](x0+h, arg, n);

        // Calculating y
        for (j = 0; j<n; j++) y[j][i+1] = y[j][i] + (h/6.0) * (k1[j] + 2*k2[j] + 2*k3[j] + k4[j]);

        x0 = x0 + h;

        free(arg);
    }

    free(k1);
    free(k2);
    free(k3);
    free(k4);
}


//void RungaKuttaSystem(RmFunction *f, double x0, const double *y0, double x, double *y, const int n, double h)
//{
//    //h = 0.000001;
//    memcpy(y, y0, sizeof(double)*n);
//    if (fabs(x-x0) < fabs(h)) return;

//    double *k1 = (double*) malloc( sizeof(double) * n );
//    double *k2 = (double*) malloc( sizeof(double) * n );
//    double *k3 = (double*) malloc( sizeof(double) * n );
//    double *k4 = (double*) malloc( sizeof(double) * n );
//    double *yc = (double*) malloc( sizeof(double) * n );

//    if (x0 > x) h = -fabs(h);

//    while (fabs(x-x0)>=(fabs(h)))
//    {
//        int i=0;
//        // Calculating k1 vector
//        for (i = 0; i<n; i++) k1[i] = f[i](x0, y, n);

//        // Calculating k2 vector
//        for (i = 0; i<n; i++) yc[i] = y[i] + (h/2.0) * k1[i];
//        for (i = 0; i<n; i++) k2[i] = f[i](x0+h/2.0, yc, n);

//        // Calculating k3 vector
//        for (i = 0; i<n; i++) yc[i] = y[i] + (h/2.0) * k2[i];
//        for (i = 0; i<n; i++) k3[i] = f[i](x0+h/2.0, yc, n);

//        // Calculating k1 vector
//        for (i = 0; i<n; i++) yc[i] = y[i] + h * k3[i];
//        for (i = 0; i<n; i++) k4[i] = f[i](x0+h, yc, n);

//        // Calculating y
//        for (i = 0; i<n; i++) y[i] = y[i] + (h/6.0) * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);

//        x0 = x0 + h;
//    }

//    free(yc);
//    free(k1);
//    free(k2);
//    free(k3);
//    free(k4);
//}
