#include "methods.h"

double RungaKutta(R2Function f, double y0, double x0, double x, double h)
{
    while (x0 <= x)
    {
        double k1 = f(x0, y0);
        double k2 = f(x0+h/2.0, y0+(h/2.0)*k1);
        double k3 = f(x0+h/2.0, y0+(h/2.0)*k2);
        double k4 = f(x0+h, y0+h*k3);

        y0 = y0 + (h/6) * (k1 + 2*k2 + 2*k3 + k4);
        x0 = x0 + h;
    }
    return y0;
}

void RungaKuttaSystem(RmFunction *f, double x0, const double *y0, double x, double *y, const int n, double h)
{
	h = 0.000001;
	memcpy(y, y0, sizeof(double)*n);
    if (fabs(x-x0) < fabs(h)) return;

    double *k1 = (double*) malloc( sizeof(double) * n );
    double *k2 = (double*) malloc( sizeof(double) * n );
    double *k3 = (double*) malloc( sizeof(double) * n );
    double *k4 = (double*) malloc( sizeof(double) * n );
    double *yc = (double*) malloc( sizeof(double) * n );

    if (x0 > x) h = -fabs(h);

    while (fabs(x-x0)>=(fabs(h)))
    {
        int i=0;
        // Calculating k1 vector
        for (i = 0; i<n; i++) k1[i] = f[i](x0, y, n);

        // Calculating k2 vector
        for (i = 0; i<n; i++) yc[i] = y[i] + (h/2.0) * k1[i];
        for (i = 0; i<n; i++) k2[i] = f[i](x0+h/2.0, yc, n);

        // Calculating k3 vector
        for (i = 0; i<n; i++) yc[i] = y[i] + (h/2.0) * k2[i];
        for (i = 0; i<n; i++) k3[i] = f[i](x0+h/2.0, yc, n);

        // Calculating k1 vector
        for (i = 0; i<n; i++) yc[i] = y[i] + h * k3[i];
        for (i = 0; i<n; i++) k4[i] = f[i](x0+h, yc, n);

        // Calculating y
        for (i = 0; i<n; i++) y[i] = y[i] + (h/6.0) * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);

        x0 = x0 + h;
    }

    free(yc);
    free(k1);
    free(k2);
    free(k3);
    free(k4);
}

/**
 * @breaf Метод Эйлера. Системой дифференциальных уравнений
 * @param f
 * @param x независимый аргумент
 * @param y0
 * @param x
 * @param h
 * @return
 */
double EulerMethod(R2Function f, double x0, double y0, double x, double h)
{
    if (fabs(x-x0) < fabs(h)) return y0;

    if (x0 > x) h = -h;

    while (fabs(x-x0) > fabs(h))
    {
        y0 = y0 + h * f(x0, y0);
        x0 = x0 + h;
    }
    return y0;
}

/**
 * @breaf Метод Эйлера.
 * @param f
 * @param x0 независимый аргумент
 * @param y0
 * @param x
 * @param y
 * @param n
 * @param h
 * @return
 */
void EulerMethodSystem(RmFunction *f, double x0, const double *y0, double x, double *y, int n, double h)
{
    memcpy(y, y0, sizeof(double)*n);
    if (fabs(x-x0) < fabs(h)) return;

    if (x0 > x) h = -h;
    int i=0;
    while (fabs(x-x0) > fabs(h))
    {
        double yc[n];
        memcpy(yc, y, sizeof(double)*n);
        for (i=0; i<n; i++)
            y[i] = yc[i] + h*f[i](x0, yc, n);
        x0 = x0 + h;
    }
}

double integeral_trapezoidal_rule1(double *y, double *x, int n)
{
    int i=0;
    double sum = 0.0;

    for (i=0; i<(n-1); i++)
    {
        sum += (( y[i+1] + y[i] ) * ( x[i+1] - x[i] )) / 2.0;
    }
    return sum;
}

double integeral_trapezoidal_rule2(R1Function f, double *x, int n)
{
    int i=0;
    double sum = 0.0;

    for (i=0; i<(n-1); i++)
    {
        sum += (( f(x[i+1]) + f(x[i]) ) * ( x[i+1] - x[i] )) / 2.0;
    }
    return sum;
}

