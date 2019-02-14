#include "nlode1o.h"
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

double NonLinearODE1stOrder::f(double, double, unsigned int) const { return NAN; }

double NonLinearODE1stOrder::f(double, const DoubleVector &, unsigned int, unsigned int) const { return NAN; }

void NonLinearODE1stOrder::cauchyProblem(double x0, double y0, DoubleVector &y, OdeSolverMethod method, Direction direction)
{
    switch (method)
    {
    case OdeSolverMethod::RK2:
        calculateRK2(x0, y0, y, direction);
        break;
    case OdeSolverMethod::RK4:
        calculateRK4(x0, y0, y, direction);
        break;
    case OdeSolverMethod::EULER:
        calculateEuler(x0, y0, y, direction);
        break;
    case OdeSolverMethod::EULER_MOD:
        calculateEulerMod(x0, y0, y, direction);
        break;
    }
}

void NonLinearODE1stOrder::calculateRK2(double x0, double y0, DoubleVector &y, Direction direction)
{
    unsigned int N = mgrid.dimension().size();
    double h = mgrid.dimension().step();

    y.clear();
    y.resize(N+1);

    double k1 = 0.0;
    double k2 = 0.0;

    if (direction == L2R)
    {
        double h2 = h/2.0;
        double xn = x0;
        double yn = y0;

        y[0] = yn;
        for (unsigned int n=1; n<=N; n++)
        {
            k1 = f(xn,    yn,       n-1);
            k2 = f(xn+h2, yn+h2*k1, n-1);
            yn += h*k2;
            xn += h;
            y[n] = yn;
        }
    }

    if (direction == R2L)
    {
        double h2 = h/2.0;
        double xn = x0;
        double yn = y0;

        y[N] = yn;
        for (unsigned int n=N-1; n!=UINT32_MAX; n--)
        {
            k1 = f(xn,    yn,       n+1);
            k2 = f(xn-h2, yn-h2*k1, n+1);
            yn -= h*k2;
            xn -= h;
            y[n] = yn;
        }
    }
}

void NonLinearODE1stOrder::calculateRK4(double x0, double y0, DoubleVector &y, Direction direction)
{
    unsigned int N = mgrid.dimension().size();
    double h = mgrid.dimension().step();

    y.clear();
    y.resize(N+1);

    double k1 = 0.0;
    double k2 = 0.0;
    double k3 = 0.0;
    double k4 = 0.0;

    if (direction == L2R)
    {
        double h2 = h/2.0;
        double h6 = h/6.0;
        double xn = x0;
        double yn = y0;

        y[0] = yn;
        for (unsigned int n=1; n<=N; n++)
        {
            k1 = f(xn,    yn,       n-1);
            k2 = f(xn+h2, yn+h2*k1, n-1);
            k3 = f(xn+h2, yn+h2*k2, n-1);
            k4 = f(xn+h,  yn+h*k3,  n-1);

            yn += h6*(k1+2.0*k2+2.0*k3+k4);
            xn += h;

            y[n] = yn;
        }
    }

    if (direction == R2L)
    {
        double h2 = h/2.0;
        double h6 = h/6.0;
        double xn = x0;
        double yn = y0;

        y[N] = yn;
        for (unsigned int n=N-1; n!=UINT32_MAX; n--)
        {
            k1 = f(xn,    yn,      n+1);
            k2 = f(xn-h2, yn-h2*k1,n+1);
            k3 = f(xn-h2, yn-h2*k2,n+1);
            k4 = f(xn-h,  yn-h*k3, n+1);

            yn -= h6*(k1+2.0*k2+2.0*k3+k4);
            xn -= h;

            y[n] = yn;
        }
    }
}

void NonLinearODE1stOrder::calculateEuler(double x0, double y0, DoubleVector &y, Direction direction)
{
    unsigned int N = mgrid.dimension().size();
    double h = mgrid.dimension().step();

    y.clear();
    y.resize(N+1);

    if (direction == L2R)
    {
        double xn = x0;
        double yn = y0;

        y[0] = yn;
        for (unsigned int n=1; n<=N; n++)
        {
            yn = yn + h*f(xn, yn, n-1);
            y[n] = yn;
            xn += h;
        }
    }

    if (direction == R2L)
    {
        double xn = x0;
        double yn = y0;

        y[N] = yn;
        for (unsigned int n=N-1; n!=UINT32_MAX; n--)
        {
            yn = yn - h*f(xn, yn, n+1);
            y[n] = yn;
            xn -= h;
        }
    }
}

void NonLinearODE1stOrder::calculateEulerMod(double x0, double y0, DoubleVector &y, Direction direction)
{
    unsigned int N = mgrid.dimension().size();
    double h = mgrid.dimension().step();

    y.clear();
    y.resize(N+1);

    if (direction == L2R)
    {
        double xn = x0;
        double yn = y0;

        y[0] = yn;
        for (unsigned int n=1; n<=N; n++)
        {
            yn = yn + h*f(xn, yn, n-1);
            y[n] = yn;
            xn += h;
        }
    }

    if (direction == R2L)
    {
        double xn = x0;
        double yn = y0;

        y[N] = yn;
        for (unsigned int n=N-1; n!=UINT32_MAX; n--)
        {
            yn = yn - h*f(xn, yn, n+1);
            y[n] = yn;
            xn -= h;
        }
    }
}

void NonLinearODE1stOrder::cauchyProblem(double x0, const DoubleVector &y0, std::vector<DoubleVector> &ry, OdeSolverMethod method, Direction direction)
{
    switch (method)
    {
    case RK2:
        calculateRK2(x0, y0, ry, direction);
        break;
    case RK4:
        calculateRK4(x0, y0, ry, direction);
        break;
    case EULER:
        calculateEuler(x0, y0, ry, direction);
        break;
    case EULER_MOD:
        calculateEulerMod(x0, y0, ry, direction);
        break;
    }
}

void NonLinearODE1stOrder::calculateRK2(double x0, const DoubleVector &y0, std::vector<DoubleVector> &ry, Direction direction)
{
    const int min = grid().dimension().min();
    const int max = grid().dimension().max();
    const int N = grid().dimension().size();
    const double h = grid().dimension().step();
    const unsigned int n = y0.length();

    ry.clear();
    ry.resize(n);
    for (unsigned int i=0; i<n; i++) ry[i].resize(N+1);

    double *k1 = (double*) malloc( sizeof(double) * n );
    double *k2 = (double*) malloc( sizeof(double) * n );

    if (direction == L2R)
    {
        double h2 = h/2.0;

        double xn = x0;

        DoubleVector yn(n);

        /* initializing */
        for (unsigned int j=0; j<n; j++) ry[j][0] = y0[j];

        for (unsigned int i=min+1; i<=max; i++)
        {
            // k1
            for (unsigned int j=0; j<n; j++) yn[j] = ry[j][(i-1)-min];
            for (unsigned int j=0; j<n; j++) k1[j] = f(xn, yn, i-1, j);

            // k2
            for (unsigned int j=0; j<n; j++) yn[j] = ry[j][(i-1)-min]+h2*k1[j];
            for (unsigned int j=0; j<n; j++) k2[j] = f(xn+h2, yn, i-1, j);

            for (unsigned int j=0; j<n; j++) ry[j][i-min] = ry[j][(i-1)-min] + h * k2[j];

            xn += h;
        }
    }

    if (direction == R2L)
    {
        double h2 = h/2.0;

        double xn = x0;

        DoubleVector yn(n);

        /* initializing */
        for (unsigned int j=0; j<n; j++) ry[j][N] = y0[j];

        for (unsigned int i=N-1; i!=UINT32_MAX; i--)
        {
            // k1
            for (unsigned int j=0; j<n; j++) yn[j] = ry[j][i+1];
            for (unsigned int j=0; j<n; j++) k1[j] = f(xn, yn, i+1, j);

            // k2
            for (unsigned int j=0; j<n; j++) yn[j] = ry[j][i+1]-h2*k1[j];
            for (unsigned int j=0; j<n; j++) k2[j] = f(xn-h2, yn, i+1, j);

            for (unsigned int j=0; j<n; j++) ry[j][i] = ry[j][i+1] - h * k2[j];

            xn -= h;
        }
    }

    free(k2);
    free(k1);
}

void NonLinearODE1stOrder::calculateRK4(double x0, const DoubleVector &y0, std::vector<DoubleVector> &ry, Direction direction)
{
    const int min = grid().dimension().min();
    const int max = grid().dimension().max();
    const int N = grid().dimension().size();
    const double h = grid().dimension().step();
    const unsigned int n = y0.length();


    ry.clear();
    ry.resize(n);
    for (unsigned int i=0; i<n; i++) ry[i].resize(N+1);

    double *k1 = ( double* ) malloc( sizeof(double) * n );
    double *k2 = ( double* ) malloc( sizeof(double) * n );
    double *k3 = ( double* ) malloc( sizeof(double) * n );
    double *k4 = ( double* ) malloc( sizeof(double) * n );

    if (direction == L2R)
    {
        double h2 = h/2.0;
        double h6 = h/6.0;

        double xn = x0;

        DoubleVector yn(n);

        /* initializing */
        for (unsigned int j=0; j<n; j++) ry[j][0] = y0[j];

        for (unsigned int i=min; i<max; i++)
        {
            // k1
            for (unsigned int j=0; j<n; j++) yn[j] = ry[j][i-min];
            for (unsigned int j=0; j<n; j++) k1[j] = f(xn, yn, i, j);

            // k2
            for (unsigned int j=0; j<n; j++) yn[j] = ry[j][i-min]+h2*k1[j];
            for (unsigned int j=0; j<n; j++) k2[j] = f(xn+h2, yn, i, j);

            // k3
            for (unsigned int j=0; j<n; j++) yn[j] = ry[j][i-min]+h2*k2[j];
            for (unsigned int j=0; j<n; j++) k3[j] = f(xn+h2, yn, i, j);

            // k4
            for (unsigned int j=0; j<n; j++) yn[j] = ry[j][i-min]+h*k3[j];
            for (unsigned int j=0; j<n; j++) k4[j] = f(xn+h, yn, i, j);

            for (unsigned int j=0; j<n; j++) ry[j][i-min+1] = ry[j][i-min] + h6 * (k1[j] + 2*k2[j] + 2*k3[j] + k4[j]);

            xn += h;
        }
    }

    if (direction == R2L)
    {
        double h2 = h/2.0;
        double h6 = h/6.0;

        double xn = x0;

        DoubleVector yn(n);

        /* initializing */
        for (unsigned int j=0; j<n; j++) ry[j][N] = y0[j];

        for (unsigned int k=max; k!=min+0; k--)
        {
            // k1
            for (unsigned int j=0; j<n; j++) yn[j] = ry[j][k-min];
            for (unsigned int j=0; j<n; j++) k1[j] = f(xn, yn, k, j);

            // k2
            for (unsigned int j=0; j<n; j++) yn[j] = ry[j][k-min]-h2*k1[j];
            for (unsigned int j=0; j<n; j++) k2[j] = f(xn-h2, yn, k, j);

            // k3
            for (unsigned int j=0; j<n; j++) yn[j] = ry[j][k-min]-h2*k2[j];
            for (unsigned int j=0; j<n; j++) k3[j] = f(xn-h2, yn, k, j);

            // k4
            for (unsigned int j=0; j<n; j++) yn[j] = ry[j][k-min]-h*k3[j];
            for (unsigned int j=0; j<n; j++) k4[j] = f(xn-h, yn, k, j);

            for (unsigned int j=0; j<n; j++) ry[j][k-min-1] = ry[j][k-min] - h6 * (k1[j] + 2*k2[j] + 2*k3[j] + k4[j]);

            xn -= h;
        }
    }


    free(k4);
    free(k3);
    free(k2);
    free(k1);
}

void NonLinearODE1stOrder::calculateEuler(double x0, const DoubleVector &y0, std::vector<DoubleVector> &ry, Direction direction)
{
    const int min = grid().dimension().min();
    const int max = grid().dimension().max();
    const int N = grid().dimension().size();
    const double h = grid().dimension().step();
    const unsigned int n = y0.length();

    ry.clear();
    ry.resize(n);
    for (unsigned int i=0; i<n; i++) ry[i].resize(N+1);

    if (direction == L2R)
    {
        double xn = x0;
        DoubleVector yn(n);
        /* initializing */
        for (unsigned int j=0; j<n; j++) ry[j][0] = y0[j];

        for (unsigned int i=min+1; i<=max; i++)
        {
            for (unsigned int j=0; j<n; j++) yn[j] = ry[j][(i-1)-min];
            for (unsigned int j=0; j<n; j++) ry[j][i-min] = ry[j][(i-1)-min] + h * f(xn, yn, i-1, j);

            xn += h;
        }
    }

    if (direction == R2L)
    {
        double xn = x0;
        DoubleVector yn(n);
        /* initializing */
        for (unsigned int j=0; j<n; j++) ry[j][N] = y0[j];

        for (unsigned int i=N-1; i!=UINT32_MAX; i--)
        {
            for (unsigned int j=0; j<n; j++) yn[j] = ry[j][i+1];
            for (unsigned int j=0; j<n; j++) ry[j][i] = ry[j][i+1] - h * f(xn, yn, i+1, j);

            xn -= h;
        }
    }
}

void NonLinearODE1stOrder::calculateEulerMod(double x0, const DoubleVector &y0, std::vector<DoubleVector> &ry, Direction direction)
{
    const int min = grid().dimension().min();
    const int max = grid().dimension().max();
    const int N = grid().dimension().size();
    const double h = grid().dimension().step();
    const unsigned int n = y0.length();


    ry.clear();
    ry.resize(n);
    for (unsigned int i=0; i<n; i++) ry[i].resize(N+1);

    if (direction == L2R)
    {
        double xn = x0;
        DoubleVector yn(n);
        /* initializing */
        for (unsigned int j=0; j<n; j++) ry[j][0] = y0[j];

        for (unsigned int i=min+1; i<=max; i++)
        {
            for (unsigned int j=0; j<n; j++) yn[j] = ry[j][(i-1)-min];
            for (unsigned int j=0; j<n; j++) ry[j][i-min] = ry[j][(i-1)-min] + h * f(xn, yn, i-1, j);

            xn += h;
        }
    }

    if (direction == R2L)
    {
        double xn = x0;
        DoubleVector yn(n);
        /* initializing */
        for (unsigned int j=0; j<n; j++) ry[j][N] = y0[j];

        for (unsigned int i=N-1; i!=UINT32_MAX; i--)
        {
            for (unsigned int j=0; j<n; j++) yn[j] = ry[j][i+1];
            for (unsigned int j=0; j<n; j++) ry[j][i] = ry[j][i+1] - h * f(xn, yn, i+1, j);

            xn -= h;
        }
    }
}

void NonLinearODE1stOrder::cauchyProblem(double x0, const DoubleVector &y0, DoubleVector &ry, OdeSolverMethod method, Direction direction)
{
    switch (method)
    {
    case RK2:
        //calculateRK2(x0, y0, y, direction);
        break;
    case RK4:
        calculateRK4(x0, y0, ry, direction);
        break;
    case EULER:
        //calculateEuler(x0, y0, y, direction);
        break;
    case EULER_MOD:
        //calculateEulerMod(x0, y0, y, direction);
        break;
    }
}

void NonLinearODE1stOrder::calculateRK4(double x0, const DoubleVector &y0, DoubleVector &ry, Direction direction)
{
    const int min = grid().dimension().min();
    const int max = grid().dimension().max();
    const int N = grid().dimension().size();
    const double h = grid().dimension().step();
    const unsigned int n = y0.length();

    ry.clear();
    ry.resize(n, N+1);

    double *k1 = (double *)malloc(sizeof(double)*n);
    double *k2 = (double *)malloc(sizeof(double)*n);
    double *k3 = (double *)malloc(sizeof(double)*n);
    double *k4 = (double *)malloc(sizeof(double)*n);

    if (direction == L2R)
    {
        double h2 = h/2.0;
        double h6 = h/6.0;

        double xn = x0;

        DoubleVector yn(n);

        /* initializing */
        for (unsigned int j=0; j<n; j++) ry[j] = y0[j];

        for (unsigned int i=min+1; i<=max; i++)
        {
            // k1
            for (unsigned int j=0; j<n; j++) yn[j] = ry[j];
            for (unsigned int j=0; j<n; j++) k1[j] = f(xn, yn, i-1, j);

            // k2
            for (unsigned int j=0; j<n; j++) yn[j] = ry[j]+h2*k1[j];
            for (unsigned int j=0; j<n; j++) k2[j] = f(xn+h2, yn, i-1, j);

            // k3
            for (unsigned int j=0; j<n; j++) yn[j] = ry[j]+h2*k2[j];
            for (unsigned int j=0; j<n; j++) k3[j] = f(xn+h2, yn, i-1, j);

            // k4
            for (unsigned int j=0; j<n; j++) yn[j] = ry[j]+h*k3[j];
            for (unsigned int j=0; j<n; j++) k4[j] = f(xn+h, yn, i-1, j);

            for (unsigned int j=0; j<n; j++) ry[j] += h6 * (k1[j] + 2*k2[j] + 2*k3[j] + k4[j]);

            xn += h;
        }
    }

    if (direction == R2L)
    {
        double h2 = h/2.0;
        double h6 = h/6.0;

        double xn = x0;

        DoubleVector yn(n);

        /* initializing */
        for (unsigned int j=0; j<n; j++) ry[j] = y0[j];

        for (unsigned int i=N-1; i!=UINT32_MAX; i--)
        {
            // k1
            for (unsigned int j=0; j<n; j++) yn[j] = ry[j];
            for (unsigned int j=0; j<n; j++) k1[j] = f(xn, yn, i+1, j);

            // k2
            for (unsigned int j=0; j<n; j++) yn[j] = ry[j]-h2*k1[j];
            for (unsigned int j=0; j<n; j++) k2[j] = f(xn-h2, yn, i+1, j);

            // k3
            for (unsigned int j=0; j<n; j++) yn[j] = ry[j]-h2*k2[j];
            for (unsigned int j=0; j<n; j++) k3[j] = f(xn-h2, yn, i+1, j);

            // k4
            for (unsigned int j=0; j<n; j++) yn[j] = ry[j]-h*k3[j];
            for (unsigned int j=0; j<n; j++) k4[j] = f(xn-h, yn, i+1, j);

            for (unsigned int j=0; j<n; j++) ry[j] -= h6 * (k1[j] + 2*k2[j] + 2*k3[j] + k4[j]);

            xn -= h;
        }
    }

    free(k4);
    free(k3);
    free(k2);
    free(k1);
}

