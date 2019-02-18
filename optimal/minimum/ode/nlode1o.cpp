#include "nlode1o.h"
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

double FirstOrderNonLinearODE::f(double, double, unsigned int) const { return NAN; }

double FirstOrderNonLinearODE::f(double, const DoubleVector &, unsigned int, unsigned int) const { return NAN; }

void FirstOrderNonLinearODE::cauchyProblem(double x0, double y0, DoubleVector &y, OdeSolverMethod method, Direction direction)
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

void FirstOrderNonLinearODE::calculateRK2(double x0, double y0, DoubleVector &y, Direction direction)
{
    int N = dimension().size();
    double h = dimension().step();

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

void FirstOrderNonLinearODE::calculateRK4(double x0, double y0, DoubleVector &y, Direction direction)
{
    int N = dimension().size();
    double h = dimension().step();

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

void FirstOrderNonLinearODE::calculateEuler(double x0, double y0, DoubleVector &y, Direction direction)
{
    unsigned int N = dimension().size();
    double h = dimension().step();

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

void FirstOrderNonLinearODE::calculateEulerMod(double x0, double y0, DoubleVector &y, Direction direction)
{
    unsigned int N = dimension().size();
    double h = dimension().step();

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

void FirstOrderNonLinearODE::cauchyProblem(double x0, const DoubleVector &y0, std::vector<DoubleVector> &ry, OdeSolverMethod method, Direction direction)
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

void FirstOrderNonLinearODE::calculateRK2(double x0, const DoubleVector &y0, std::vector<DoubleVector> &ry, Direction direction)
{
    const int min = dimension().min();
    const int max = dimension().max();
    const unsigned int N = static_cast<unsigned int>( dimension().size() );
    const double h = dimension().step();
    const unsigned int m = count();

    ry.clear(); ry.resize(N+1); for (unsigned int n=0; n<=N; n++) ry[n].resize(m);

    double *k1 = static_cast<double*>( malloc(sizeof(double)*m) );
    double *k2 = static_cast<double*>( malloc(sizeof(double)*m) );

    if (direction == L2R)
    {
        double h2 = h/2.0;
        ry[0] = y0;

//        double xn = x0;

        DoubleVector yn(m);

        /* initializing */
//        for (unsigned int j=0; j<n; j++) ry[j][0] = y0[j];

        for (int i=min+1; i<=max; i++)
        {
            const unsigned int mi = static_cast<unsigned int>(i-min);
            const PointNodeODE node1((i-1)*h, i-1);
            const PointNodeODE node2((i-1)*h+h2, i-1);

            // k1
            //for (unsigned int j=0; j<n; j++) yn[j] = ry[j][(i-1)-min];
            for (unsigned int j=0; j<m; j++) k1[j] = f(node1, ry[mi-1], j+1);

            // k2
            for (unsigned int j=0; j<m; j++) yn[j] = ry[mi-1][j]+h2*k1[j];
            for (unsigned int j=0; j<m; j++) k2[j] = f(node2, yn, j+1);

            for (unsigned int j=0; j<m; j++) ry[mi][j] = ry[mi-1][j] + h * k2[j];
        }
    }

    const unsigned int n = y0.length();

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

void FirstOrderNonLinearODE::calculateRK4(double x0, const DoubleVector &y0, std::vector<DoubleVector> &ry, Direction direction)
{
    const int min = dimension().min();
    const int max = dimension().max();
    const unsigned int N = static_cast<unsigned int>( dimension().size() );
    const double h = dimension().step();
    const unsigned int m = count();

    ry.resize(N+1); for (unsigned int n=0; n<=N; n++) ry[n].resize(N+1);

    double *k1 = static_cast<double*>( malloc( sizeof(double)*m ) );
    double *k2 = static_cast<double*>( malloc( sizeof(double)*m ) );
    double *k3 = static_cast<double*>( malloc( sizeof(double)*m ) );
    double *k4 = static_cast<double*>( malloc( sizeof(double)*m ) );

    if (direction == L2R)
    {
        double h2 = h/2.0;
        double h6 = h/6.0;
        ry[0] = y0;

        DoubleVector yn(m);

        for (int i=min+1; i<=max; i++)
        {
            const unsigned int mi = static_cast<unsigned int>(i-min);
            const PointNodeODE node1((i-1)*h, i-1);
            const PointNodeODE node2((i-1)*h+h2, i-1);
            const PointNodeODE node3((i-1)*h+h2, i-1);
            const PointNodeODE node4((i-1)*h+h, i-1);

            // k1
            for (unsigned int j=0; j<m; j++) k1[j] = f(node1, ry[mi-1], j+1);

            // k2
            for (unsigned int j=0; j<m; j++) yn[j] = ry[mi-1][j]+h2*k1[j];
            for (unsigned int j=0; j<m; j++) k2[j] = f(node2, yn, j+1);

            // k3
            for (unsigned int j=0; j<m; j++) yn[j] = ry[mi-1][j]+h2*k2[j];
            for (unsigned int j=0; j<m; j++) k3[j] = f(node3, yn, j+1);

            // k4
            for (unsigned int j=0; j<m; j++) yn[j] = ry[mi-1][j]+h*k3[j];
            for (unsigned int j=0; j<m; j++) k4[j] = f(node4, yn, j+1);

            for (unsigned int j=0; j<m; j++) ry[mi][j] = ry[mi-1][j] + h6 * (k1[j] + 2*k2[j] + 2*k3[j] + k4[j]);
        }
    }

    const unsigned int n = y0.length();
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

void FirstOrderNonLinearODE::calculateEuler(double x0, const DoubleVector &y0, std::vector<DoubleVector> &ry, Direction direction)
{
    const int min = dimension().min();
    const int max = dimension().max();
    const unsigned int N = static_cast<unsigned int>( dimension().size() );
    const double h = dimension().step();
    const unsigned int m = count();

    ry.resize(N+1); for (unsigned int n=0; n<=N; n++) ry[n].resize(m);

    if (direction == L2R)
    {
        ry[0] = y0;
        for (int i=min+1; i<=max; i++)
        {
            const unsigned int mi = static_cast<unsigned int>(i-min);
            const PointNodeODE node((i-1)*h, i-1);
            for (unsigned int j=0; j<m; j++)
                ry[mi][j] = ry[mi-1][j] + h * f(node, ry[mi-1], j+1);
        }
    }

    if (direction == R2L)
    {
        const unsigned int n = y0.length();
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

void FirstOrderNonLinearODE::calculateEulerMod(double x0, const DoubleVector &y0, std::vector<DoubleVector> &ry, Direction direction)
{
    const int min = dimension().min();
    const int max = dimension().max();
    const unsigned int N = static_cast<unsigned int>( dimension().size() );
    const double h = dimension().step();
    const unsigned int m = count();

    ry.resize(N+1); for (unsigned int n=0; n<=N; n++) ry[n].resize(m);

    if (direction == L2R)
    {
        ry[0] = y0;
        for (int i=min+1; i<=max; i++)
        {
            const unsigned int mi = static_cast<unsigned int>(i-min);
            const PointNodeODE node((i-1)*h, i-1);
            for (unsigned int j=0; j<m; j++)
                ry[mi][j] = ry[mi-1][j] + h * f(node, ry[mi-1], j+1);
        }
    }

    if (direction == R2L)
    {
        const unsigned int n = y0.length();
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

void FirstOrderNonLinearODE::cauchyProblem(double x0, const DoubleVector &y0, DoubleVector &ry, OdeSolverMethod method, Direction direction)
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

void FirstOrderNonLinearODE::calculateRK4(double x0, const DoubleVector &y0, DoubleVector &ry, Direction direction)
{
    const int min = dimension().min();
    const int max = dimension().max();
    const int N = dimension().size();
    const double h = dimension().step();
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

