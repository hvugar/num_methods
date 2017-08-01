#include "cauchyp.h"

CauchyProblem1stOrder::CauchyProblem1stOrder(const Dimension &grid) : mgrid(grid)
{}

void CauchyProblem1stOrder::calculate(double x0, double y0, DoubleVector &y, CauchyProblem1stOrder::Method method, Direction direction)
{
    switch (method)
    {
    case RK2:
        calculateRK2(x0, y0, y, direction);
        break;
    case RK4:
        calculateRK4(x0, y0, y, direction);
        break;
    case EULER:
        calculateEuler(x0, y0, y, direction);
        break;
    case EULER_MOD:
        calculateEulerMod(x0, y0, y, direction);
        break;
    default:
        break;
    }
}

void CauchyProblem1stOrder::calculateRK2(double x0, double y0, DoubleVector &y, Direction direction)
{
    unsigned int N = mgrid.sizeN();
    double h = mgrid.step();

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

void CauchyProblem1stOrder::calculateRK4(double x0, double y0, DoubleVector &y, Direction direction)
{
    unsigned int N = mgrid.sizeN();
    double h = mgrid.step();

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

void CauchyProblem1stOrder::calculateEuler(double x0, double y0, DoubleVector &y, Direction direction)
{
    unsigned int N = mgrid.sizeN();
    double h = mgrid.step();

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

void CauchyProblem1stOrder::calculateEulerMod(double x0, double y0, DoubleVector &y, Direction direction)
{
    unsigned int N = mgrid.sizeN();
    double h = mgrid.step();

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

CauchyProblemM1stOrder::CauchyProblemM1stOrder(const Dimension &grid) : mgrid(grid)
{}

void CauchyProblemM1stOrder::calculateCP(double x0, const DoubleVector &y0, DoubleMatrix &y, Method method, Direction direction)
{
    switch (method)
    {
    case RK2:
        calculateRK2(x0, y0, y, direction);
        break;
    case RK4:
        calculateRK4(x0, y0, y, direction);
        break;
    case EULER:
        calculateEuler(x0, y0, y, direction);
        break;
    case EULER_MOD:
        calculateEulerMod(x0, y0, y, direction);
        break;
    default:
        break;
    }
}

void CauchyProblemM1stOrder::calculateCP(double x0, const DoubleVector &y0, DoubleVector &ry, Method method, Direction direction)
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
    default:
        break;
    }
}


void CauchyProblemM1stOrder::calculateRK2(double x0, const DoubleVector &y0, DoubleMatrix &ry, Direction direction)
{
    unsigned int minN = mgrid.minN();
    unsigned int maxN = mgrid.maxN();
    unsigned int N = mgrid.sizeN();
    double h = mgrid.step();
    unsigned int n = y0.size();

    ry.clear();
    ry.resize(n, N+1);

    double *k1 = (double *)malloc(sizeof(double)*n);
    double *k2 = (double *)malloc(sizeof(double)*n);

    if (direction == L2R)
    {
        double h2 = h/2.0;

        double xn = x0;

        DoubleVector yn(n);

        /* initializing */
        for (unsigned int j=0; j<n; j++) ry[j][0] = y0[j];

        for (unsigned int i=minN+1; i<=maxN; i++)
        {
            // k1
            for (unsigned int j=0; j<n; j++) yn[j] = ry[j][(i-1)-minN];
            for (unsigned int j=0; j<n; j++) k1[j] = f(xn, yn, i-1, j);

            // k2
            for (unsigned int j=0; j<n; j++) yn[j] = ry[j][(i-1)-minN]+h2*k1[j];
            for (unsigned int j=0; j<n; j++) k2[j] = f(xn+h2, yn, i-1, j);

            for (unsigned int j=0; j<n; j++) ry[j][i-minN] = ry[j][(i-1)-minN] + h * k2[j];

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

void CauchyProblemM1stOrder::calculateRK4(double x0, const DoubleVector &y0, DoubleMatrix &ry, Direction direction)
{
    unsigned int minN = mgrid.minN();
    unsigned int maxN = mgrid.maxN();
    unsigned int N = mgrid.sizeN();
    double h = mgrid.step();
    unsigned int n = y0.size();

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
        for (unsigned int j=0; j<n; j++) ry[j][0] = y0[j];

        for (unsigned int i=minN+1; i<=maxN; i++)
        {
            // k1
            for (unsigned int j=0; j<n; j++) yn[j] = ry[j][(i-1)-minN];
            for (unsigned int j=0; j<n; j++) k1[j] = f(xn, yn, i-1, j);

            // k2
            for (unsigned int j=0; j<n; j++) yn[j] = ry[j][(i-1)-minN]+h2*k1[j];
            for (unsigned int j=0; j<n; j++) k2[j] = f(xn+h2, yn, i-1, j);

            // k3
            for (unsigned int j=0; j<n; j++) yn[j] = ry[j][(i-1)-minN]+h2*k2[j];
            for (unsigned int j=0; j<n; j++) k3[j] = f(xn+h2, yn, i-1, j);

            // k4
            for (unsigned int j=0; j<n; j++) yn[j] = ry[j][(i-1)-minN]+h*k3[j];
            for (unsigned int j=0; j<n; j++) k4[j] = f(xn+h, yn, i-1, j);

            for (unsigned int j=0; j<n; j++) ry[j][i-minN] = ry[j][(i-1)-minN] + h6 * (k1[j] + 2*k2[j] + 2*k3[j] + k4[j]);

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

        for (unsigned int i=N-1; i!=UINT32_MAX; i--)
        {
            // k1
            for (unsigned int j=0; j<n; j++) yn[j] = ry[j][i+1];
            for (unsigned int j=0; j<n; j++) k1[j] = f(xn, yn, i+1, j);

            // k2
            for (unsigned int j=0; j<n; j++) yn[j] = ry[j][i+1]-h2*k1[j];
            for (unsigned int j=0; j<n; j++) k2[j] = f(xn-h2, yn, i+1, j);

            // k3
            for (unsigned int j=0; j<n; j++) yn[j] = ry[j][i+1]-h2*k2[j];
            for (unsigned int j=0; j<n; j++) k3[j] = f(xn-h2, yn, i+1, j);

            // k4
            for (unsigned int j=0; j<n; j++) yn[j] = ry[j][i+1]-h*k3[j];
            for (unsigned int j=0; j<n; j++) k4[j] = f(xn-h, yn, i+1, j);

            for (unsigned int j=0; j<n; j++) ry[j][i] = ry[j][i+1] - h6 * (k1[j] + 2*k2[j] + 2*k3[j] + k4[j]);

            xn -= h;
        }
    }


    free(k4);
    free(k3);
    free(k2);
    free(k1);
}

void CauchyProblemM1stOrder::calculateEuler(double x0, const DoubleVector &y0, DoubleMatrix &ry, Direction direction)
{
    unsigned int N = mgrid.sizeN();
    double h = mgrid.step();
    unsigned int n = y0.size();

    ry.clear();
    ry.resize(n, N+1);

    if (direction == L2R)
    {
        double xn = x0;
        DoubleVector yn(n);
        /* initializing */
        for (unsigned int j=0; j<n; j++) ry[j][0] = y0[j];

        for (unsigned int i=1; i<=N; i++)
        {
            for (unsigned int j=0; j<n; j++) yn[j] = ry[j][i-1];
            for (unsigned int j=0; j<n; j++) ry[j][i] = ry[j][i-1] + h * f(xn, yn, i-1, j);

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

void CauchyProblemM1stOrder::calculateEulerMod(double x0, const DoubleVector &y0, DoubleMatrix &ry, Direction direction)
{
    unsigned int N = mgrid.sizeN();
    double h = mgrid.step();
    unsigned int n = y0.size();

    ry.clear();
    ry.resize(n, N+1);

    if (direction == L2R)
    {
        double xn = x0;
        DoubleVector yn(n);
        /* initializing */
        for (unsigned int j=0; j<n; j++) ry[j][0] = y0[j];

        for (unsigned int i=1; i<=N; i++)
        {
            for (unsigned int j=0; j<n; j++) yn[j] = ry[j][i-1];
            for (unsigned int j=0; j<n; j++) ry[j][i] = ry[j][i-1] + h * f(xn, yn, i-1, j);

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

void CauchyProblemM1stOrder::calculateRK4(double x0, const DoubleVector &y0, DoubleVector &ry, Direction direction)
{
    unsigned int minN = mgrid.minN();
    unsigned int maxN = mgrid.maxN();
    unsigned int N = mgrid.sizeN();
    double h = mgrid.step();
    unsigned int n = y0.size();

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

        for (unsigned int i=minN+1; i<=maxN; i++)
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

const Dimension& CauchyProblemM1stOrder::grid() const
{
    return mgrid;
}
