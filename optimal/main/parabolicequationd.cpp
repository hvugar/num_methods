#include "parabolicequationd.h"

void ParabolicEquationD::calculateU(DoubleMatrix &u, double hx, double ht, unsigned int N, unsigned int M, double a) const
{
    u.clear();
    u.resize(M+1, N+1);

    DoubleVector da(N-1);
    DoubleVector db(N-1);
    DoubleVector dc(N-1);
    DoubleVector dd(N-1);
    DoubleVector rx(N-1);

    double alpha = -(a*a*ht)/(hx*hx);
    double beta  = 1.0 + (2.0*a*a*ht)/(hx*hx);

    /* initial condition */
    for (unsigned int n=0; n<=N; n++) u[0][n] = initial(n, hx);

    for (unsigned int m=1; m<=M; m++)
    {
        /* boundary condition */
        u[m][0] = boundary(Left, m, ht);
        u[m][N] = boundary(Right, m, ht);

        for (unsigned int n=1; n<=N-1; n++)
        {
            da[n-1] = alpha;
            db[n-1] = beta;
            dc[n-1] = alpha;
            dd[n-1] = u[m-1][n] + ht*f(n, m, hx, ht);
        }

        da[0]   = 0.0;
        dc[N-2] = 0.0;

        dd[0]   -= alpha * u[m][0];
        dd[N-2] -= alpha * u[m][N];

        tomasAlgorithm(da.data(), db.data(), dc.data(), dd.data(), rx.data(), rx.size());

        for (unsigned int n=1; n<=N-1; n++)
        {
            u[m][n] = rx[n-1];
        }
    }

    da.clear();
    db.clear();
    dc.clear();
    dd.clear();
    rx.clear();
}
