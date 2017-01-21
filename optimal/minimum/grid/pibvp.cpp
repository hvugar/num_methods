#include "pibvp.h"

void ParabolicIBVP::gridMethod(DoubleMatrix &u, SweepMethodDirection direction)
{
    Grid g = grid();
    double ht = g.ht;
    unsigned int M = g.M;
    double hx = g.hx1;
    unsigned int N = g.N1;

    //double k = ht/(hx*hx);

    u.clear();
    u.resize(M+1, N+1);

    double *da = (double*) malloc(sizeof(double)*(N-1));
    double *db = (double*) malloc(sizeof(double)*(N-1));
    double *dc = (double*) malloc(sizeof(double)*(N-1));
    double *dd = (double*) malloc(sizeof(double)*(N-1));
    double *rx = (double*) malloc(sizeof(double)*(N-1));

    for (unsigned int n=0; n<=N; n++) u[0][n] = initial(n);

    for (unsigned int m=1; m<=M; m++)
    {
        u[m][0] = boundary(m, Left);
        u[m][N] = boundary(m, Right);

        for (unsigned int n=1; n<=N-1; n++)
        {
            double alpha = -a(n,m)*(ht/(hx*hx));
            double betta = 1.0 - 2.0*alpha;

            da[n-1] = alpha;
            db[n-1] = betta;
            dc[n-1] = alpha;
            dd[n-1] = u[m-1][n] + ht * f(n, m);
        }

        da[0]   = 0.0;
        dc[N-2] = 0.0;

        double alpha0 = -a(0,m)*(ht/(hx*hx));
        dd[0] -= alpha0 * u[m][0];

        double alphaN = -a(N,m)*(ht/(hx*hx));
        dd[N-2] -= alphaN * u[m][N];

        //tomasAlgorithm(da, db, dc, dd, rx, N-1);
        if (direction == ForwardSweep)
            tomasAlgorithmL2R(da, db, dc, dd, rx, N-1);
        else
            tomasAlgorithmR2L(da, db, dc, dd, rx, N-1);

        for (unsigned int n=1; n<=N-1; n++) u[m][n] = rx[n-1];
    }

    free(da);
    free(db);
    free(dc);
    free(dd);
    free(rx);
}

void ParabolicIBVP::calculateN4L2RD(DoubleMatrix &u)
{
    Grid g = grid();
    double ht = g.ht;
    unsigned int M = g.M;
    double hx = g.hx1;
    unsigned int N = g.N1;

    double a = 1.0;
    const unsigned int k = 4;
    double alpha = (ht*a*a)/(24.0*hx*hx);

    u.clear();
    u.resize(M+1, N+1);

    /*
    double D[k+1][k+1] =
    {
        {+70.0, -208.0, +228.0, -112.0, +22.0},
        {+22.0, -40.0,  +12.0,  +8.0,   -2.0},
        {-2.0,  +32.0,  -60.0,  +32.0,  -2.0},
        {-2.0,  +8.0,   +12.0,  -40.0,  +22.0},
        {+22.0, -112.0, +228.0, -208.0, +70.0}
    };
    */

    DoubleMatrix A(k, k, 0.0);
    DoubleVector b(k, 0.0);
    DoubleVector x(k, 0.0);
    DoubleMatrix ems(N-k, k);

    /* initial condition */
    for (unsigned int n=0; n<=N; n++) u.at(0,n) = initial(n);

    /* border conditions */
    for (unsigned int m=1; m<=M; m++)
    {
        u.at(m,0) = boundary(m, Left);
        u.at(m,N) = boundary(m, Right);
    }

    for (unsigned int m=1; m<=M; m++)
    {
        A[0][0] = -40.0*alpha - 1.0;
        A[0][1] = +12.0*alpha;
        A[0][2] = +8.0*alpha;
        A[0][3] = -2.0*alpha;
        b[0]    = -u.at(m-1,1) - (+22.0*alpha)*u.at(m,0) - ht*f(1,m);

        A[0][1] /= A[0][0];
        A[0][2] /= A[0][0];
        A[0][3] /= A[0][0];
        b[0]    /= A[0][0];
        A[0][0] = 1.0;

        ems.at(0,0) = A[0][1];
        ems.at(0,1) = A[0][2];
        ems.at(0,2) = A[0][3];
        ems.at(0,3) = b[0];

        for (unsigned int n=1; n<=N-(k+1); n++)
        {
            double g1 = +70.0*alpha-1.0;
            double g2 = -208.0*alpha;
            double g3 = +228.0*alpha;
            double g4 = -112.0*alpha;
            double g5 = +22.0*alpha;
            double g0 = u.at(m-1,n) + ht*f(n,m);

            g2 /= -g1;
            g3 /= -g1;
            g4 /= -g1;
            g5 /= -g1;
            g0 /= -g1;
            g1 = 1.0;

            A[0][0] = A[0][1] + g2;
            A[0][1] = A[0][2] + g3;
            A[0][2] = A[0][3] + g4;
            A[0][3] = g5;
            b[0]    = b[0] - g0;
            \
            A[0][1] /= A[0][0];
            A[0][2] /= A[0][0];
            A[0][3] /= A[0][0];
            b[0]    /= A[0][0];
            A[0][0] = 1.0;

            ems.at(n,0) = A[0][1];
            ems.at(n,1) = A[0][2];
            ems.at(n,2) = A[0][3];
            ems.at(n,3) = b[0];
        }

        A[1][0] = +22.0*alpha;
        A[1][1] = -40.0*alpha - 1.0;
        A[1][2] = +12.0*alpha;
        A[1][3] = +8.0*alpha;
        b[1]    = -u.at(m-1,N-3) - (-2.0*alpha)*u.at(m,N) - ht*f(N-3,m);

        A[2][0] = -2.0*alpha;
        A[2][1] = +32.0*alpha;
        A[2][2] = -60.0*alpha - 1.0;
        A[2][3] = +32.0*alpha;
        b[2]    = -u.at(m-1,N-2) - (-2.0*alpha)*u.at(m,N) - ht*f(N-2,m);

        A[3][0] = -2.0*alpha;
        A[3][1] = +8.0*alpha;
        A[3][2] = +12.0*alpha;
        A[3][3] = -40.0*alpha - 1.0;
        b[3]    = -u.at(m-1,N-1) - (+22.0*alpha)*u.at(m,N) - ht*f(N-1,m);

        GaussianElimination(A, b, x);

        u.at(m, N-1) = x.at(3);
        u.at(m, N-2) = x.at(2);
        u.at(m, N-3) = x.at(1);
        u.at(m, N-4) = x.at(0);
        for (unsigned int i=N-(k+1); i>=1; i--)
        {
            u.at(m,i) = -ems.at(i-1,0)*u.at(m,i+1) - ems.at(i-1,1)*u.at(m,i+2) - ems.at(i-1,2)*u.at(m,i+3) + ems.at(i-1,3);
        }
    }

    ems.clear();
    x.clear();
    b.clear();
    A.clear();
}
