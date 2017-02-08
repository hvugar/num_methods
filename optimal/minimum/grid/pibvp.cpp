#include "pibvp.h"

void ParabolicIBVP::gridMethod(DoubleVector &u, SweepMethodDirection direction)
{

}

void ParabolicIBVP::gridMethod(DoubleMatrix &u, SweepMethodDirection direction)
{
    typedef void (*t_algorithm)(const double*, const double*, const double*, const double*, double*, unsigned int);
    t_algorithm algorithm = &tomasAlgorithm;
    if (direction == ForwardSweep) algorithm = &tomasAlgorithmL2R;
    if (direction == BackwardSweep) algorithm = &tomasAlgorithmR2L;

    Dimension time = mtimeDimension;
    Dimension dim1 = mspaceDimension.at(0);

    double ht = time.step();
    unsigned int minM = time.minN();
    unsigned int maxM = time.maxN();
    unsigned int M = maxM-minM;

    double hx = dim1.step();
    unsigned int minN = dim1.minN();
    unsigned int maxN = dim1.maxN();
    unsigned int N = maxN-minN;

    double h = ht/(hx*hx);

    u.clear();
    u.resize(M+1, N+1);

    double *ka = (double*) malloc(sizeof(double)*(N-1));
    double *kb = (double*) malloc(sizeof(double)*(N-1));
    double *kc = (double*) malloc(sizeof(double)*(N-1));
    double *kd = (double*) malloc(sizeof(double)*(N-1));
    double *rx = (double*) malloc(sizeof(double)*(N-1));

    /* initial condition */
    SpaceNode isn;
    for (unsigned int n=0; n<=N; n++)
    {
        isn.i = n+minN;
        isn.x = isn.i*hx;
        u[0][n] = initial(isn);
    }

    SpaceNode lsn;
    lsn.i = minN;
    lsn.x = minN*hx;

    SpaceNode rsn;
    rsn.i = maxN;
    rsn.x = maxN*hx;

    TimeNode tn;
    for (unsigned int m=1; m<=M; m++)
    {
        tn.i = m+minM;
        tn.t = tn.i*ht;

        for (unsigned int n=1; n<=N-1; n++)
        {
            isn.i = n+minN;
            isn.x = isn.i*hx;

            double alpha = -a(isn,tn)*h;
            double betta = 1.0 - 2.0*alpha;

            ka[n-1] = alpha;
            kb[n-1] = betta;
            kc[n-1] = alpha;
            kd[n-1] = u[m-1][n] + ht * f(isn, tn);
        }

        ka[0]   = 0.0;
        kc[N-2] = 0.0;

        /* border conditions */
        u[m][0] = boundary(lsn, tn, Left);
        u[m][N] = boundary(rsn, tn, Right);

        kd[0]   += a(lsn,tn) * h * u[m][0];
        kd[N-2] += a(rsn,tn) * h * u[m][N];

        (*algorithm)(ka, kb, kc, kd, rx, N-1);

        for (unsigned int n=1; n<=N-1; n++) u[m][n] = rx[n-1];
    }

    free(ka);
    free(kb);
    free(kc);
    free(kd);
    free(rx);
}

void ParabolicIBVP::calculateN2L2RD(DoubleMatrix &u)
{
    Dimension time = mtimeDimension;
    Dimension dim1 = mspaceDimension.at(0);

    double ht = time.step();
    unsigned int minM = time.minN();
    unsigned int maxM = time.maxN();
    unsigned int M = maxM - minM;

    double hx = dim1.step();
    unsigned int minN = dim1.minN();
    unsigned int maxN = dim1.maxN();
    unsigned int N = maxN - minN;

    const unsigned int k = 2;
    double h = ht/(hx*hx);

    u.clear();
    u.resize(M+1, N+1);

    DoubleMatrix A(k, k, 0.0);
    DoubleVector b(k, 0.0);
    DoubleVector x(k, 0.0);
    DoubleMatrix ems(N-k, k);

    /* initial condition */
    SpaceNode isn;
    for (unsigned int n=0; n<=N; n++)
    {
        isn.i = n+minN;
        isn.x = isn.i*hx;
        u[0][n] = initial(isn);
    }

    TimeNode tn;
    SpaceNode lsn;
    SpaceNode rsn;
    lsn.i = minN; lsn.x = minN*hx;
    rsn.i = maxN; rsn.x = maxN*hx;
    for (unsigned int m=1; m<=M; m++)
    {
        tn.i = m+minM;
        tn.t = tn.i*ht;

        /* border conditions */
        u[m][0] = boundary(lsn, tn, Left);
        u[m][N] = boundary(rsn, tn, Right);

        /* n=1 */
        isn.i = minN+1;
        isn.x = isn.i*hx;

        double alpha = a(isn,tn)*h;
        A[0][0] = -2.0*alpha - 1.0;
        A[0][1] = alpha;
        b[0]    = -u[m-1][1] - alpha*u[m][0] - ht*f(isn,tn);

        A[0][1] /= A[0][0];
        b[0]    /= A[0][0];
        A[0][0] = 1.0;

        ems[0][0] = A[0][1];
        ems[0][1] = b[0];

        for (unsigned int n=2; n<=N-k; n++)
        {
            isn.i = n+minN;
            isn.x = isn.i*hx;

            alpha = a(isn,tn)*h;

            double g1 = alpha;
            double g2 = -2.0*alpha-1.0;
            double g3 = alpha;
            double fi = -u[m-1][n] - ht*f(isn,tn);

            g2 /= -g1;
            g3 /= -g1;
            fi /= +g1;
            g1  = 1.0;

            A[0][0] = A[0][1] + g2;
            A[0][1] = g3;
            b[0]    = b[0] - fi;
            \
            A[0][1] /= A[0][0];
            b[0]    /= A[0][0];
            A[0][0] = 1.0;

            ems[n-1][0] = A[0][1];
            ems[n-1][1] = b[0];
        }

        isn.i = maxN-1;
        isn.x = isn.i*hx;
        alpha = a(isn,tn)*h;
        A[1][0] = alpha;
        A[1][1] = -2.0*alpha - 1.0;
        b[1]    = -u[m-1][N-1] - alpha*u[m][N] - ht*f(isn,tn);

        GaussianElimination(A, b, x);

        u[m][N-1] = x[1];
        u[m][N-2] = x[0];
        for (unsigned int n=N-(k+1); n>=1; n--)
        {
            u[m][n] = -ems[n-1][0]*u[m][n+1]+ems[n-1][1];
        }
    }

    ems.clear();
    x.clear();
    b.clear();
    A.clear();
}

void ParabolicIBVP::calculateN4L2RD(DoubleMatrix &u)
{
    Dimension time = mtimeDimension;
    Dimension dim1 = mspaceDimension.at(0);

    double ht = time.step();
    unsigned int minM = time.minN();
    unsigned int maxM = time.maxN();
    unsigned int M = maxM - minM;

    double hx = dim1.step();
    unsigned int minN = dim1.minN();
    unsigned int maxN = dim1.maxN();
    unsigned int N = maxN - minN;

    const unsigned int k = 4;
    double h = ht/(24.0*hx*hx);

    u.clear();
    u.resize(M+1, N+1);

    DoubleMatrix A(k, k, 0.0);
    DoubleVector b(k, 0.0);
    DoubleVector x(k, 0.0);
    DoubleMatrix ems(N-k, k);

    /* initial condition */
    SpaceNode isn;
    for (unsigned int n=0; n<=N; n++)
    {
        isn.i = n+minN;
        isn.x = isn.i*hx;
        u[0][n] = initial(isn);
    }

    TimeNode tn;
    SpaceNode lsn;
    SpaceNode rsn;
    lsn.i = minN; lsn.x = minN*hx;
    rsn.i = maxN; rsn.x = maxN*hx;
    for (unsigned int m=1; m<=M; m++)
    {
        tn.i = m+minM;
        tn.t = tn.i*ht;

        /* border conditions */
        u[m][0] = boundary(lsn, tn, Left);
        u[m][N] = boundary(rsn, tn, Right);

        /* n=1 */
        isn.i = minN+1;
        isn.x = isn.i*hx;

        /* n=1 */
        double alpha = a(isn,tn)*h;
        A[0][0] = -40.0*alpha - 1.0;
        A[0][1] = +12.0*alpha;
        A[0][2] = +8.0*alpha;
        A[0][3] = -2.0*alpha;
        b[0]    = -u[m-1][1] - (+22.0*alpha)*u[m][0] - ht*f(isn,tn);

        A[0][1] /= A[0][0];
        A[0][2] /= A[0][0];
        A[0][3] /= A[0][0];
        b[0]    /= A[0][0];
        A[0][0] = 1.0;

        ems[0][0] = A[0][1];
        ems[0][1] = A[0][2];
        ems[0][2] = A[0][3];
        ems[0][3] = b[0];

        for (unsigned int n=1; n<=N-(k+1); n++)
        {
            isn.i = n+minN;
            isn.x = isn.i*hx;

            alpha = a(isn,tn)*h;

            double g1 = +70.0*alpha-1.0;
            double g2 = -208.0*alpha;
            double g3 = +228.0*alpha;
            double g4 = -112.0*alpha;
            double g5 = +22.0*alpha;
            double fi = -u[m-1][n] - ht*f(isn,tn);

            g2 /= -g1;
            g3 /= -g1;
            g4 /= -g1;
            g5 /= -g1;
            fi /= +g1;
            g1 = 1.0;

            A[0][0] = A[0][1] + g2;
            A[0][1] = A[0][2] + g3;
            A[0][2] = A[0][3] + g4;
            A[0][3] = g5;
            b[0]    = b[0] - fi;
            \
            A[0][1] /= A[0][0];
            A[0][2] /= A[0][0];
            A[0][3] /= A[0][0];
            b[0]    /= A[0][0];
            A[0][0] = 1.0;

            ems[n][0] = A[0][1];
            ems[n][1] = A[0][2];
            ems[n][2] = A[0][3];
            ems[n][3] = b[0];
        }

        isn.i = maxN-3;
        isn.x = isn.i*hx;
        alpha = a(isn,tn)*h;
        A[1][0] = +22.0*alpha;
        A[1][1] = -40.0*alpha - 1.0;
        A[1][2] = +12.0*alpha;
        A[1][3] = +8.0*alpha;
        b[1]    = -u[m-1][N-3] - (-2.0*alpha)*u[m][N] - ht*f(isn,tn);

        isn.i = maxN-2;
        isn.x = isn.i*hx;
        alpha = a(isn,tn)*h;
        A[2][0] = -2.0*alpha;
        A[2][1] = +32.0*alpha;
        A[2][2] = -60.0*alpha - 1.0;
        A[2][3] = +32.0*alpha;
        b[2]    = -u[m-1][N-2] - (-2.0*alpha)*u[m][N] - ht*f(isn,tn);

        isn.i = maxN-1;
        isn.x = isn.i*hx;
        A[3][0] = -2.0*alpha;
        A[3][1] = +8.0*alpha;
        A[3][2] = +12.0*alpha;
        A[3][3] = -40.0*alpha - 1.0;
        b[3]    = -u[m-1][N-1] - (+22.0*alpha)*u[m][N] - ht*f(isn,tn);

        GaussianElimination(A, b, x);

        u[m][N-1] = x[3];
        u[m][N-2] = x[2];
        u[m][N-3] = x[1];
        u[m][N-4] = x[0];
        for (unsigned int n=N-(k+1); n>=1; n--)
        {
            u[m][n] = -ems[n-1][0]*u[m][n+1]
                      -ems[n-1][1]*u[m][n+2]
                      -ems[n-1][2]*u[m][n+3]
                      +ems[n-1][3];
        }
    }

    ems.clear();
    x.clear();
    b.clear();
    A.clear();
}

void ParabolicIBVP::calculateMVD(DoubleMatrix &u) const
{
    Dimension time = mtimeDimension;
    Dimension dim1 = mspaceDimension.at(0);
    Dimension dim2 = mspaceDimension.at(1);

    double ht = time.step();
    double h1 = dim1.step();
    double h2 = dim2.step();
    unsigned int M = time.maxN();
    unsigned int N1 = dim1.maxN();
    unsigned int N2 = dim2.maxN();

    double a1 = 1.0;
    double a2 = 1.0;

    //cleaning matrix
    u.resize(N2+1, N1+1);

    DoubleMatrix uh(N2+1, N1+1);

    DoubleVector da1(N1-1);
    DoubleVector db1(N1-1);
    DoubleVector dc1(N1-1);
    DoubleVector dd1(N1-1);
    DoubleVector rx1(N1-1);

    DoubleVector da2(N2-1);
    DoubleVector db2(N2-1);
    DoubleVector dc2(N2-1);
    DoubleVector dd2(N2-1);
    DoubleVector rx2(N2-1);

    double x1_a = -(a1*ht)/(2.0*h1*h1);
    double x1_b  = 1.0 + (a1*ht)/(h1*h1);
    double x1_c = (a2*ht)/(2.0*h2*h2);

    double x2_a = -(a2*ht)/(2.0*h2*h2);
    double x2_b  = 1.0 + (a2*ht)/(h2*h2);
    double x2_c = (a1*ht)/(2.0*h1*h1);

    for (unsigned int j=0; j<=N2; j++)
    {
        SpaceNode sn;
        sn.j = j;
        sn.y = j*h2;
        for (unsigned int i=0; i<=N1; i++)
        {
            sn.i = i;
            sn.x = i*h1;
            u[j][i] = initial(sn);
        }
    }

    TimeNode tn;
    SpaceNode sn;
    for (unsigned int k=1; k<=M; k++)
    {
        tn.i = 2*k-1;
        tn.t = 0.5*(2*k-1)*ht;
        // Approximation to x1 direction
        for (unsigned int j=1; j<N2; j++)
        {
            sn.j = j;
            sn.y = j*h2;
            for (unsigned int i=1; i<N1; i++)
            {
                sn.i = i;
                sn.x = i*h1;

                da1[i-1] = x1_a;
                db1[i-1] = x1_b;
                dc1[i-1] = x1_a;
                dd1[i-1] = x1_c*(u[j-1][i] - 2.0*u[j][i] + u[j+1][i]) + u[j][i] + (ht/2.0) * f(sn, tn);
            }

            da1[0]     = 0.0;
            dc1[N1-2]  = 0.0;

            SpaceNode sn0;
            sn0.i = 0;
            sn0.x = 0.0;
            sn0.j = j;
            sn0.y = j*h2;
            uh[j][0]  = boundary(sn0, tn);
            SpaceNode snN;
            snN.i = N1;
            snN.x = N1*h1;
            snN.j = j;
            snN.y = j*h2;
            uh[j][N1] = boundary(snN, tn);

            dd1[0]    -= x1_a * uh[j][0];
            dd1[N1-2] -= x1_a * uh[j][N1];

            tomasAlgorithm(da1.data(), db1.data(), dc1.data(), dd1.data(), rx1.data(), rx1.size());

            for (unsigned int i=1; i<N1; i++)
            {
                uh[j][i] = rx1[i-1];
            }
        }

        for (unsigned int i=0; i<=N1; i++)
        {
            SpaceNode sn0;
            sn0.i = i;
            sn0.x = i*h1;
            sn0.j = 0;
            sn0.y = 0.0;
            uh[0][i]  = boundary(sn0, tn);
            SpaceNode snN;
            snN.i = i;
            snN.x = i*h1;
            snN.j = N2;
            snN.y = N2*h2;
            uh[N2][i] = boundary(snN, tn);
        }

        tn.i = k;
        tn.t = k*ht;
        // Approximation to x2 direction
        for (unsigned int i=1; i<N1; i++)
        {
            sn.i = i;
            sn.x = i*h1;
            for (unsigned int j=1; j<N2; j++)
            {
                sn.j = j;
                sn.y = j*h2;
                da2[j-1] = x2_a;
                db2[j-1] = x2_b;
                dc2[j-1] = x2_a;
                dd2[j-1] = x2_c*(uh[j][i-1] - 2.0*uh[j][i] + uh[j][i+1]) + uh[j][i] + (ht/2.0) * f(sn, tn);
            }
            da2[0]     = 0.0;
            dc2[N2-2]  = 0.0;

            SpaceNode sn0;
            sn0.i = i;
            sn0.x = i*h1;
            sn0.j = 0;
            sn0.y = 0.0;
            u[0][i]  = boundary(sn0, tn);
            SpaceNode snN;
            snN.i = i;
            snN.x = i*h1;
            snN.j = N2;
            snN.y = N2*h2;
            u[N2][i] = boundary(snN, tn);

            dd2[0]    -= x2_a * u[0][i];
            dd2[N2-2] -= x2_a * u[N2][i];

            tomasAlgorithm(da2.data(), db2.data(), dc2.data(), dd2.data(), rx2.data(), rx2.size());

            for (unsigned int j=1; j<N2; j++)
            {
                u[j][i] = rx2[j-1];
            }
        }

        for (unsigned int j=0; j<=N2; j++)
        {
            SpaceNode sn0;
            sn0.i = 0;
            sn0.x = 0*h1;
            sn0.j = j;
            sn0.y = j*h2;
            u[j][0]  = boundary(sn0, tn);
            SpaceNode snN;
            snN.i = N1;
            snN.x = N1*h1;
            snN.j = j;
            snN.y = j*h2;
            u[j][N1] = boundary(snN, tn);
        }
    }

    da1.clear();
    db1.clear();
    dc1.clear();
    dd1.clear();
    rx1.clear();

    da2.clear();
    db2.clear();
    dc2.clear();
    dd2.clear();
    rx2.clear();
}
