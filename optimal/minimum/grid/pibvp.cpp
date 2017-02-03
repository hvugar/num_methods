#include "pibvp.h"

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
        u[m][0] = boundary(lsn, tn);
        u[m][N] = boundary(rsn, tn);

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
    for (unsigned int n=0; n<=N; n++) u[0][n] = initial(n+minN);

    for (unsigned int m=1; m<=M; m++)
    {
        unsigned int m1 = m+minM;

        /* border conditions */
        u[m][0] = boundary(m1, Left);
        u[m][N] = boundary(m1, Right);

        /* n=1 */
        double alpha = a(minN+1,m1)*h;
        A[0][0] = -2.0*alpha - 1.0;
        A[0][1] = alpha;
        b[0]    = -u[m-1][1] - alpha*u[m][0] - ht*f(minN+1,m1);

        A[0][1] /= A[0][0];
        b[0]    /= A[0][0];
        A[0][0] = 1.0;

        ems[0][0] = A[0][1];
        ems[0][1] = b[0];

        for (unsigned int n=2; n<=N-k; n++)
        {
            alpha = a(n+minN,m1)*h;

            double g1 = alpha;
            double g2 = -2.0*alpha-1.0;
            double g3 = alpha;
            double fi = -u.at(m-1,n) - ht*f(n+minN,m1);

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

        alpha = a(maxN-1,m1)*h;
        A[1][0] = alpha;
        A[1][1] = -2.0*alpha - 1.0;
        b[1]    = -u[m-1][N-1] - alpha*u[m][N] - ht*f(maxN-1,m1);

        GaussianElimination(A, b, x);

        u[m][N-1] = x.at(1);
        u[m][N-2] = x.at(0);
        for (unsigned int n=N-(k+1); n>=1; n--)
        {
            u[m][n] = -ems[n-1][0]*u[m][n+1]+ems[n-1][1];
        }
        //        for (unsigned int n=N-k; n>=2; n--)
        //        {
        //            double d0 = alpha;
        //            double d1 = -2.0*alpha-1.0;
        //            double d2 = alpha;
        //            double fi = -u.at(m-1,n) - ht*f(n,m);
        //            u.at(m,n-1) = -d1*u.at(m,n) - d2*u.at(m,n+1) + fi;
        //            u.at(m,n-1) /= d0;
        //        }
    }

    ems.clear();
    x.clear();
    b.clear();
    A.clear();
}

void ParabolicIBVP::calculateN2R2LD(DoubleMatrix &u)
{
    C_UNUSED(u);
}

void ParabolicIBVP::calculateN4L2RD(DoubleMatrix &u)
{
    Dimension time = mtimeDimension;
    Dimension dim1 = mspaceDimension.at(0);

    double ht = time.step();
    //unsigned int M1 = time.minN();
    unsigned int M = time.maxN();

    double hx = dim1.step();
    //unsigned int N1 = dim1.minN();
    unsigned int N = dim1.maxN();

    double a = 1.0;
    const unsigned int k = 4;
    double alpha = (ht*a*a)/(24.0*hx*hx);

    u.clear();
    u.resize(M+1, N+1);

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

void ParabolicIBVP::calculateN4R2LD(DoubleMatrix &u)
{
    C_UNUSED(u);
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

    double x1_a = -(a1*a1*ht)/(2.0*h1*h1);
    double x1_b  = 1.0 + (a1*a1*ht)/(h1*h1);
    double x1_c = (a2*a2*ht)/(2.0*h2*h2);

    double x2_a = -(a2*a2*ht)/(2.0*h2*h2);
    double x2_b  = 1.0 + (a2*a2*ht)/(h2*h2);
    double x2_c = (a1*a1*ht)/(2.0*h1*h1);

    for (unsigned int j=0; j<=N2; j++)
    {
        for (unsigned int i=0; i<=N1; i++)
        {
            u[j][i] = initial(i, j);
        }
    }

    for (unsigned int k=1; k<=M; k++)
    {
        // Approximation to x1 direction
        unsigned int k1 = 2*k-1;
        for (unsigned int j=1; j<N2; j++)
        {
            for (unsigned int i=1; i<N1; i++)
            {
                da1[i-1] = x1_a;
                db1[i-1] = x1_b;
                dc1[i-1] = x1_a;
                dd1[i-1] = x1_c*(u[j-1][i] - 2.0*u[j][i] + u[j+1][i]) + u[j][i] + (ht/2.0) * f(i, j, k1);
            }

            da1[0]     = 0.0;
            dc1[N1-2]  = 0.0;

            uh[j][0]  = boundary(0, j, k1);
            uh[j][N1] = boundary(N1, j, k1);

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
            uh[0][i]  = boundary(i, 0, k1);
            uh[N2][i] = boundary(i, N2, k1);
        }

        // Approximation to x2 direction
        unsigned int k2 = 2*k;
        for (unsigned int i=1; i<N1; i++)
        {
            for (unsigned int j=1; j<N2; j++)
            {
                da2[j-1] = x2_a;
                db2[j-1] = x2_b;
                dc2[j-1] = x2_a;
                dd2[j-1] = x2_c*(uh[j][i-1] - 2.0*uh[j][i] + uh[j][i+1]) + uh[j][i] + (ht/2.0) * f(i, j, k2);
            }
            da2[0]     = 0.0;
            dc2[N2-2]  = 0.0;

            u[0][i]  = boundary(i, 0, k2);
            u[N2][i] = boundary(i, N2, k2);

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
            u[j][0]  = boundary(0, j, k2);
            u[j][N1] = boundary(N1, j, k2);
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

void ParabolicIBVP::calculateMVD2(DoubleMatrix &u) const
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

    double *a01 = (double*)malloc(sizeof(double)*(N1-1));
    double *b01 = (double*)malloc(sizeof(double)*(N1-1));
    double *c01 = (double*)malloc(sizeof(double)*(N1-1));
    double *d01 = (double*)malloc(sizeof(double)*(N1-1));
    double *rx1 = (double*)malloc(sizeof(double)*(N1-1));

    double *a02 = (double*)malloc(sizeof(double)*(N2-1));
    double *b02 = (double*)malloc(sizeof(double)*(N2-1));
    double *c02 = (double*)malloc(sizeof(double)*(N2-1));
    double *d02 = (double*)malloc(sizeof(double)*(N2-1));
    double *rx2 = (double*)malloc(sizeof(double)*(N2-1));

    double a11 = -(a1*a1*ht)/(2.0*h1*h1);
    double b11  = 1.0 + (a1*a1*ht)/(h1*h1);
    double c11 = (a2*a2*ht)/(2.0*h2*h2);

    double a21 = -(a2*a2*ht)/(2.0*h2*h2);
    double b21  = 1.0 + (a2*a2*ht)/(h2*h2);
    double c21 = (a1*a1*ht)/(2.0*h1*h1);

    for (unsigned int j=0; j<=N2; j++)
    {
        for (unsigned int i=0; i<=N1; i++)
        {
            u[j][i] = initial(i, j);
        }
    }

    for (unsigned int k=1; k<=M; k++)
    {
        // Approximation to x1 direction
        unsigned int k1 = 2*k-1;
        for (unsigned int j=1; j<N2; j++)
        {
            for (unsigned int i=1; i<N1; i++)
            {
                a01[i-1] = a11;
                b01[i-1] = b11;
                c01[i-1] = a11;
                d01[i-1] = c11*(u[j-1][i] - 2.0*u[j][i] + u[j+1][i]) + u[j][i] + (ht/2.0) * f(i, j, k1);
            }

            a01[0]     = 0.0;
            c01[N1-2]  = 0.0;

            uh[j][0]  = boundary(0, j, k1);
            uh[j][N1] = boundary(N1, j, k1);

            d01[0]    -= a11 * uh[j][0];
            d01[N1-2] -= a11 * uh[j][N1];

            tomasAlgorithm(a01, b01, c01, d01, rx1, N1-1);

            for (unsigned int i=1; i<N1; i++)
            {
                uh[j][i] = rx1[i-1];
            }
        }

        for (unsigned int i=0; i<=N1; i++)
        {
            uh[0][i]  = boundary(i, 0, k1);
            uh[N2][i] = boundary(i, N2, k1);
        }

        // Approximation to x2 direction
        unsigned int k2 = 2*k;
        for (unsigned int i=1; i<N1; i++)
        {
            for (unsigned int j=1; j<N2; j++)
            {
                a02[j-1] = a21;
                b02[j-1] = b21;
                c02[j-1] = a21;
                d02[j-1] = c21*(uh[j][i-1] - 2.0*uh[j][i] + uh[j][i+1]) + uh[j][i] + (ht/2.0) * f(i, j, k2);
            }
            a02[0]     = 0.0;
            c02[N2-2]  = 0.0;

            u[0][i]  = boundary(i, 0, k2);
            u[N2][i] = boundary(i, N2, k2);

            d02[0]    -= a21 * u[0][i];
            d02[N2-2] -= a21 * u[N2][i];

            tomasAlgorithm(a02, b02, c02, d02, rx2, N2-1);

            for (unsigned int j=1; j<N2; j++)
            {
                u[j][i] = rx2[j-1];
            }
        }

        for (unsigned int j=0; j<=N2; j++)
        {
            u[j][0]  = boundary(0, j, k2);
            u[j][N1] = boundary(N1, j, k2);
        }
    }

    free(rx1);
    free(d01);
    free(c01);
    free(b01);
    free(a01);
}

void ParabolicIBVP::calculateMVD1(DoubleMatrix &u) const
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

void ParabolicIBVP::calculateN2L2RD1(DoubleMatrix &u)
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
    lsn.i = 0; lsn.x = 0.0;
    rsn.i = N; rsn.x = N*hx;
    for (unsigned int m=1; m<=M; m++)
    {
        tn.i = m+minM;
        tn.t = tn.i*ht;

        /* border conditions */
        u[m][0] = boundary(lsn, tn);
        u[m][N] = boundary(rsn, tn);

        /* n=1 */
        isn.i = minN+1;
        isn.x = isn.i*hx;

        double alpha = a(minN+1,tn)*h;
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

        A[1][0] = alpha;
        A[1][1] = -2.0*alpha - 1.0;
        isn.i = N-1;
        isn.x = (N-1)*hx;
        b[1]    = -u[m-1][N-1] - alpha*u[m][N] - ht*f(isn,tn);

        GaussianElimination(A, b, x);

        u[m][N-1] = x[1];
        u[m][N-2] = x[0];
        for (unsigned int n=N-(k+1); n>=1; n--)
        {
            u[m][n] = -ems[n-1][0]*u[m][n+1]+ems[n-1][1];
        }
        //        for (unsigned int n=N-k; n>=2; n--)
        //        {
        //            double d0 = alpha;
        //            double d1 = -2.0*alpha-1.0;
        //            double d2 = alpha;
        //            double fi = -u.at(m-1,n) - ht*f(n,m);
        //            u.at(m,n-1) = -d1*u.at(m,n) - d2*u.at(m,n+1) + fi;
        //            u.at(m,n-1) /= d0;
        //        }
    }

    ems.clear();
    x.clear();
    b.clear();
    A.clear();
}

void ParabolicIBVP::calculateN4L2RD1(DoubleMatrix &u)
{
    Dimension time = mtimeDimension;
    Dimension dim1 = mspaceDimension.at(0);

    double ht = time.step();
    //unsigned int M1 = time.minN();
    unsigned int M = time.maxN();

    double hx = dim1.step();
    //unsigned int N1 = dim1.minN();
    unsigned int N = dim1.maxN();

    double a = 1.0;
    const unsigned int k = 4;
    double alpha = (ht*a*a)/(24.0*hx*hx);

    u.clear();
    u.resize(M+1, N+1);

    DoubleMatrix A(k, k, 0.0);
    DoubleVector b(k, 0.0);
    DoubleVector x(k, 0.0);
    DoubleMatrix ems(N-k, k);

    /* initial condition */
    SpaceNode sn;
    for (unsigned int n=0; n<=N; n++)
    {
        sn.i = n;
        sn.x = n*hx;
        u.at(0,n) = initial(sn);
    }

    /* border conditions */
    TimeNode tn;
    for (unsigned int m=1; m<=M; m++)
    {
        tn.i = m;
        tn.t = m*ht;

        SpaceNode left;
        left.i = 0; left.x = 0.0;
        u[m][0] = boundary(left, tn);

        SpaceNode right;
        right.i = N; right.x = N*hx;
        u[m][N] = boundary(right, tn);
    }

    for (unsigned int m=1; m<=M; m++)
    {
        tn.i = m;
        tn.t = m*ht;

        A[0][0] = -40.0*alpha - 1.0;
        A[0][1] = +12.0*alpha;
        A[0][2] = +8.0*alpha;
        A[0][3] = -2.0*alpha;
        sn.i = 1;
        sn.x = hx;
        b[0]    = -u.at(m-1,1) - (+22.0*alpha)*u.at(m,0) - ht*f(sn,tn);

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
            sn.i = n;
            sn.x = n*hx;

            double g1 = +70.0*alpha-1.0;
            double g2 = -208.0*alpha;
            double g3 = +228.0*alpha;
            double g4 = -112.0*alpha;
            double g5 = +22.0*alpha;
            double g0 = u.at(m-1,n) + ht*f(sn,tn);

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
        sn.i = N-3;
        sn.x = (N-3)*hx;
        b[1]    = -u.at(m-1,N-3) - (-2.0*alpha)*u.at(m,N) - ht*f(sn,tn);

        A[2][0] = -2.0*alpha;
        A[2][1] = +32.0*alpha;
        A[2][2] = -60.0*alpha - 1.0;
        A[2][3] = +32.0*alpha;
        sn.i = N-2;
        sn.x = (N-2)*hx;
        b[2]    = -u.at(m-1,N-2) - (-2.0*alpha)*u.at(m,N) - ht*f(sn,tn);

        A[3][0] = -2.0*alpha;
        A[3][1] = +8.0*alpha;
        A[3][2] = +12.0*alpha;
        A[3][3] = -40.0*alpha - 1.0;
        sn.i = N-1;
        sn.x = (N-1)*hx;
        b[3]    = -u.at(m-1,N-1) - (+22.0*alpha)*u.at(m,N) - ht*f(sn,tn);

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

void ParabolicIBVP::gridMethod2(DoubleMatrix &u, SweepMethodDirection direction)
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
    for (unsigned int n=0; n<=N; n++) u[0][n] = initial(n+minN);

    for (unsigned int m=1; m<=M; m++)
    {
        unsigned int m1 = m+minM;
        for (unsigned int n=1; n<=N-1; n++)
        {
            double alpha = -a(n+minN,m1)*h;
            double betta = 1.0 - 2.0*alpha;

            ka[n-1] = alpha;
            kb[n-1] = betta;
            kc[n-1] = alpha;
            kd[n-1] = u[m-1][n] + ht * f(n+minN, m1);
        }

        ka[0]   = 0.0;
        kc[N-2] = 0.0;

        u[m][0] = boundary(m1, Left);
        u[m][N] = boundary(m1, Right);
        //double alpha0 = -a(0,m)*h;
        //kd[0] -= alpha0 * u[m][0];

        //double alphaN = -a(N,m)*h;
        //kd[N-2] -= alphaN * u[m][N];

        kd[0]   += a(0,m1)   *h*u[m][0];
        kd[N-2] += a(maxN,m1)*h*u[m][N];

        //tomasAlgorithm(da, db, dc, dd, rx, N-1);
        //if (direction == ForwardSweep)
        //    tomasAlgorithmL2R(ka, kb, kc, kd, rx, N-1);
        //else
        //    tomasAlgorithmR2L(ka, kb, kc, kd, rx, N-1);
        (*algorithm)(ka, kb, kc, kd, rx, N-1);

        for (unsigned int n=1; n<=N-1; n++) u[m][n] = rx[n-1];
    }

    free(ka);
    free(kb);
    free(kc);
    free(kd);
    free(rx);
}
