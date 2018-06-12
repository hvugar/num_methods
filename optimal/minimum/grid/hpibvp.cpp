#include "hpibvp.h"

void HeatEquationIBVP::gridMethod(DoubleVector &u, double a, SweepMethodDirection direction)
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
    u.resize(N+1);

    double *ka = (double*) malloc(sizeof(double)*(N-1));
    double *kb = (double*) malloc(sizeof(double)*(N-1));
    double *kc = (double*) malloc(sizeof(double)*(N-1));
    double *kd = (double*) malloc(sizeof(double)*(N-1));
    double *rx = (double*) malloc(sizeof(double)*(N-1));

    /* initial condition */
    SpaceNodePDE isn;
    for (unsigned int n=0; n<=N; n++)
    {
        isn.i = n+minN;
        isn.x = isn.i*hx;
        u[n] = initial(isn);
    }
    layerInfo(u, 0);

    SpaceNodePDE lsn;
    lsn.i = minN;
    lsn.x = minN*hx;

    SpaceNodePDE rsn;
    rsn.i = maxN;
    rsn.x = maxN*hx;

    TimeNodePDE tn;
    for (unsigned int m=1; m<=M; m++)
    {
        tn.i = m+minM;
        tn.t = tn.i*ht;

        for (unsigned int n=1; n<=N-1; n++)
        {
            isn.i = n+minN;
            isn.x = isn.i*hx;

            double alpha = -a*a*h;
            double betta = 1.0 - 2.0*alpha;

            ka[n-1] = alpha;
            kb[n-1] = betta;
            kc[n-1] = alpha;
            kd[n-1] = u[n] + ht * f(isn, tn);
        }

        ka[0]   = 0.0;
        kc[N-2] = 0.0;

        /* border conditions */
        u[0] = boundary(lsn, tn, Left);
        u[N] = boundary(rsn, tn, Right);

        kd[0]   += a*a * h * u[0];
        kd[N-2] += a*a * h * u[N];

        (*algorithm)(ka, kb, kc, kd, rx, N-1);

        for (unsigned int n=1; n<=N-1; n++) u[n] = rx[n-1];
        layerInfo(u, m);
    }

    free(ka);
    free(kb);
    free(kc);
    free(kd);
    free(rx);
}

void HeatEquationIBVP::gridMethod1(DoubleVector &u, double a)
{
    Dimension time = mtimeDimension;
    Dimension dimX = mspaceDimension.at(0);

    double ht = time.step();
    unsigned int M = time.sizeN();

    double hx = dimX.step();
    unsigned int N =dimX.sizeN();

    u.clear();
    u.resize(N+1);

    double *ka = (double*) malloc(sizeof(double)*(N+1));
    double *kb = (double*) malloc(sizeof(double)*(N+1));
    double *kc = (double*) malloc(sizeof(double)*(N+1));
    double *kd = (double*) malloc(sizeof(double)*(N+1));
    double *rx = (double*) malloc(sizeof(double)*(N+1));

    /* initial condition */
    for (unsigned int n=0; n<=N; n++) u[n] = 0.0;
    layerInfo(u, 0);

    double q = 1.0;
    for (unsigned int m=1; m<=M; m++)
    {
        for (unsigned int n=0; n<=N; n++)
        {
            ka[n] = -(a*a*ht)/(hx*hx);
            kb[n] = 1.0 + 2.0*(a*a*ht)/(hx*hx);
            kc[n] = -(a*a*ht)/(hx*hx);

            kd[n] = u[n];

            //if (45 <= n && n <= 55)
            {
                kd[n] += ht * q * (1.0/sqrt(2.0*M_PI*hx*hx))*exp(-(n*hx-0.5)/(2.0*hx*hx));
            }
        }

        ka[0] = 0.0;
        kc[0] = -2.0*(a*a*ht)/(hx*hx);

        ka[N] = -2.0*(a*a*ht)/(hx*hx);
        kc[N] = 0.0;

        tomasAlgorithm(ka, kb, kc, kd, rx, N+1);

        for (unsigned int n=0; n<=N; n++) u[n] = rx[n];
        layerInfo(u, m);
    }

    free(ka);
    free(kb);
    free(kc);
    free(kd);
    free(rx);
}

void HeatEquationIBVP2D::calculateU(DoubleMatrix &u, double a, double alpha, double lambda)
{
    Dimension xd = spaceDimension(Dimension::DimensionX);
    Dimension yd = spaceDimension(Dimension::DimensionY);
    Dimension td = timeDimension();

    unsigned int N = xd.sizeN();
    unsigned int M = yd.sizeN();
    unsigned int L = td.sizeN();
    double hx = xd.step();
    double hy = yd.step();
    double ht = td.step();

    u.clear();
    u.resize(M+1, N+1);

    DoubleMatrix uh(M+1, N+1);

    SpaceNodePDE sn;
    //------------------------------------- initial conditions -------------------------------------//
    for (unsigned int m=0; m<=M; m++)
    {
        sn.j = m; sn.y = m*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = n; sn.x = n*hx;
            u[m][n] = initial(sn);
        }
    }
    layerInfo(u, 0);
    IPrinter::printMatrix(u);
    IPrinter::printSeperatorLine();
    //------------------------------------- initial conditions -------------------------------------//

    //double a2_ht__hx2 = ((a*a*ht)/(hx*hx));
    //double a2_ht__hy2 = ((a*a*ht)/(hy*hy));
    //double a2_lambda_ht__hy = (a*a*lambda*ht)/(hy);
    //double a2_lambda_ht__hx = (a*a*lambda*ht)/(hx);

    double *a1X = (double *) malloc(sizeof(double)*(N+1));
    double *b1X = (double *) malloc(sizeof(double)*(N+1));
    double *c1X = (double *) malloc(sizeof(double)*(N+1));
    double *d1X = (double *) malloc(sizeof(double)*(N+1));
    double *x1X = (double *) malloc(sizeof(double)*(N+1));

    double *a1Y = (double *) malloc(sizeof(double)*(M+1));
    double *b1Y = (double *) malloc(sizeof(double)*(M+1));
    double *c1Y = (double *) malloc(sizeof(double)*(M+1));
    double *d1Y = (double *) malloc(sizeof(double)*(M+1));
    double *x1Y = (double *) malloc(sizeof(double)*(M+1));

    TimeNodePDE tn;
    for (unsigned int l=1; l<=L; l++)
    {
        //------------------------------------- approximatin to x direction conditions -------------------------------------//
        tn.i = l;
        tn.t = l*ht - 0.5*ht;
        //--------------------------------------------------------------------------//
        for (unsigned int m=0; m<=M; m++)
        {
            sn.j = m; sn.y = m*hy;

            for (unsigned int n=0; n<=N; n++)
            {
                sn.i = n; sn.x = n*hx;

                d1X[n] = 2.0*u[m][n] + alpha*ht*env0(sn, tn) + ht*f(sn, tn);

                if (m==0)       d1X[n] += ((a*a*ht))*((u[0][n]   - 2.0*u[1][n]   + u[2][n])/(hy*hy));
                if (m>0 && m<M) d1X[n] += ((a*a*ht))*((u[m-1][n] - 2.0*u[m][n]   + u[m+1][n])/(hy*hy));
                if (m==M)       d1X[n] += ((a*a*ht))*((u[M-2][n] - 2.0*u[M-1][n] + u[M][n])/(hy*hy));

                if (n == 0)
                {
                    a1X[0] = 0.0;
                    b1X[0] = +2.0 + 2.0*((a*a*ht)/(hx*hx)) + alpha*ht - 2.0*(a*a*lambda*ht/hx);
                    c1X[0] = -2.0*((a*a*ht)/(hx*hx));
                    d1X[0] -= 2.0*(a*a*lambda*ht/hx)*env1(sn, tn);
                }
                else if (n == N)
                {
                    a1X[N] = -2.0*((a*a*ht)/(hx*hx));
                    b1X[N] = +2.0 + 2.0*((a*a*ht)/(hx*hx)) + alpha*ht - 2.0*(a*a*lambda*ht/hx);
                    c1X[N] = 0.0;
                    d1X[N] -= 2.0*(a*a*lambda*ht/hx)*env1(sn, tn);
                }
                else // n=1,...,N-1; m=1,...,M-1.
                {
                    a1X[n] = -(a*a*ht)/(hx*hx);
                    b1X[n] = +2.0 + 2.0*((a*a*ht)/(hx*hx)) + alpha*ht;
                    c1X[n] = -(a*a*ht)/(hx*hx);
                }
            }
            tomasAlgorithm(a1X, b1X, c1X, d1X, x1X, N+1);
            for (unsigned int n=0; n<=N; n++) uh[m][n] = x1X[n];
        }
        if (l==1)
        {
            IPrinter::printMatrix(uh);
            IPrinter::printSeperatorLine();
        }
        //------------------------------------- approximatin to x direction conditions -------------------------------------//

        //------------------------------------- approximatin to y direction conditions -------------------------------------//
        tn.i = l;
        tn.t = l*ht;
        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = n; sn.x = n*hx;

            for (unsigned int m=0; m<=M; m++)
            {
                sn.j = m; sn.y = m*hy;

                d1Y[m] = 2.0*uh[m][n] + alpha*ht*env0(sn, tn) + ht*f(sn, tn);

                if (n==0)       d1Y[m] += ((a*a*ht))*((uh[m][0]   - 2.0*uh[m][1]   + uh[m][2])/(hx*hx));
                if (n>0 && n<N) d1Y[m] += ((a*a*ht))*((uh[m][n-1] - 2.0*uh[m][n]   + uh[m][n+1])/(hx*hx));
                if (n==N)       d1Y[m] += ((a*a*ht))*((uh[m][N-2] - 2.0*uh[m][N-1] + uh[m][N])/(hx*hx));

                if (m == 0)
                {
                    a1Y[0] = 0.0;
                    b1Y[0] = +2.0 + 2.0*((a*a*ht)/(hy*hy)) + alpha*ht - 2.0*(a*a*lambda*ht/hy);
                    c1Y[0] = -2.0*((a*a*ht)/(hy*hy));
                    d1Y[0] -= 2.0*(a*a*lambda*ht/hy)*env1(sn, tn);
                }
                else if (m == M)
                {
                    a1Y[M] = -2.0*((a*a*ht)/(hy*hy));
                    b1Y[M] = +2.0 + 2.0*((a*a*ht)/(hy*hy)) + alpha*ht - 2.0*(a*a*lambda*ht/hy);
                    c1Y[M] = 0.0;
                    d1Y[M] -= 2.0*(a*a*lambda*ht/hy)*env1(sn, tn);
                }
                else // n=1,...,N-1; m=1,...,M-1.
                {
                    a1Y[m] = -((a*a*ht)/(hy*hy));
                    b1Y[m] = +2.0 + 2.0*((a*a*ht)/(hy*hy)) + alpha*ht;
                    c1Y[m] = -((a*a*ht)/(hy*hy));
                }
            }
            tomasAlgorithm(a1Y, b1Y, c1Y, d1Y, x1Y, M+1);
            for (unsigned int m=0; m<=M; m++) u[m][n] = x1Y[m];
        }
        //------------------------------------- approximatin to y direction conditions -------------------------------------//

        layerInfo(u, l);
        if (l==1)
        {
            IPrinter::printMatrix(u);
            IPrinter::printSeperatorLine();
        }
    }
    uh.clear();

    free(x1X);
    free(d1X);
    free(c1X);
    free(b1X);
    free(a1X);

    free(x1Y);
    free(d1Y);
    free(c1Y);
    free(b1Y);
    free(a1Y);
}
