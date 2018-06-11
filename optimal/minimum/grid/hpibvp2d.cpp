#include "hpibvp2d.h"

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
    //------------------------------------- initial conditions -------------------------------------//

    double a2_ht__hx2 = ((a*a*ht)/(hx*hx));
    double a2_ht__hy2 = ((a*a*ht)/(hy*hy));
    double a2_lambda_ht__hy = (a*a*lambda*ht)/(hy);
    double a2_lambda_ht__hx = (a*a*lambda*ht)/(hx);

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

                if (m==0)       d1X[n] += a2_ht__hy2*(u[0][n]   - 2.0*u[1][n]   + u[2][n]);
                if (m>0 && m<M) d1X[n] += a2_ht__hy2*(u[m-1][n] - 2.0*u[m][n]   + u[m+1][n]);
                if (m==M)       d1X[n] += a2_ht__hy2*(u[M-2][n] - 2.0*u[M-1][m] + u[M][n]);

                if (n == 0)
                {
                    a1X[0] = 0.0;
                    b1X[0] = +2.0 + 2.0*a2_ht__hx2 + alpha*ht - 2.0*a2_lambda_ht__hx - 2.0*lambda*a*a*ht/hx;
                    c1X[0] = -2.0*a2_ht__hx2;
                    d1X[0] += 2.0*a2_lambda_ht__hx*env1(sn, tn);
                }
                else if (n == N)
                {
                    a1X[N] = -2.0*a2_ht__hy2;
                    b1X[N] = +2.0 + 2.0*a2_ht__hx2 + alpha*ht - 2.0*a2_lambda_ht__hy;
                    c1X[N] = 0.0;
                    d1X[N] += 2.0*a2_lambda_ht__hx*env1(sn, tn);
                }
                else
                {
                    a1Y[m] = -a2_ht__hx2;
                    b1Y[m] = 2.0 + 2.0*a2_ht__hx2 + alpha*ht;
                    c1Y[m] = -a2_ht__hx2;
                }
            }
            tomasAlgorithm(a1X, b1X, c1X, d1X, x1X, N+1);
            for (unsigned int n=0; n<=N; n++) uh[m][n] = x1X[n];
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

                if (n==0)       d1Y[m] += a2_ht__hx2*(uh[m][0]   - 2.0*uh[m][1]   + uh[m][2]);
                if (n>0 && n<N) d1X[n] += a2_ht__hx2*(uh[m][n-1] - 2.0*uh[m][n]   + uh[m][n+1]);
                if (n==N)       d1X[n] += a2_ht__hx2*(uh[m][N-2] - 2.0*uh[m][N-1] + uh[m][N]);

                if (m == 0)
                {
                    a1Y[0] = 0.0;
                    b1Y[0] = +2.0 + 2.0*a2_ht__hy2 + alpha*ht - 2.0*a2_lambda_ht__hy;
                    c1Y[0] = -2.0*a2_ht__hy2;
                    d1Y[0] += 2.0*a2_lambda_ht__hy*env1(sn, tn);
                }
                else if (m == M)
                {
                    a1Y[N] = -2.0*a2_ht__hy2;
                    b1Y[N] = +2.0 + 2.0*a2_ht__hy2 + alpha*ht - 2.0*a2_lambda_ht__hy;
                    c1Y[N] = 0.0;
                    d1Y[N] += 2.0*a2_lambda_ht__hy*env1(sn, tn);
                }
                else
                {
                    a1Y[m] = -a2_ht__hy2;
                    b1Y[m] = +2.0 + 2.0*a2_ht__hy2 + alpha*ht;
                    c1Y[m] = -a2_ht__hy2;
                }
            }
            tomasAlgorithm(a1Y, b1Y, c1Y, d1Y, x1Y, M+1);
            for (unsigned int m=0; m<=M; m++) u[m][n] = x1Y[m];
        }
        //------------------------------------- approximatin to y direction conditions -------------------------------------//

        layerInfo(u, l);
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

