#include "problem2hnm.h"

void Problem2HNM::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    Problem2HNM prob;
    prob.a = 1.0;
    prob.lambda = 0.01;
    DoubleMatrix u;
    prob.solveEquation2D2(Dimension(0.001, 0, 1000), Dimension(0.01, 0, 100), Dimension(0.01, 0, 100), u);

    //DoubleVector u;
    //prob.solveEquation1D1(Dimension(0.001, 0, 1000), Dimension(0.01, 0, 100), u);
}

double U(const SpaceNodePDE &sn, const TimeNodePDE &tn)
{
    /********* 1D *********/
    //return sn.x + tn.t;
    //return sn.x*sn.x + tn.t;
    //return sn.x*sn.x*sn.x + tn.t;
    //return sn.x*sn.x*sn.x*sn.x + tn.t;

    //return sn.x*sn.x*sn.x + tn.t*tn.t;//works
    //return sn.x*sn.x*sn.x + tn.t*tn.t*tn.t;

    /********* 2D *********/
    //return sn.x*sn.x + sn.y*sn.y + tn.t;                    // solveEquation2D1 solveEquation2D2 solveEquation2D3 works
    //return sn.x*sn.x*sn.x + sn.y*sn.y*sn.y + tn.t;          // solveEquation2D2 solveEquation2D3 works
    //return sn.x*sn.x + sn.y*sn.y + tn.t*tn.t;               // solveEquation2D1 solveEquation2D2 solveEquation2D3 works
    return sn.x*sn.x*sn.x + sn.y*sn.y*sn.y + tn.t*tn.t;     // solveEquation2D2 solveEquation2D3 works
    //return sn.x*sn.x*sn.x + sn.y*sn.y*sn.y + tn.t*tn.t*tn.t;//
}

double Problem2HNM::initial1(const SpaceNodePDE &sn) const
{
    TimeNodePDE tn; tn.i = 0; tn.t = 0.0; return U(sn, tn);
}

double Problem2HNM::initial2(const SpaceNodePDE &) const
{
    //return 1.0;
    //return 1.0;
    //return 1.0;
    //return 1.0;

    //return 0.0;
    //return 0.0;

    //return 1.0;
    //return 1.0;
    //return 0.0;
    return 0.0;
    //return 0.0;
}

double Problem2HNM::boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    return U(sn, tn);
}

double Problem2HNM::f(const SpaceNodePDE &sn UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM) const
{
    //return 0.0 - 0.0 + lambda;
    //return 0.0 - 2.0*a*a + lambda;
    //return 0.0 - 6.0*a*a*sn.x + lambda;
    //return 0.0 - 24.0*a*a*sn.x*sn.x + lambda;

    //return 2.0 - 6.0*a*a*sn.x + lambda*2.0*tn.t;
    //return 6.0*tn.t - 6.0*a*a*sn.x + lambda*3.0*tn.t*tn.t;

    //return 0.0 - 4.0*a*a + lambda;
    //return 0.0 - (6.0*a*a*sn.x + 6.0*a*a*sn.y) + lambda;
    //return 2.0 - 4.0*a*a + lambda*2.0*tn.t;
    return 2.0 - (6.0*a*a*sn.x + 6.0*a*a*sn.y) + lambda*2.0*tn.t;
    //return 6.0*tn.t - (6.0*a*a*sn.x + 6.0*a*a*sn.y) + lambda*3.0*tn.t*tn.t;
}

void Problem2HNM::layerInfo(DoubleMatrix &u)
{
    IPrinter::printMatrix(u);
}

void Problem2HNM::layerInfo(DoubleVector &u)
{
    IPrinter::printVector(u);
}

void Problem2HNM::solveEquation2D1(const Dimension &time, const Dimension &dimx, const Dimension &dimy, DoubleMatrix &u UNUSED_PARAM)
{
    const unsigned int N = dimx.sizeN();
    const unsigned int M = dimy.sizeN();
    const unsigned int L = time.sizeN();

    const double hx = dimx.step();
    const double hy = dimy.step();
    const double ht = time.step();

    DoubleMatrix u00(M+1, N+1);
    DoubleMatrix u05(M+1, N+1);
    DoubleMatrix u10(M+1, N+1);
    DoubleMatrix u15(M+1, N+1);
    DoubleMatrix u20(M+1, N+1);

    SpaceNodePDE sn;
    TimeNodePDE tn;

    /**************************************************** initial contions ****************************************************/

    for (unsigned int m=0; m<=M; m++)
    {
        sn.j = m; sn.y = m*hy;
        for(unsigned int n=0; n<=N; n++)
        {
            sn.i = n; sn.x = n*hx;
            u00[m][n] = initial1(sn);
        }
    }
    layerInfo(u00);
    IPrinter::printSeperatorLine();

    /**************************************************** border conditions ***************************************************/

    TimeNodePDE tn05; tn05.i = 1; tn05.t = ht*0.5;
    TimeNodePDE tn10; tn10.i = 1; tn10.t = ht;

    SpaceNodePDE sn0, sn1;
    sn0.i = 0; sn0.x = 0.0;
    sn1.i = N; sn1.x = N*hx;
    for (unsigned int m=0; m<=M; m++)
    {
        sn0.j = m; sn0.y = m*hy; u05[m][0] = boundary(sn0, tn05); u10[m][0] = boundary(sn0, tn10);
        sn1.j = m; sn1.y = m*hy; u05[m][N] = boundary(sn1, tn05); u10[m][N] = boundary(sn1, tn10);
    }

    sn0.j = 0;  sn0.y = 0.0;
    sn1.j = M; sn1.y = M*hy;
    for (unsigned int n=0; n<=N; n++)
    {
        sn0.i = n; sn0.x = n*hx; u05[0][n] = boundary(sn0, tn05); u10[0][n] = boundary(sn0, tn10);
        sn1.i = n; sn1.x = n*hx; u05[M][n] = boundary(sn1, tn05); u10[M][n] = boundary(sn1, tn10);
    }

    /**************************************************** border conditions ***************************************************/

    for (unsigned int m=1; m<=M-1; m++)
    {
        sn.j = m; sn.y = m*hy;
        for(unsigned int n=1; n<=N-1; n++)
        {
            sn.i = n; sn.x = n*hx;

            double sum = 0.0;

            double aa__hxhx = (a*a)/(hx*hx);
            double aa__hyhy = (a*a)/(hy*hy);
            sum += aa__hxhx*(u00[m][n-1]-2.0*u00[m][n]+u00[m][n+1]);
            sum += aa__hyhy*(u00[m-1][n]-2.0*u00[m][n]+u00[m+1][n]);
            sum -= lambda*initial2(sn);
            sum += f(sn, tn);

            u05[m][n] = u00[m][n] + initial2(sn)*ht*0.5 + 0.125*ht*ht*sum;
            u10[m][n] = u00[m][n] + initial2(sn)*ht     + 0.500*ht*ht*sum;
        }
    }
    layerInfo(u05);
    IPrinter::printSeperatorLine();
    layerInfo(u10);
    IPrinter::printSeperatorLine();
    IPrinter::printSeperatorLine();

    /**************************************************** initial contions ****************************************************/

    double *ax = (double *) malloc(sizeof(double)*(N-1));
    double *bx = (double *) malloc(sizeof(double)*(N-1));
    double *cx = (double *) malloc(sizeof(double)*(N-1));
    double *dx = (double *) malloc(sizeof(double)*(N-1));
    double *rx = (double *) malloc(sizeof(double)*(N-1));

    double *ay = (double *) malloc(sizeof(double)*(M-1));
    double *by = (double *) malloc(sizeof(double)*(M-1));
    double *cy = (double *) malloc(sizeof(double)*(M-1));
    double *dy = (double *) malloc(sizeof(double)*(M-1));
    double *ry = (double *) malloc(sizeof(double)*(M-1));

    for (unsigned int l=2; l<=L; l++)
    {
        TimeNodePDE tn05; tn05.i = l; tn05.t = l*ht-0.5*ht;
        TimeNodePDE tn10; tn10.i = l; tn10.t = l*ht;

        /**************************************************** border conditions ***************************************************/

        SpaceNodePDE sn0, sn1;
        sn0.i = 0; sn0.x = 0.0;
        sn1.i = N; sn1.x = N*hx;
        for (unsigned int m=0; m<=M; m++)
        {
            sn0.j = m; sn0.y = m*hy; u15[m][0] = boundary(sn0, tn05); u20[m][0] = boundary(sn0, tn10);
            sn1.j = m; sn1.y = m*hy; u15[m][N] = boundary(sn1, tn05); u20[m][N] = boundary(sn1, tn10);
        }

        sn0.j = 0; sn0.y = 0.0;
        sn1.j = M; sn1.y = M*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn0.i = n; sn0.x = n*hx; u15[0][n] = boundary(sn0, tn05); u20[0][n] = boundary(sn0, tn10);
            sn1.i = n; sn1.x = n*hx; u15[M][n] = boundary(sn1, tn05); u20[M][n] = boundary(sn1, tn10);
        }

        /**************************************************** border conditions ***************************************************/

        /**************************************************** x direction apprx ***************************************************/

        for (unsigned int m=1; m<=M-1; m++)
        {
            sn.j = m; sn.y = m*hy;

            for(unsigned int n=1; n<=N-1; n++)
            {
                sn.i = n; sn.x = n*hx;

                ax[n-1] = -(a*a*ht*ht)/(4.0*hx*hx);
                bx[n-1] = 1.0+(a*a*ht*ht)/(2.0*hx*hx)+(lambda*ht*0.75);
                cx[n-1] = -(a*a*ht*ht)/(4.0*hx*hx);

                dx[n-1] = ((a*a*ht*ht)/(4.0*hy*hy))*(u10[m-1][n]-2.0*u10[m][n]+u10[m+1][n]);
                dx[n-1] += (0.25*lambda*ht)*(4.0*u10[m][n]-u05[m][n]);
                dx[n-1] += 2.0*u10[m][n]-u05[m][n];
                dx[n-1] += 0.25*ht*ht*f(sn,tn05);
            }

            dx[0]   -= -(a*a*ht*ht)/(4.0*hx*hx) * u15[m][0];
            dx[N-2] -= -(a*a*ht*ht)/(4.0*hx*hx) * u15[m][N];

            ax[0] = cx[N-2] = 0.0;

            tomasAlgorithm(ax, bx, cx, dx, rx, N-1);
            for (unsigned int n=1; n<=N-1; n++) u15[m][n] = rx[n-1];
        }

        //layerInfo(u15);
        //IPrinter::printSeperatorLine();

        /**************************************************** x direction apprx ***************************************************/

        /**************************************************** y direction apprx ***************************************************/

        for(unsigned int n=1; n<=N-1; n++)
        {
            sn.i = n; sn.x = n*hx;

            for (unsigned int m=1; m<=M-1; m++)
            {
                sn.j = m; sn.y = m*hy;

                ay[m-1] = -(a*a*ht*ht)/(4.0*hy*hy);
                by[m-1] = 1.0+(a*a*ht*ht)/(2.0*hy*hy)+(0.75*lambda*ht);
                cy[m-1] = -(a*a*ht*ht)/(4.0*hy*hy);

                dy[m-1] = ((a*a*ht*ht)/(4.0*hx*hx))*(u15[m-1][n]-2.0*u15[m][n]+u15[m+1][n]);
                dy[m-1] += (0.25*lambda*ht)*(4.0*u15[m][n]-u10[m][n]);
                dy[m-1] += 2.0*u15[m][n]-u10[m][n];
                dy[m-1] += 0.25*ht*ht*f(sn,tn10);
            }

            dy[0]   -= -(a*a*ht*ht)/(4.0*hy*hy) * u20[0][n];
            dy[M-2] -= -(a*a*ht*ht)/(4.0*hy*hy) * u20[M][n];

            ay[0] = cy[M-2] = 0.0;

            tomasAlgorithm(ay, by, cy, dy, ry, M-1);
            for (unsigned int m=1; m<=M-1; m++) u20[m][n] = ry[m-1];
        }

        //layerInfo(u20);
        //IPrinter::printSeperatorLine();

        /**************************************************** y direction apprx ***************************************************/

        for (unsigned int m=0; m<=M; m++)
        {
            for (unsigned int n=0; n<=N; n++)
            {
                u00[m][n] = u10[m][n];
                u05[m][n] = u15[m][n];
                u10[m][n] = u20[m][n];
            }
        }
    }

    layerInfo(u20);
    IPrinter::printSeperatorLine();
}

void Problem2HNM::solveEquation2D2(const Dimension &time, const Dimension &dimx, const Dimension &dimy, DoubleMatrix &u UNUSED_PARAM)
{
    const unsigned int N = dimx.sizeN();
    const unsigned int M = dimy.sizeN();
    const unsigned int L = time.sizeN();

    const double hx = dimx.step();
    const double hy = dimy.step();
    const double ht = time.step();

    DoubleMatrix u00(M+1, N+1);
    DoubleMatrix u05(M+1, N+1);
    DoubleMatrix u10(M+1, N+1);
    DoubleMatrix u15(M+1, N+1);
    DoubleMatrix u20(M+1, N+1);

    SpaceNodePDE sn;
    TimeNodePDE tn;

    /**************************************************** initial contions ****************************************************/

    for (unsigned int m=0; m<=M; m++)
    {
        sn.j = m; sn.y = m*hy;
        for(unsigned int n=0; n<=N; n++)
        {
            sn.i = n; sn.x = n*hx;
            u00[m][n] = initial1(sn);
        }
    }
    layerInfo(u00);
    IPrinter::printSeperatorLine();

    /**************************************************** border conditions ***************************************************/

    TimeNodePDE tn05; tn05.i = 1; tn05.t = ht*0.5;
    TimeNodePDE tn10; tn10.i = 1; tn10.t = ht;

    SpaceNodePDE sn0, sn1;
    sn0.i = 0; sn0.x = 0.0;
    sn1.i = N; sn1.x = N*hx;
    for (unsigned int m=0; m<=M; m++)
    {
        sn0.j = m; sn0.y = m*hy; u05[m][0] = boundary(sn0, tn05); u10[m][0] = boundary(sn0, tn10);
        sn1.j = m; sn1.y = m*hy; u05[m][N] = boundary(sn1, tn05); u10[m][N] = boundary(sn1, tn10);
    }

    sn0.j = 0;  sn0.y = 0.0;
    sn1.j = M; sn1.y = M*hy;
    for (unsigned int n=0; n<=N; n++)
    {
        sn0.i = n; sn0.x = n*hx; u05[0][n] = boundary(sn0, tn05); u10[0][n] = boundary(sn0, tn10);
        sn1.i = n; sn1.x = n*hx; u05[M][n] = boundary(sn1, tn05); u10[M][n] = boundary(sn1, tn10);
    }

    /**************************************************** border conditions ***************************************************/

    for (unsigned int m=1; m<=M-1; m++)
    {
        sn.j = m; sn.y = m*hy;
        for(unsigned int n=1; n<=N-1; n++)
        {
            sn.i = n; sn.x = n*hx;

            double sum = 0.0;

            double aa__hxhx = (a*a)/(hx*hx);
            double aa__hyhy = (a*a)/(hy*hy);
            sum += aa__hxhx*(u00[m][n-1]-2.0*u00[m][n]+u00[m][n+1]);
            sum += aa__hyhy*(u00[m-1][n]-2.0*u00[m][n]+u00[m+1][n]);
            sum -= lambda*initial2(sn);
            sum += f(sn, tn);

            u05[m][n] = u00[m][n] + initial2(sn)*ht*0.5 + 0.125*ht*ht*sum;
            u10[m][n] = u00[m][n] + initial2(sn)*ht     + 0.500*ht*ht*sum;
        }
    }
    layerInfo(u05);
    IPrinter::printSeperatorLine();
    layerInfo(u10);
    IPrinter::printSeperatorLine();
    IPrinter::printSeperatorLine();

    /**************************************************** initial contions ****************************************************/

    double *ax = (double *) malloc(sizeof(double)*(N-1));
    double *bx = (double *) malloc(sizeof(double)*(N-1));
    double *cx = (double *) malloc(sizeof(double)*(N-1));
    double *dx = (double *) malloc(sizeof(double)*(N-1));
    double *rx = (double *) malloc(sizeof(double)*(N-1));

    double *ay = (double *) malloc(sizeof(double)*(M-1));
    double *by = (double *) malloc(sizeof(double)*(M-1));
    double *cy = (double *) malloc(sizeof(double)*(M-1));
    double *dy = (double *) malloc(sizeof(double)*(M-1));
    double *ry = (double *) malloc(sizeof(double)*(M-1));

    for (unsigned int l=2; l<=L; l++)
    {
        TimeNodePDE tn05; tn05.i = l; tn05.t = l*ht-0.5*ht;
        TimeNodePDE tn10; tn10.i = l; tn10.t = l*ht;

        /**************************************************** border conditions ***************************************************/

        SpaceNodePDE sn0, sn1;
        sn0.i = 0; sn0.x = 0.0;
        sn1.i = N; sn1.x = N*hx;
        for (unsigned int m=0; m<=M; m++)
        {
            sn0.j = m; sn0.y = m*hy; u15[m][0] = boundary(sn0, tn05); u20[m][0] = boundary(sn0, tn10);
            sn1.j = m; sn1.y = m*hy; u15[m][N] = boundary(sn1, tn05); u20[m][N] = boundary(sn1, tn10);
        }

        sn0.j = 0; sn0.y = 0.0;
        sn1.j = M; sn1.y = M*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn0.i = n; sn0.x = n*hx; u15[0][n] = boundary(sn0, tn05); u20[0][n] = boundary(sn0, tn10);
            sn1.i = n; sn1.x = n*hx; u15[M][n] = boundary(sn1, tn05); u20[M][n] = boundary(sn1, tn10);
        }

        /**************************************************** border conditions ***************************************************/

        /**************************************************** x direction apprx ***************************************************/

        for (unsigned int m=1; m<=M-1; m++)
        {
            sn.j = m; sn.y = m*hy;

            for(unsigned int n=1; n<=N-1; n++)
            {
                sn.i = n; sn.x = n*hx;

                ax[n-1] = -(a*a*ht*ht)/(hx*hx);
                bx[n-1] = 8.0+2.0*((a*a*ht*ht)/(hx*hx))+3.0*(lambda*ht);
                cx[n-1] = -(a*a*ht*ht)/(hx*hx);

                dx[n-1] = ((a*a*ht*ht)/(hy*hy))*(u10[m-1][n]-2.0*u10[m][n]+u10[m+1][n]);
                dx[n-1] += (lambda*ht)*(4.0*u10[m][n]-u05[m][n]);
                dx[n-1] += 20.0*u10[m][n]-16.0*u05[m][n]+4.0*u00[m][n];
                dx[n-1] += ht*ht*f(sn,tn05);
            }

            dx[0]   -= -(a*a*ht*ht)/(hx*hx) * u15[m][0];
            dx[N-2] -= -(a*a*ht*ht)/(hx*hx) * u15[m][N];

            ax[0] = cx[N-2] = 0.0;

            tomasAlgorithm(ax, bx, cx, dx, rx, N-1);
            for (unsigned int n=1; n<=N-1; n++) u15[m][n] = rx[n-1];
        }

        //layerInfo(u15);
        //IPrinter::printSeperatorLine();

        /**************************************************** x direction apprx ***************************************************/

        /**************************************************** y direction apprx ***************************************************/

        for(unsigned int n=1; n<=N-1; n++)
        {
            sn.i = n; sn.x = n*hx;

            for (unsigned int m=1; m<=M-1; m++)
            {
                sn.j = m; sn.y = m*hy;

                ay[m-1] = -(a*a*ht*ht)/(hy*hy);
                by[m-1] = 8.0+(2.0*a*a*ht*ht)/(hy*hy)+(3.0*lambda*ht);
                cy[m-1] = -(a*a*ht*ht)/(hy*hy);

                dy[m-1] = ((a*a*ht*ht)/(hx*hx))*(u15[m][n-1]-2.0*u15[m][n]+u15[m][n+1]);
                dy[m-1] += (lambda*ht)*(4.0*u15[m][n]-u10[m][n]);
                dy[m-1] += 20.0*u15[m][n]-16.0*u10[m][n]+4.0*u05[m][n];
                dy[m-1] += ht*ht*f(sn,tn10);
            }

            dy[0]   -= -(a*a*ht*ht)/(hy*hy) * u20[0][n];
            dy[M-2] -= -(a*a*ht*ht)/(hy*hy) * u20[M][n];

            ay[0] = cy[M-2] = 0.0;

            tomasAlgorithm(ay, by, cy, dy, ry, M-1);
            for (unsigned int m=1; m<=M-1; m++) u20[m][n] = ry[m-1];
        }

        //layerInfo(u20);
        //IPrinter::printSeperatorLine();

        /**************************************************** y direction apprx ***************************************************/

        for (unsigned int m=0; m<=M; m++)
        {
            for (unsigned int n=0; n<=N; n++)
            {
                u00[m][n] = u10[m][n];
                u05[m][n] = u15[m][n];
                u10[m][n] = u20[m][n];
            }
        }
    }

    layerInfo(u20);
    IPrinter::printSeperatorLine();
}

void Problem2HNM::solveEquation2D3(const Dimension &time, const Dimension &dimx, const Dimension &dimy, DoubleMatrix &u UNUSED_PARAM)
{
    const unsigned int N = dimx.sizeN();
    const unsigned int M = dimy.sizeN();
    const unsigned int L = time.sizeN();

    const double hx = dimx.step();
    const double hy = dimy.step();
    const double ht = time.step();

    DoubleMatrix u00(M+1, N+1);
    DoubleMatrix u05(M+1, N+1);
    DoubleMatrix u10(M+1, N+1);
    DoubleMatrix u15(M+1, N+1);
    DoubleMatrix u20(M+1, N+1);

    SpaceNodePDE sn;
    TimeNodePDE tn;

    /**************************************************** initial contions ****************************************************/

    for (unsigned int m=0; m<=M; m++)
    {
        sn.j = m; sn.y = m*hy;
        for(unsigned int n=0; n<=N; n++)
        {
            sn.i = n; sn.x = n*hx;
            u00[m][n] = initial1(sn);
        }
    }
    layerInfo(u00);
    IPrinter::printSeperatorLine();

    /**************************************************** border conditions ***************************************************/

    TimeNodePDE tn05; tn05.i = 1; tn05.t = ht*0.5;
    TimeNodePDE tn10; tn10.i = 1; tn10.t = ht;

    SpaceNodePDE sn0, sn1;
    sn0.i = 0; sn0.x = 0.0;
    sn1.i = N; sn1.x = N*hx;
    for (unsigned int m=0; m<=M; m++)
    {
        sn0.j = m; sn0.y = m*hy; u05[m][0] = boundary(sn0, tn05); u10[m][0] = boundary(sn0, tn10);
        sn1.j = m; sn1.y = m*hy; u05[m][N] = boundary(sn1, tn05); u10[m][N] = boundary(sn1, tn10);
    }

    sn0.j = 0;  sn0.y = 0.0;
    sn1.j = M; sn1.y = M*hy;
    for (unsigned int n=0; n<=N; n++)
    {
        sn0.i = n; sn0.x = n*hx; u05[0][n] = boundary(sn0, tn05); u10[0][n] = boundary(sn0, tn10);
        sn1.i = n; sn1.x = n*hx; u05[M][n] = boundary(sn1, tn05); u10[M][n] = boundary(sn1, tn10);
    }

    /**************************************************** border conditions ***************************************************/

    for (unsigned int m=1; m<=M-1; m++)
    {
        sn.j = m; sn.y = m*hy;
        for(unsigned int n=1; n<=N-1; n++)
        {
            sn.i = n; sn.x = n*hx;

            double sum = 0.0;

            double aa__hxhx = (a*a)/(hx*hx);
            double aa__hyhy = (a*a)/(hy*hy);
            sum += aa__hxhx*(u00[m][n-1]-2.0*u00[m][n]+u00[m][n+1]);
            sum += aa__hyhy*(u00[m-1][n]-2.0*u00[m][n]+u00[m+1][n]);
            sum -= lambda*initial2(sn);
            sum += f(sn, tn);

            u05[m][n] = u00[m][n] + initial2(sn)*ht*0.5 + 0.125*ht*ht*sum;
            u10[m][n] = u00[m][n] + initial2(sn)*ht     + 0.500*ht*ht*sum;
        }
    }
    layerInfo(u05);
    IPrinter::printSeperatorLine();
    layerInfo(u10);
    IPrinter::printSeperatorLine();
    IPrinter::printSeperatorLine();

    /**************************************************** initial contions ****************************************************/

    double *ax = (double *) malloc(sizeof(double)*(N-1));
    double *bx = (double *) malloc(sizeof(double)*(N-1));
    double *cx = (double *) malloc(sizeof(double)*(N-1));
    double *dx = (double *) malloc(sizeof(double)*(N-1));
    double *rx = (double *) malloc(sizeof(double)*(N-1));

    double *ay = (double *) malloc(sizeof(double)*(M-1));
    double *by = (double *) malloc(sizeof(double)*(M-1));
    double *cy = (double *) malloc(sizeof(double)*(M-1));
    double *dy = (double *) malloc(sizeof(double)*(M-1));
    double *ry = (double *) malloc(sizeof(double)*(M-1));

    for (unsigned int l=2; l<=L; l++)
    {
        TimeNodePDE tn05; tn05.i = l; tn05.t = l*ht-0.5*ht;
        TimeNodePDE tn10; tn10.i = l; tn10.t = l*ht;

        /**************************************************** border conditions ***************************************************/

        SpaceNodePDE sn0, sn1;
        sn0.i = 0; sn0.x = 0.0;
        sn1.i = N; sn1.x = N*hx;
        for (unsigned int m=0; m<=M; m++)
        {
            sn0.j = m; sn0.y = m*hy; u15[m][0] = boundary(sn0, tn05); u20[m][0] = boundary(sn0, tn10);
            sn1.j = m; sn1.y = m*hy; u15[m][N] = boundary(sn1, tn05); u20[m][N] = boundary(sn1, tn10);
        }

        sn0.j = 0; sn0.y = 0.0;
        sn1.j = M; sn1.y = M*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn0.i = n; sn0.x = n*hx; u15[0][n] = boundary(sn0, tn05); u20[0][n] = boundary(sn0, tn10);
            sn1.i = n; sn1.x = n*hx; u15[M][n] = boundary(sn1, tn05); u20[M][n] = boundary(sn1, tn10);
        }

        /**************************************************** border conditions ***************************************************/

        /**************************************************** x direction apprx ***************************************************/

        for (unsigned int m=1; m<=M-1; m++)
        {
            sn.j = m; sn.y = m*hy;

            for(unsigned int n=1; n<=N-1; n++)
            {
                sn.i = n; sn.x = n*hx;

                ax[n-1] = -((a*a*ht*ht)/(hx*hx));
                bx[n-1] = +8.0+2.0*((a*a*ht*ht)/(hx*hx))+(11.0/3.0)*(lambda*ht);
                cx[n-1] = -((a*a*ht*ht)/(hx*hx));

                dx[n-1] = 0.0;
                dx[n-1] += +((a*a*ht*ht)/(hy*hy))*(u10[m-1][n]-2.0*u10[m][n]+u10[m+1][n]);
                dx[n-1] += +4.0*(5.0*u10[m][n]-4.0*u05[m][n]+u00[m][n]);
                dx[n-1] += +((lambda*ht)/3.0)*(18.0*u10[m][n]-9.0*u05[m][n]+2.0*u00[m][n]);
                dx[n-1] += ht*ht*f(sn,tn05);
            }

            dx[0]   -= -((a*a*ht*ht)/(hx*hx)) * u15[m][0];
            dx[N-2] -= -((a*a*ht*ht)/(hx*hx)) * u15[m][N];

            ax[0] = cx[N-2] = 0.0;

            tomasAlgorithm(ax, bx, cx, dx, rx, N-1);
            for (unsigned int n=1; n<=N-1; n++) u15[m][n] = rx[n-1];
        }

        /**************************************************** x direction apprx ***************************************************/

        /**************************************************** y direction apprx ***************************************************/

        for(unsigned int n=1; n<=N-1; n++)
        {
            sn.i = n; sn.x = n*hx;

            for (unsigned int m=1; m<=M-1; m++)
            {
                sn.j = m; sn.y = m*hy;

                ay[m-1] = -((a*a*ht*ht)/(hy*hy));
                by[m-1] = +8.0+2.0*((a*a*ht*ht)/(hy*hy))+(11.0/3.0)*(lambda*ht);
                cy[m-1] = -((a*a*ht*ht)/(hy*hy));

                dy[m-1] = 0.0;
                dy[m-1] += ((a*a*ht*ht)/(hx*hx))*(u15[m][n-1]-2.0*u15[m][n]+u15[m][n+1]);
                dy[m-1] += 4.0*(5.0*u15[m][n]-4.0*u10[m][n]+u05[m][n]);
                dy[m-1] += ((lambda*ht)/3.0)*(18.0*u15[m][n]-9.0*u10[m][n]+2.0*u05[m][n]);
                dy[m-1] += ht*ht*f(sn,tn10);
            }

            dy[0]   -= -((a*a*ht*ht)/(hy*hy)) * u20[0][n];
            dy[M-2] -= -((a*a*ht*ht)/(hy*hy)) * u20[M][n];

            ay[0] = cy[M-2] = 0.0;

            tomasAlgorithm(ay, by, cy, dy, ry, M-1);
            for (unsigned int m=1; m<=M-1; m++) u20[m][n] = ry[m-1];
        }

        /**************************************************** y direction apprx ***************************************************/

        for (unsigned int m=0; m<=M; m++)
        {
            for (unsigned int n=0; n<=N; n++)
            {
                u00[m][n] = u10[m][n];
                u05[m][n] = u15[m][n];
                u10[m][n] = u20[m][n];
            }
        }
    }

    layerInfo(u20);
    IPrinter::printSeperatorLine();
}

void Problem2HNM::solveEquation1D1(const Dimension &time, const Dimension &dimx, DoubleVector &u UNUSED_PARAM)
{
    const unsigned int N = dimx.sizeN();
    const unsigned int L = time.sizeN();

    const double hx = dimx.step();
    const double ht = time.step();

    DoubleVector u00(N+1);
    DoubleVector u10(N+1);
    DoubleVector u20(N+1);

    SpaceNodePDE sn;
    TimeNodePDE tn;

    /**************************************************** initial contions ****************************************************/

    for(unsigned int n=0; n<=N; n++)
    {
        sn.i = n; sn.x = n*hx;
        u00[n] = initial1(sn);
    }
    layerInfo(u00);
    IPrinter::printSeperatorLine();

    /**************************************************** border conditions ***************************************************/

    tn.i = 1; tn.t = ht;
    SpaceNodePDE sn0, sn1;
    sn0.i = 0; sn0.x = 0.0;
    sn1.i = N; sn1.x = N*hx;
    u10[0] = boundary(sn0, tn);
    u10[N] = boundary(sn1, tn);

    /**************************************************** border conditions ***************************************************/

    double aa__hxhx = (a*a)/(hx*hx);
    for(unsigned int n=1; n<=N-1; n++)
    {
        sn.i = n; sn.x = n*hx;
        double sum = 0.0;
        sum += aa__hxhx*(u00[n-1]-2.0*u00[n]+u00[n+1]);
        sum -= lambda*initial2(sn);
        sum += f(sn, tn);
        u10[n] = u00[n] + initial2(sn)*ht + 0.5*ht*ht*sum;
    }
    layerInfo(u10);
    IPrinter::printSeperatorLine();
    IPrinter::printSeperatorLine();

    /**************************************************** initial contions ****************************************************/

    double *ax = (double *) malloc(sizeof(double)*(N-1));
    double *bx = (double *) malloc(sizeof(double)*(N-1));
    double *cx = (double *) malloc(sizeof(double)*(N-1));
    double *dx = (double *) malloc(sizeof(double)*(N-1));
    double *rx = (double *) malloc(sizeof(double)*(N-1));

    /**************************************************** calculating layers *************************************************/

    for (unsigned int l=2; l<=L; l++)
    {
        /*** border ***/
        tn.i = l; tn.t = l*ht;

        sn0.i = 0; sn0.x = 0.0;
        sn1.i = N; sn1.x = N*hx;
        u20[0] = boundary(sn0, tn);
        u20[N] = boundary(sn1, tn);

        /*** border ***/

        for(unsigned int n=1; n<=N-1; n++)
        {
            sn.i = n; sn.x = n*hx;

            ax[n-1] = -(a*a*ht*ht)/(hx*hx);
            bx[n-1] = 1.0+(2.0*a*a*ht*ht)/(hx*hx)+(lambda*ht*1.5);
            cx[n-1] = -(a*a*ht*ht)/(hx*hx);

            dx[n-1] = 2.0*u10[n]-u00[n];
            dx[n-1] += (lambda*ht*0.5)*(4.0*u10[n]-u00[n]);
            dx[n-1] += ht*ht*f(sn,tn);
        }

        dx[0]   -= -(a*a*ht*ht)/(hx*hx) * u20[0];
        dx[N-2] -= -(a*a*ht*ht)/(hx*hx) * u20[N];

        ax[0] = cx[N-2] = 0.0;

        tomasAlgorithm(ax, bx, cx, dx, rx, N-1);
        for (unsigned int n=1; n<=N-1; n++) u20[n] = rx[n-1];

        //layerInfo(u20);
        //IPrinter::printSeperatorLine();

        for (unsigned int n=0; n<=N; n++)
        {
            u00[n] = u10[n];
            u10[n] = u20[n];
        }
    }

    /**************************************************** calculating layers *************************************************/

    layerInfo(u20);
    IPrinter::printSeperatorLine();

    free(ax);
    free(bx);
    free(cx);
    free(dx);
    free(rx);
}










