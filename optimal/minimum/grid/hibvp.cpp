#include "hibvp.h"

IHyperbolicIBVP::~IHyperbolicIBVP() {}

CCIHyperbolicIBVP::~CCIHyperbolicIBVP() {}

void CCIHyperbolicIBVP::calculateD1V2(DoubleVector &u, double a, double lambda) const
{
    const Dimension &dimx = spaceDimension(Dimension::DimensionX);
    const Dimension &time = timeDimension();
    const unsigned int N = static_cast<unsigned int>(dimx.size());
    const unsigned int M = static_cast<unsigned int>(time.size());
    const double hx = dimx.step();
    const double ht = time.step();

    const double alpha = -lambda*(a*a)*((ht*ht)/(hx*hx));
    const double betta = +1.0 + 2.0*lambda*(a*a)*((ht*ht)/(hx*hx));
    const double gamma = +(1.0-2.0*lambda)*(a*a)*((ht*ht)/(hx*hx));
    const double theta = +lambda*(a*a)*((ht*ht)/(hx*hx));
    const double ht_ht = ht*ht;

    u.clear();
    u.resize(N+1);

    DoubleVector u0(N+1);
    DoubleVector u1(N+1);

    double *da = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *db = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *dc = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *dd = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *rx = static_cast<double*>(malloc(sizeof(double)*(N-1)));

    for (unsigned int n=0; n<N-1; n++)
    {
        da[n] = dc[n] = alpha;
        db[n] = betta;
    }
    da[0] = 0.0; dc[N-2] = 0.0;

    SpaceNodePDE sn;
    TimeNodePDE tn0, tn1, tn2;

    for (unsigned int n=0; n<=N; n++)
    {
        sn.i = n; sn.x = n*hx;
        u0[n] = initial1(sn);
    }
    layerInfo(u0, 0);

    tn0.i = 0; tn0.t = tn0.i*ht;
    tn1.i = 1; tn1.t = tn1.i*ht;
    sn.i = 0; sn.x = 0*hx;
    u1[0] = boundary(sn, tn1);
    for (unsigned int n=1; n<=N-1; n++)
    {
        sn.i = n; sn.x = n*hx;
        u1[n] = u0[n] + ht*initial2(sn) + 0.5*ht*ht*(a*a*((u0[n-1]-2.0*u0[n]+u0[n+1]))/(hx*hx)+f(sn, tn0));
    }
    sn.i = N; sn.x = N*hx;
    u1[N] = boundary(sn, tn1);
    layerInfo(u1, 1);

    for (unsigned int m=2; m<=M; m++)
    {
        tn2.i = m+0; tn2.t = tn2.i*ht;
        tn1.i = m-1; tn1.t = tn1.i*ht;
        tn0.i = m-2; tn0.t = tn0.i*ht;

        sn.i = 0; sn.x = 0*hx;
        u[0] = boundary(sn, tn2);
        sn.i = N; sn.x = N*hx;
        u[N] = boundary(sn, tn2);
        for (unsigned int n=1; n<=N-1; n++)
        {
            sn.i = n; sn.x = n*hx;
            dd[n-1] = gamma*(u1[n-1]-2.0*u1[n]+u1[n+1])
                    + theta*(u0[n-1]-2.0*u0[n]+u0[n+1])
                    + 2.0*u1[n] - u0[n];
            //dd[n-1] += ht_ht*(lambda*f(sn, tn2) + (1.0-2.0*lambda)*f(sn, tn1) + lambda*f(sn, tn0));
            dd[n-1] += ht_ht*f(sn, tn1);
        }
        dd[0]   -= alpha*u[0];
        dd[N-2] -= alpha*u[N];

        tomasAlgorithm(da, db, dc, dd, rx, N-1);

        for (unsigned int n=1; n<=N-1; n++) u[n] = rx[n-1];

        for (unsigned int n=0; n<=N; n++)
        {
            u0[n] = u1[n];
            u1[n] = u[n];
        }

        layerInfo(u, m);
    }

    free(rx);
    free(dd);
    free(dc);
    free(db);
    free(da);
}

void CCIHyperbolicIBVP::calculateD1V1(DoubleVector &u, double a) const
{
    const Dimension &dimx = spaceDimension(Dimension::DimensionX);
    const Dimension &time = timeDimension();
    const unsigned int N = static_cast<const unsigned int>(dimx.size());
    const unsigned int M = static_cast<const unsigned int>(time.size());
    const double hx = dimx.step();
    const double ht = time.step();

    const double alpha = -(a*a)*((ht*ht)/(hx*hx));
    const double betta = +(1.0 + 2.0*(a*a)*((ht*ht)/(hx*hx)));
    const double ht_ht = ht*ht;

    u.clear();
    u.resize(N+1);

    DoubleVector u0(N+1);
    DoubleVector u1(N+1);

    double *da = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *db = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *dc = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *dd = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *rx = static_cast<double*>(malloc(sizeof(double)*(N-1)));

    for (unsigned int n=0; n<N-1; n++)
    {
        da[n] = dc[n] = alpha;
        db[n] = betta;
    }
    da[0] = 0.0; dc[N-2] = 0.0;

    SpaceNodePDE sn;
    TimeNodePDE tn0, tn1;

    for (unsigned int n=0; n<=N; n++)
    {
        sn.i = n; sn.x = n*hx;
        u0[n] = initial1(sn);
    }
    layerInfo(u0, 0);

    tn0.i = 0; tn0.t = tn0.i*ht;
    tn1.i = 1; tn1.t = tn1.i*ht;

    sn.i = 0; sn.x = 0*hx; u1[0] = boundary(sn, tn1);
    for (unsigned int n=1; n<=N-1; n++)
    {
        sn.i = n; sn.x = n*hx;
        u1[n] = u0[n] + ht*initial2(sn) + 0.5*ht*ht*(a*a*((u0[n-1]-2.0*u0[n]+u0[n+1]))/(hx*hx)+f(sn, tn0));
    }
    sn.i = N; sn.x = N*hx; u1[N] = boundary(sn, tn1);
    layerInfo(u0, 1);

    TimeNodePDE tn;
    for (unsigned int m=2; m<=M; m++)
    {
        tn.i = m; tn.t = tn.i*ht;

        sn.i = 0; sn.x = 0*hx; u[0] = boundary(sn, tn);
        sn.i = N; sn.x = N*hx; u[N] = boundary(sn, tn);

        for (unsigned int n=1; n<=N-1; n++)
        {
            sn.i = n; sn.x = n*hx;
            dd[n-1] = 2.0*u1[n] - u0[n] + ht_ht*f(sn, tn);
        }
        dd[0]   -= alpha*u[0];
        dd[N-2] -= alpha*u[N];

        tomasAlgorithm(da, db, dc, dd, rx, N-1);

        for (unsigned int n=1; n<=N-1; n++) u[n] = rx[n-1];

        for (unsigned int n=0; n<=N; n++)
        {
            u0[n] = u1[n];
            u1[n] = u[n];
        }
        layerInfo(u, m);
    }
}

void CCIHyperbolicIBVP::calculateD2V1(DoubleMatrix &u, double a) const
{
    const Dimension &dimX = spaceDimension(Dimension::DimensionX);
    const Dimension &dimY = spaceDimension(Dimension::DimensionX);
    const Dimension &time = timeDimension();

    const unsigned int N = static_cast<unsigned int>( dimX.size() );
    const unsigned int M = static_cast<unsigned int>( dimY.size() );
    const unsigned int L = static_cast<unsigned int>( time.size() );

    const double hx = dimX.step();
    const double hy = dimY.step();
    const double ht = time.step();

    const double m_aa_htht__hxhx = -(a*a)*((ht*ht)/(hx*hx));
    const double b_aa_htht__hxhx = +(2.0 + 2.0*(a*a)*((ht*ht)/(hx*hx)));
    const double p_aa_htht__hyhy = +(a*a)*((ht*ht)/(hy*hy));

    const double m_aa_htht__hyhy = -(a*a)*((ht*ht)/(hy*hy));
    const double b_aa_htht__hyhy = +(2.0 + 2.0*(a*a)*((ht*ht)/(hy*hy)));
    const double p_aa_htht__hxhx = +(a*a)*((ht*ht)/(hx*hx));
    const double ht_ht = ht*ht;

    const double aa__hxhx = (a*a)/(hx*hx);
    const double aa__hyhy = (a*a)/(hy*hy);

    double *ax = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *bx = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *cx = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *dx = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *rx = static_cast<double*>(malloc(sizeof(double)*(N-1)));

    for (unsigned int n=0; n<=N-2; n++)
    {
        ax[n] = m_aa_htht__hxhx;
        bx[n] = b_aa_htht__hxhx;
        cx[n] = m_aa_htht__hxhx;
    }
    ax[0] = 0.0; cx[N-2] = 0.0;

    double *ay = static_cast<double*>(malloc(sizeof(double)*(M-1)));
    double *by = static_cast<double*>(malloc(sizeof(double)*(M-1)));
    double *cy = static_cast<double*>(malloc(sizeof(double)*(M-1)));
    double *dy = static_cast<double*>(malloc(sizeof(double)*(M-1)));
    double *ry = static_cast<double*>(malloc(sizeof(double)*(M-1)));

    for (unsigned int m=0; m<=M-2; m++)
    {
        ay[m] = m_aa_htht__hyhy;
        by[m] = b_aa_htht__hyhy;
        cy[m] = m_aa_htht__hyhy;
    }
    ay[0] = 0.0; cy[M-2] = 0.0;

    DoubleMatrix u00(M+1, N+1);
    DoubleMatrix u10(M+1, N+1);
    DoubleMatrix u15(M+1, N+1);
    DoubleMatrix u20(M+1, N+1);

    /***********************************************************************************************/

    SpaceNodePDE sn;
    for (unsigned int m=0; m<=M; m++)
    {
        sn.j = m; sn.y = sn.j*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = n; sn.x = sn.i*hx;
            u00[m][n] = initial1(sn);
        }
    }
    layerInfo(u00, 0);

    /***********************************************************************************************/

    TimeNodePDE tn00; tn00.i = 0; tn00.t = 0.0;
    TimeNodePDE tn10; tn10.i = 1; tn10.t = ht;
    SpaceNodePDE sn0;
    SpaceNodePDE sn1;
    sn0.i = 0; sn0.x = 0.0; sn1.i = N; sn1.x = hx*N;
    for (unsigned int m=0; m<=M; m++)
    {
        sn0.j = m; sn0.y = m*hy; u10[m][0] = boundary(sn0, tn10);
        sn1.j = m; sn1.y = m*hy; u10[m][N] = boundary(sn1, tn10);
    }
    sn0.j = 0; sn0.y = 0.0; sn1.j = M; sn1.y = hy*M;
    for (unsigned int n=0; n<=N; n++)
    {
        sn0.i = n; sn0.x = n*hx; u10[0][n] = boundary(sn0, tn10);
        sn1.i = n; sn1.x = n*hx; u10[M][n] = boundary(sn1, tn10);
    }

    for (unsigned int m=1; m<=M-1; m++)
    {
        sn.j = m; sn.y = sn.j*hy;
        for (unsigned int n=1; n<=N-1; n++)
        {
            sn.i = n; sn.x = sn.i*hx;
            u10[m][n] = u00[m][n] + ht*initial2(sn);
            u10[m][n] += 0.5*ht*ht*(aa__hxhx*(u00[m][n-1]-2.0*u00[m][n]+u00[m][n+1])
                    +aa__hyhy*(u00[m-1][n]-2.0*u00[m][n]+u00[m+1][n])
                    +f(sn,tn00));
        }
    }
    layerInfo(u10, 1);

    /***********************************************************************************************/

    for (unsigned int ln=1; ln<=L-1; ln++)
    {
        TimeNodePDE tn;   tn.i = ln;   tn.t = tn.i*ht;
        TimeNodePDE tn15; tn15.i = ln; tn15.t = ln*ht+0.5*ht;
        TimeNodePDE tn20; tn20.i = ln; tn20.t = ln*ht+1.0*ht;

        SpaceNodePDE sn0;
        SpaceNodePDE sn1;

        /**************************************************** border conditions ***************************************************/

        sn0.i = 0; sn0.x = 0.0;
        sn1.i = N; sn1.x = hx*N;
        for (unsigned int m=0; m<=M; m++)
        {
            sn0.j = m; sn0.y = m*hy;
            u15[m][0] = boundary(sn0, tn15);
            u20[m][0] = boundary(sn0, tn20);

            sn1.j = m; sn1.y = m*hy;
            u15[m][N] = boundary(sn1, tn15);
            u20[m][N] = boundary(sn1, tn20);
        }

        sn0.j = 0; sn0.y = 0.0;
        sn1.j = M; sn1.y = hy*M;
        for (unsigned int n=0; n<=N; n++)
        {
            sn0.i = n; sn0.x = n*hx;
            u15[0][n] = boundary(sn0, tn15);
            u20[0][n] = boundary(sn0, tn20);

            sn1.i = n; sn1.x = n*hx;
            u15[M][n] = boundary(sn1, tn15);
            u20[M][n] = boundary(sn1, tn20);
        }

        /**************************************************** border conditions ***************************************************/

        /**************************************************** x direction apprx ***************************************************/

        for (unsigned int m=1; m<=M-1; m++)
        {
            sn.j = m; sn.y = m*hy;

            for (unsigned int n=1; n<=N-1; n++)
            {
                sn.i = n; sn.x = n*hx;

                dx[n-1] = 0.0;
                //if (m>0 && m<M)  dx[n-1] = p_aa_htht__hyhy*(u10[m-1][n] - 2.0*u10[m][n]   + u10[m+1][n]);
                //else if (m == 0) dx[n-1] = p_aa_htht__hyhy*(u10[0][n]   - 2.0*u10[1][n]   + u10[2][n]);
                //else if (m == M) dx[n-1] = p_aa_htht__hyhy*(u10[M-2][n] - 2.0*u10[M-1][n] + u10[M][n]);
                dx[n-1] += p_aa_htht__hyhy*(u10[m-1][n] - 2.0*u10[m][n]   + u10[m+1][n]);
                dx[n-1] += 2.0*u10[m][n] + (u10[m][n]-u00[m][n]) + ht_ht*f(sn, tn);
            }
            dx[0]   -= m_aa_htht__hxhx * u15[m][0];
            dx[N-2] -= m_aa_htht__hxhx * u15[m][N];
            tomasAlgorithm(ax, bx, cx, dx, rx, N-1);
            for (unsigned int n=1; n<=N-1; n++) u15[m][n] = rx[n-1];
        }

        /**************************************************** x direction apprx ***************************************************/

        /**************************************************** y direction apprx ***************************************************/

        for (unsigned int n=1; n<=N-1; n++)
        {
            sn.i = n; sn.x = n*hx;
            for (unsigned int m=1; m<=M-1; m++)
            {
                sn.j = m; sn.y = m*hy;
                dy[m-1] = 0.0;
                //if (n>0 && n<N) dy[m-1] = p_aa_htht__hxhx*(u15[m][n-1] - 2.0*u15[m][n]   + u15[m][n+1]);
                //else if (n==0)  dy[m-1] = p_aa_htht__hxhx*(u15[m][0]   - 2.0*u15[m][1]   + u15[m][2]);
                //else if (n==N)  dy[m-1] = p_aa_htht__hxhx*(u15[m][N-2] - 2.0*u15[m][N-1] + u15[m][N]);
                dy[m-1] += p_aa_htht__hxhx*(u15[m][n-1] - 2.0*u15[m][n]   + u15[m][n+1]);
                dy[m-1] += 2.0*u15[m][n] + (u10[m][n]-u00[m][n]) + ht_ht*f(sn, tn);
            }
            dy[0]   -= m_aa_htht__hyhy * u20[0][n];
            dy[M-2] -= m_aa_htht__hyhy * u20[M][n];
            tomasAlgorithm(ay, by, cy, dy, ry, M-1);
            for (unsigned int m=1; m<=M-1; m++) u20[m][n] = ry[m-1];
        }

        /**************************************************** y direction apprx ***************************************************/

        //layerInfo(u15, 2*ln-1);

        for (unsigned int m=0; m<=M; m++)
        {
            for (unsigned int n=0; n<=N; n++)
            {
                u00[m][n] = u10[m][n];
                u10[m][n] = u20[m][n];
            }
        }

        layerInfo(u20, ln+1);
    }

    u00.clear();
    u10.clear();
    u15.clear();
    u20.clear();

    free(ry);
    free(dy);
    free(cy);
    free(by);
    free(ay);

    free(rx);
    free(dx);
    free(cx);
    free(bx);
    free(ax);
}

void CCIHyperbolicIBVP::calculateD2V2(DoubleMatrix &u, double a, double lambda) const
{}

CC1IHyperbolicIBVP::~CC1IHyperbolicIBVP()
{}

void CC1IHyperbolicIBVP::initial_calculate(DoubleMatrix &u00, DoubleMatrix &u10, unsigned int N, double hx, unsigned int M, double hy, double ht, double a, double alpha) const
{
    const double aa__hxhx = (a*a)/(hx*hx);
    const double aa__hyhy = (a*a)/(hy*hy);
    const double htht_05 = 0.5*ht*ht;

    /***********************************************************************************************/

    SpaceNodePDE sn;
    for (unsigned int m=0; m<=M; m++)
    {
        sn.j = static_cast<int>(m); sn.y = sn.j*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = static_cast<int>(n); sn.x = sn.i*hx;
            u00[m][n] = initial1(sn);
        }
    }
    layerInfo(u00, 0);

    /***********************************************************************************************/

    TimeNodePDE tn00; tn00.i = 0; tn00.t = 0.0;
    TimeNodePDE tn10; tn10.i = 1; tn10.t = ht;
    SpaceNodePDE sn0;
    SpaceNodePDE sn1;
    sn0.i = 0; sn0.x = 0.0; sn1.i = static_cast<int>(N); sn1.x = N*hx;
    for (unsigned int m=0; m<=M; m++)
    {
        sn0.j = static_cast<int>(m); sn0.y = m*hy; u10[m][0] = boundary(sn0, tn10);
        sn1.j = static_cast<int>(m); sn1.y = m*hy; u10[m][N] = boundary(sn1, tn10);
    }
    sn0.j = 0; sn0.y = 0.0; sn1.j = static_cast<int>(M); sn1.y = M*hy;
    for (unsigned int n=0; n<=N; n++)
    {
        sn0.i = static_cast<int>(n); sn0.x = n*hx; u10[0][n] = boundary(sn0, tn10);
        sn1.i = static_cast<int>(n); sn1.x = n*hx; u10[M][n] = boundary(sn1, tn10);
    }

    for (unsigned int m=1; m<=M-1; m++)
    {
        sn.j = static_cast<int>(m); sn.y = m*hy;
        for (unsigned int n=1; n<=N-1; n++)
        {
            sn.i = static_cast<int>(n); sn.x = n*hx;
            u10[m][n] = u00[m][n] + ht*initial2(sn);
            u10[m][n] += htht_05*(aa__hxhx*(u00[m][n-1]-2.0*u00[m][n]+u00[m][n+1])+aa__hyhy*(u00[m-1][n]-2.0*u00[m][n]+u00[m+1][n])+f(sn,tn00)-alpha*initial2(sn));
        }
    }
    layerInfo(u10, 1);

    /***********************************************************************************************/
}

void CC1IHyperbolicIBVP::border__calculate(DoubleMatrix &u15, DoubleMatrix &u20, unsigned int N, double hx, unsigned int M, double hy, const TimeNodePDE &tn15, const TimeNodePDE &tn20) const
{
    SpaceNodePDE sn0;
    SpaceNodePDE sn1;

    sn0.i = 0; sn0.x = 0.0; sn1.i = static_cast<int>(N); sn1.x = hx*N;
    for (unsigned int m=0; m<=M; m++)
    {
        sn0.j = static_cast<int>(m); sn0.y = m*hy; u15[m][0] = boundary(sn0, tn15); u20[m][0] = boundary(sn0, tn20);
        sn1.j = static_cast<int>(m); sn1.y = m*hy; u15[m][N] = boundary(sn1, tn15); u20[m][N] = boundary(sn1, tn20);
    }

    sn0.j = 0; sn0.y = 0.0; sn1.j = static_cast<int>(M); sn1.y = hy*M;
    for (unsigned int n=0; n<=N; n++)
    {
        sn0.i = static_cast<int>(n); sn0.x = n*hx; u15[0][n] = boundary(sn0, tn15); u20[0][n] = boundary(sn0, tn20);
        sn1.i = static_cast<int>(n); sn1.x = n*hx; u15[M][n] = boundary(sn1, tn15); u20[M][n] = boundary(sn1, tn20);
    }
}

void CC1IHyperbolicIBVP::initial_calculate(DoubleVector &u00, DoubleVector &u10, unsigned int N, double hx, double ht, double a, double alpha) const
{
    const double aa__hxhx = (a*a)/(hx*hx);
    const double htht_05 = 0.5*ht*ht;

    /***********************************************************************************************/

    SpaceNodePDE sn; sn.j = 0; sn.y = 0.0;
    for (unsigned int n=0; n<=N; n++)
    {
        sn.i = static_cast<int>(n); sn.x = sn.i*hx;
        u00[n] = initial1(sn);
    }
    layerInfo(u00, 0);

    /***********************************************************************************************/

    TimeNodePDE tn00; tn00.i = 0; tn00.t = 0.0;
    TimeNodePDE tn10; tn10.i = 1; tn10.t = ht;
    sn.i = static_cast<int>(0); sn.x = 0.0;   u10[0] = boundary(sn, tn10);
    sn.i = static_cast<int>(N); sn.x = N*hx;  u10[N] = boundary(sn, tn10);
    for (unsigned int n=1; n<=N-1; n++)
    {
        sn.i = static_cast<int>(n); sn.x = n*hx;
        u10[n] = u00[n] + ht*initial2(sn);
        u10[n] += htht_05*(aa__hxhx*(u00[n-1]-2.0*u00[n]+u00[n+1])+f(sn,tn00));
        u10[n] += htht_05*(-alpha*initial2(sn));
    }
    layerInfo(u10, 1);

    /***********************************************************************************************/
}

void CC1IHyperbolicIBVP::explicit_calculate_D1V1(DoubleVector &, double, double) const {}

void CC1IHyperbolicIBVP::implicit_calculate_D1V1(DoubleVector &, double, double) const {}

void CC1IHyperbolicIBVP::implicit_calculate_D1V2(DoubleVector &, double, double, double) const {}

void CC1IHyperbolicIBVP::explicit_calculate_D2V1(DoubleMatrix &u, double a, double alpha) const
{
    const Dimension &dimX = spaceDimension(Dimension::DimensionX);
    const Dimension &dimY = spaceDimension(Dimension::DimensionX);
    const Dimension &time = timeDimension();

    const unsigned int N = static_cast<unsigned int>( dimX.size() );
    const unsigned int M = static_cast<unsigned int>( dimY.size() );
    const unsigned int L = static_cast<unsigned int>( time.size() );

    const double hx = dimX.step();
    const double hy = dimY.step();
    const double ht = time.step();

    const double alpha_ht_05 = alpha*ht*0.5;
    const double inv__alpha_ht = 1.0/(1.0 + alpha_ht_05);
    const double aa_htht__hxhx = ((a*a*ht*ht)/(hx*hx));
    const double aa_htht__hyhy = ((a*a*ht*ht)/(hy*hy));

    DoubleMatrix u00(M+1, N+1);
    DoubleMatrix u10(M+1, N+1);
    DoubleMatrix u20(M+1, N+1);

    /**********************************************************************/
    initial_calculate(u00, u10, N, hx, M, hy, ht, a, alpha);
    /**********************************************************************/

    for (unsigned int ln=2; ln<=L; ln++)
    {
        TimeNodePDE tn20; tn20.i = ln+0; tn20.t = tn20.i*ht;
        TimeNodePDE tn10; tn10.i = ln-1; tn10.t = tn10.i*ht;

        /**************************** border conditions ****************************/
        SpaceNodePDE sn0, sn1;
        sn0.i = static_cast<int>(0); sn0.x = 0*hx;
        sn1.i = static_cast<int>(N); sn1.x = N*hx;
        for (unsigned int m=0; m<=M; m++)
        {
            sn0.j = static_cast<int>(m); sn0.y = m*hy; u20[m][0] = boundary(sn0, tn20);
            sn1.j = static_cast<int>(m); sn1.y = m*hy; u20[m][N] = boundary(sn1, tn20);
        }

        sn0.j = static_cast<int>(0); sn0.y = 0*hy;
        sn1.j = static_cast<int>(M); sn1.y = M*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn0.i = static_cast<int>(n); sn0.x = n*hx; u20[0][n] = boundary(sn0, tn20);
            sn1.i = static_cast<int>(n); sn1.x = n*hx; u20[M][n] = boundary(sn1, tn20);
        }
        /**************************** border conditions ****************************/

        SpaceNodePDE sn;
        for (unsigned int m=1; m<=M-1; m++)
        {
            sn.j = static_cast<int>(m); sn.y = m*hy;
            for (unsigned int n=1; n<=N-1; n++)
            {
                sn.i = static_cast<int>(n); sn.x = n*hx;
                u20[m][n] = inv__alpha_ht * (aa_htht__hxhx*(u10[m][n-1]-2.0*u10[m][n]+u10[m][n+1])
                                           + aa_htht__hyhy*(u10[m-1][n]-2.0*u10[m][n]+u10[m+1][n])
                                           + alpha_ht_05*u00[m][n] + 2.0*u10[m][n] - u00[m][n] + ht*ht*f(sn,tn10));
            }
        }
        layerInfo(u20, ln);

        for (unsigned int m=0; m<=M; m++)
        {
            for (unsigned int n=0; n<=N; n++)
            {
                u00[m][n] = u10[m][n];
                u10[m][n] = u20[m][n];
            }
        }
    }

    for (unsigned int m=0; m<=M; m++)
    {
        for (unsigned int n=0; n<=N; n++)
        {
            u[m][n] = u20[m][n];
        }
    }
}

void CC1IHyperbolicIBVP::implicit_calculate_D2V1(DoubleMatrix &, double a, double alpha) const
{
    const Dimension &dimX = spaceDimension(Dimension::DimensionX);
    const Dimension &dimY = spaceDimension(Dimension::DimensionX);
    const Dimension &time = timeDimension();

    const unsigned int N = static_cast<unsigned int>( dimX.size() );
    const unsigned int M = static_cast<unsigned int>( dimY.size() );
    const unsigned int L = static_cast<unsigned int>( time.size() );

    const double hx = dimX.step();
    const double hy = dimY.step();
    const double ht = time.step();

    const double m_aa_htht__hxhx = -(a*a)*((ht*ht)/(hx*hx));
    const double b_aa_htht__hxhx = +(2.0 + 2.0*(a*a)*((ht*ht)/(hx*hx)) + alpha*ht);
    const double p_aa_htht__hyhy = +(a*a)*((ht*ht)/(hy*hy));

    const double m_aa_htht__hyhy = -(a*a)*((ht*ht)/(hy*hy));
    const double b_aa_htht__hyhy = +(2.0 + 2.0*(a*a)*((ht*ht)/(hy*hy)) + alpha*ht);
    const double p_aa_htht__hxhx = +(a*a)*((ht*ht)/(hx*hx));

    const double ht_ht = ht*ht;
    const double alpha_ht = alpha*ht;

    double *ax = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *bx = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *cx = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *dx = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *rx = static_cast<double*>(malloc(sizeof(double)*(N-1)));

    for (unsigned int n=0; n<=N-2; n++)
    {
        ax[n] = m_aa_htht__hxhx;
        bx[n] = b_aa_htht__hxhx;
        cx[n] = m_aa_htht__hxhx;
    }
    ax[0] = 0.0; cx[N-2] = 0.0;

    double *ay = static_cast<double*>(malloc(sizeof(double)*(M-1)));
    double *by = static_cast<double*>(malloc(sizeof(double)*(M-1)));
    double *cy = static_cast<double*>(malloc(sizeof(double)*(M-1)));
    double *dy = static_cast<double*>(malloc(sizeof(double)*(M-1)));
    double *ry = static_cast<double*>(malloc(sizeof(double)*(M-1)));

    for (unsigned int m=0; m<=M-2; m++)
    {
        ay[m] = m_aa_htht__hyhy;
        by[m] = b_aa_htht__hyhy;
        cy[m] = m_aa_htht__hyhy;
    }
    ay[0] = 0.0; cy[M-2] = 0.0;

    DoubleMatrix u00(M+1, N+1);
    DoubleMatrix u10(M+1, N+1);
    DoubleMatrix u15(M+1, N+1);
    DoubleMatrix u20(M+1, N+1);

    /***********************************************************************************************/
    initial_calculate(u00, u10, N, hx, M, hy, ht, a, alpha);
    /***********************************************************************************************/

    SpaceNodePDE sn;
    for (unsigned int ln=1; ln<=L-1; ln++)
    {
        TimeNodePDE tn10; tn10.i = ln; tn10.t = ln*ht;
        TimeNodePDE tn15; tn15.i = ln; tn15.t = ln*ht+0.5*ht;
        TimeNodePDE tn20; tn20.i = ln; tn20.t = ln*ht+1.0*ht;

        /**************************************************** border conditions ***************************************************/
        border__calculate(u15, u20, N, hx, M, hy, tn15, tn20);
        /**************************************************** border conditions ***************************************************/
        /**************************************************** x direction apprx ***************************************************/
        for (unsigned int m=1; m<=M-1; m++)
        {
            sn.j = static_cast<int>(m); sn.y = m*hy;
            for (unsigned int n=1; n<=N-1; n++)
            {
                sn.i = static_cast<int>(n); sn.x = n*hx;
                dx[n-1] = 0.0;
                dx[n-1] += p_aa_htht__hyhy*(u10[m-1][n] - 2.0*u10[m][n] + u10[m+1][n]);
                dx[n-1] += 2.0*u10[m][n] + (u10[m][n] - u00[m][n]);
                dx[n-1] += alpha_ht*(u10[m][n] - 0.5*(u10[m][n]-u00[m][n]));
                dx[n-1] += ht_ht*f(sn, tn10);
            }
            dx[0]   -= m_aa_htht__hxhx * u15[m][0];
            dx[N-2] -= m_aa_htht__hxhx * u15[m][N];
            tomasAlgorithm(ax, bx, cx, dx, rx, N-1);
            for (unsigned int n=1; n<=N-1; n++) u15[m][n] = rx[n-1];
        }
        /**************************************************** x direction apprx ***************************************************/
        /**************************************************** y direction apprx ***************************************************/
        for (unsigned int n=1; n<=N-1; n++)
        {
            sn.i = static_cast<int>(n); sn.x = n*hx;
            for (unsigned int m=1; m<=M-1; m++)
            {
                sn.j = static_cast<int>(m); sn.y = m*hy;
                dy[m-1] = 0.0;
                dy[m-1] += p_aa_htht__hxhx*(u15[m][n-1] - 2.0*u15[m][n] + u15[m][n+1]);
                dy[m-1] += 2.0*u15[m][n] + (u10[m][n]-u00[m][n]);
                dy[m-1] += alpha_ht*(u15[m][n] - 0.5*(u10[m][n]-u00[m][n]));
                dy[m-1] += ht_ht*f(sn, tn10);
            }
            dy[0]   -= m_aa_htht__hyhy * u20[0][n];
            dy[M-2] -= m_aa_htht__hyhy * u20[M][n];
            tomasAlgorithm(ay, by, cy, dy, ry, M-1);
            for (unsigned int m=1; m<=M-1; m++) u20[m][n] = ry[m-1];
        }
        /**************************************************** y direction apprx ***************************************************/
        for (unsigned int m=0; m<=M; m++)
        {
            for (unsigned int n=0; n<=N; n++)
            {
                u00[m][n] = u10[m][n];
                u10[m][n] = u20[m][n];
            }
        }

        layerInfo(u20, ln+1);
    }

    u00.clear();
    u10.clear();
    u15.clear();
    u20.clear();

    free(ry);
    free(dy);
    free(cy);
    free(by);
    free(ay);

    free(rx);
    free(dx);
    free(cx);
    free(bx);
    free(ax);
}

void CC1IHyperbolicIBVP::implicit_calculate_D2V3(DoubleMatrix &, double a, double alpha) const
{
    const Dimension &dimX = spaceDimension(Dimension::DimensionX);
    const Dimension &dimY = spaceDimension(Dimension::DimensionX);
    const Dimension &time = timeDimension();

    const unsigned int N = static_cast<unsigned int>( dimX.size() );
    const unsigned int M = static_cast<unsigned int>( dimY.size() );
    const unsigned int L = static_cast<unsigned int>( time.size() );

    const double hx = dimX.step();
    const double hy = dimY.step();
    const double ht = time.step();

    const double m_aa_htht__hxhx = -(a*a)*((ht*ht)/(hx*hx));
    const double b_aa_htht__hxhx = +(2.0 + 2.0*(a*a)*((ht*ht)/(hx*hx)) + 2.0*alpha*ht);
    const double p_aa_htht__hyhy = +(a*a)*((ht*ht)/(hy*hy));

    const double m_aa_htht__hyhy = -(a*a)*((ht*ht)/(hy*hy));
    const double b_aa_htht__hyhy = +(2.0 + 2.0*(a*a)*((ht*ht)/(hy*hy)) + 2.0*alpha*ht);
    const double p_aa_htht__hxhx = +(a*a)*((ht*ht)/(hx*hx));

    const double ht_ht = ht*ht;
    const double alpha_ht = alpha*ht;

    double *ax = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *bx = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *cx = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *dx = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *rx = static_cast<double*>(malloc(sizeof(double)*(N-1)));

    for (unsigned int n=0; n<=N-2; n++)
    {
        ax[n] = m_aa_htht__hxhx;
        bx[n] = b_aa_htht__hxhx;
        cx[n] = m_aa_htht__hxhx;
    }
    ax[0] = 0.0; cx[N-2] = 0.0;

    double *ay = static_cast<double*>(malloc(sizeof(double)*(M-1)));
    double *by = static_cast<double*>(malloc(sizeof(double)*(M-1)));
    double *cy = static_cast<double*>(malloc(sizeof(double)*(M-1)));
    double *dy = static_cast<double*>(malloc(sizeof(double)*(M-1)));
    double *ry = static_cast<double*>(malloc(sizeof(double)*(M-1)));

    for (unsigned int m=0; m<=M-2; m++)
    {
        ay[m] = m_aa_htht__hyhy;
        by[m] = b_aa_htht__hyhy;
        cy[m] = m_aa_htht__hyhy;
    }
    ay[0] = 0.0; cy[M-2] = 0.0;

    DoubleMatrix u00(M+1, N+1);
    DoubleMatrix u10(M+1, N+1);
    DoubleMatrix u15(M+1, N+1);
    DoubleMatrix u20(M+1, N+1);

    /***********************************************************************************************/
    initial_calculate(u00, u10, N, hx, M, hy, ht, a, alpha);
    /***********************************************************************************************/

    SpaceNodePDE sn;
    for (unsigned int ln=1; ln<=L-1; ln++)
    {
        TimeNodePDE tn10; tn10.i = ln; tn10.t = ln*ht;
        TimeNodePDE tn15; tn15.i = ln; tn15.t = ln*ht+0.5*ht;
        TimeNodePDE tn20; tn20.i = ln; tn20.t = ln*ht+1.0*ht;
        /**************************************************** border conditions ***************************************************/
        border__calculate(u15, u20, N, hx, M, hy, tn15, tn20);
        /**************************************************** border conditions ***************************************************/
        /**************************************************** x direction apprx ***************************************************/
        for (unsigned int m=1; m<=M-1; m++)
        {
            sn.j = static_cast<int>(m); sn.y = m*hy;
            for (unsigned int n=1; n<=N-1; n++)
            {
                sn.i = static_cast<int>(n); sn.x = n*hx;
                dx[n-1] = 0.0;
                dx[n-1] += p_aa_htht__hyhy*(u10[m-1][n] - 2.0*u10[m][n] + u10[m+1][n]);
                dx[n-1] += 2.0*u10[m][n] + (u10[m][n] - u00[m][n]);
                dx[n-1] += 2.0*alpha_ht*u10[m][n];
                dx[n-1] += ht_ht*f(sn, tn15);
            }
            dx[0]   -= m_aa_htht__hxhx * u15[m][0];
            dx[N-2] -= m_aa_htht__hxhx * u15[m][N];
            tomasAlgorithm(ax, bx, cx, dx, rx, N-1);
            for (unsigned int n=1; n<=N-1; n++) u15[m][n] = rx[n-1];
        }
        /**************************************************** x direction apprx ***************************************************/
        /**************************************************** y direction apprx ***************************************************/
        for (unsigned int n=1; n<=N-1; n++)
        {
            sn.i = static_cast<int>(n); sn.x = n*hx;
            for (unsigned int m=1; m<=M-1; m++)
            {
                sn.j = static_cast<int>(m); sn.y = m*hy;
                dy[m-1] = 0.0;
                dy[m-1] += p_aa_htht__hxhx*(u15[m][n-1] - 2.0*u15[m][n] + u15[m][n+1]);
                dy[m-1] += 2.0*u15[m][n] + (u10[m][n]-u00[m][n]);
                dy[m-1] += 2.0*alpha_ht*u15[m][n];
                dy[m-1] += ht_ht*f(sn, tn20);
            }
            dy[0]   -= m_aa_htht__hyhy * u20[0][n];
            dy[M-2] -= m_aa_htht__hyhy * u20[M][n];
            tomasAlgorithm(ay, by, cy, dy, ry, M-1);
            for (unsigned int m=1; m<=M-1; m++) u20[m][n] = ry[m-1];
        }
        /**************************************************** y direction apprx ***************************************************/
        for (unsigned int m=0; m<=M; m++)
        {
            for (unsigned int n=0; n<=N; n++)
            {
                u00[m][n] = u10[m][n];
                u10[m][n] = u20[m][n];
            }
        }
        layerInfo(u20, ln+1);
    }

    u00.clear();
    u10.clear();
    u15.clear();
    u20.clear();

    free(ry);
    free(dy);
    free(cy);
    free(by);
    free(ay);

    free(rx);
    free(dx);
    free(cx);
    free(bx);
    free(ax);
}

void CC1IHyperbolicIBVP::implicit_calculate_D2V2(DoubleMatrix &, double, double, double) const {}

void HyperbolicIBVP::gridMethod(DoubleMatrix &u, SweepMethodDirection direction)
{
    C_UNUSED(u);
    C_UNUSED(direction);
}

void HyperbolicIBVP::gridMethod0(DoubleMatrix &u, SweepMethodDirection direction)
{
    typedef void (*t_algorithm)(const double*, const double*, const double*, const double*, double*, unsigned int);
    t_algorithm algorithm = &tomasAlgorithm;
    if (direction == ForwardSweep) algorithm = &tomasAlgorithmL2R;
    if (direction == BackwardSweep) algorithm = &tomasAlgorithmR2L;

    Dimension time = mtimeDimension;
    Dimension dim1 = mspaceDimension.at(0);

    double ht = time.step();
    unsigned int minM = time.min();
    unsigned int maxM = time.max();
    unsigned int M = maxM-minM;

    double hx = dim1.step();
    unsigned int minN = dim1.min();
    unsigned int maxN = dim1.max();
    unsigned int N = maxN-minN;

    double lambda = 0.25;
    double ht2 = ht*ht;
    double hx2 = hx*hx;
    double h = ht2/hx2;

    u.clear();
    u.resize(M+1, N+1);

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

        u[0][n] = initial1(isn);
        u[1][n] = u[0][n] + ht*initial2(isn);

        //u[2][n] = u[0][n] + 2.0*ht*initial2(isn);

        //TimeNode tn; tn.i = 0.0; tn.t = tn.i*ht;
        //u[2][n] = u[0][n] + ht*initial2(isn) + (2.0*a(isn,tn)+f(isn,tn))*((ht*ht)/2.0);
    }
    IPrinter::printSeperatorLine();
    IPrinter::printVector(14,10,u.row(0));
    IPrinter::printVector(14,10,u.row(1));
    IPrinter::printVector(14,10,u.row(2));
    IPrinter::printSeperatorLine();

    SpaceNodePDE lsn;
    lsn.i = minN;
    lsn.x = minN*hx;

    SpaceNodePDE rsn;
    rsn.i = maxN;
    rsn.x = maxN*hx;

    TimeNodePDE tn;
    TimeNodePDE tn1;
    TimeNodePDE tn2;
    for (unsigned int m=1; m<=M-1; m++)
    {
        tn.i = m+minM+1;
        tn.t = tn.i*ht;

        tn1.i = m+minM;
        tn1.t = tn1.i*ht;

        tn2.i = m+minM-1;
        tn2.t = tn2.i*ht;

        for (unsigned int n=1; n<=N-1; n++)
        {
            isn.i = n+minN;
            isn.x = isn.i*hx;

            //double a0 = a(isn,tn);
            //double a1 = a(isn,tn1);
            //double a2 = a(isn,tn2);

            double alpha = -a(isn,tn)*lambda*h;
            double betta = 1.0 - 2.0*alpha;

            double alpha1 = (1.0 - 2.0*lambda)*a(isn,tn1)*h;
            double alpha2 = lambda*a(isn,tn2)*h;

            ka[n-1] = alpha;
            kb[n-1] = betta;
            kc[n-1] = alpha;
            kd[n-1] = (alpha1*(u[m][n+1]   - 2.0*u[m][n]   + u[m][n-1]))   + 2.0*u[m][n]
                    + (alpha2*(u[m-1][n+1] - 2.0*u[m-1][n] + u[m-1][n-1])) - u[m-1][n]
                    + (ht*ht)*f(isn, tn1);
        }

        ka[0]   = 0.0;
        kc[N-2] = 0.0;

        /* border conditions */
        u[m+1][0] = boundary(lsn, tn);
        u[m+1][N] = boundary(rsn, tn);

        kd[0]   += a(lsn,tn) * lambda * h * u[m+1][0];
        kd[N-2] += a(rsn,tn) * lambda * h * u[m+1][N];

        (*algorithm)(ka, kb, kc, kd, rx, N-1);

        for (unsigned int n=1; n<=N-1; n++) u[m+1][n] = rx[n-1];
    }

    free(ka);
    free(kb);
    free(kc);
    free(kd);
    free(rx);
}

void HyperbolicIBVP::gridMethod1(DoubleMatrix &u, SweepMethodDirection direction)
{
    typedef void (*t_algorithm)(const double*, const double*, const double*, const double*, double*, unsigned int);
    t_algorithm algorithm = &tomasAlgorithm;
    if (direction == ForwardSweep) algorithm = &tomasAlgorithmL2R;
    if (direction == BackwardSweep) algorithm = &tomasAlgorithmR2L;

    Dimension time = mtimeDimension;
    Dimension dim1 = mspaceDimension.at(0);

    double ht = time.step();
    unsigned int minM = time.min();
    unsigned int maxM = time.max();
    unsigned int M = maxM-minM;

    double hx = dim1.step();
    unsigned int minN = dim1.min();
    unsigned int maxN = dim1.max();
    unsigned int N = maxN-minN;

    double lambda = 0.25;
    double h = (ht*ht)/(hx*hx);

    u.clear();
    u.resize(M+1, N+1);

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

        u[0][n] = initial1(isn);
        u[1][n] = u[0][n] + ht*initial2(isn);
    }

    SpaceNodePDE lsn;
    lsn.i = minN;
    lsn.x = minN*hx;

    SpaceNodePDE rsn;
    rsn.i = maxN;
    rsn.x = maxN*hx;

    TimeNodePDE tn;
    TimeNodePDE tn1;
    TimeNodePDE tn2;
    for (unsigned int m=2; m<=M; m++)
    {
        tn.i = m+minM;
        tn.t = tn.i*ht;

        tn1.i = m+minM-1;
        tn1.t = tn1.i*ht;

        tn2.i = m+minM-2;
        tn2.t = tn2.i*ht;

        for (unsigned int n=1; n<=N-1; n++)
        {
            isn.i = n+minN;
            isn.x = isn.i*hx;

            double a0 = a(isn,tn);
            double a1 = a(isn,tn1);
            double a2 = a(isn,tn2);

            double alpha = -a0*lambda*h;
            double betta = 1.0 - 2.0*alpha;

            double alpha1 = (1.0 - 2.0*lambda)*a1*h;
            double alpha2 = lambda*a2*h;

            ka[n-1] = alpha;
            kb[n-1] = betta;
            kc[n-1] = alpha;
            kd[n-1] = (alpha1*(u[m-1][n+1] - 2.0*u[m-1][n] + u[m-1][n-1])) + 2.0*u[m-1][n]
                    + (alpha2*(u[m-2][n+1] - 2.0*u[m-2][n] + u[m-2][n-1])) - u[m-2][n]
                    + (ht*ht)*f(isn, tn);
        }

        ka[0]   = 0.0;
        kc[N-2] = 0.0;

        /* border conditions */
        u[m][0] = boundary(lsn, tn);
        u[m][N] = boundary(rsn, tn);

        kd[0]   += a(lsn,tn) * lambda * h * u[m][0];
        kd[N-2] += a(rsn,tn) * lambda * h * u[m][N];

        (*algorithm)(ka, kb, kc, kd, rx, N-1);

        for (unsigned int n=1; n<=N-1; n++) u[m][n] = rx[n-1];
    }

    free(ka);
    free(kb);
    free(kc);
    free(kd);
    free(rx);
}

void HyperbolicIBVP::gridMethod2(DoubleMatrix &u, SweepMethodDirection direction)
{
    typedef void (*t_algorithm)(const double*, const double*, const double*, const double*, double*, unsigned int);
    t_algorithm algorithm = &tomasAlgorithm;
    if (direction == ForwardSweep) algorithm = &tomasAlgorithmL2R;
    if (direction == BackwardSweep) algorithm = &tomasAlgorithmR2L;

    Dimension time = mtimeDimension;
    Dimension dim1 = mspaceDimension.at(0);

    double ht = time.step();
    unsigned int minM = time.min();
    unsigned int maxM = time.max();
    unsigned int M = maxM-minM;

    double hx = dim1.step();
    unsigned int minN = dim1.min();
    unsigned int maxN = dim1.max();
    unsigned int N = maxN-minN;

    double lambda = 0.25;
    double h = (ht*ht)/(hx*hx);

    u.clear();
    u.resize(M+1, N+1);

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

        u[0][n] = initial1(isn);
        u[1][n] = u[0][n] + ht*initial2(isn);
    }

    SpaceNodePDE lsn;
    lsn.i = minN;
    lsn.x = minN*hx;

    SpaceNodePDE rsn;
    rsn.i = maxN;
    rsn.x = maxN*hx;

    TimeNodePDE tn;
    TimeNodePDE tn1;
    TimeNodePDE tn2;
    for (unsigned int m=2; m<=M; m++)
    {
        tn.i = m+minM;
        tn.t = tn.i*ht;

        tn1.i = m+minM-1;
        tn1.t = tn1.i*ht;

        tn2.i = m+minM-2;
        tn2.t = tn2.i*ht;

        for (unsigned int n=1; n<=N-1; n++)
        {
            isn.i = n+minN;
            isn.x = isn.i*hx;

            double a0 = a(isn,tn);
            double a1 = a(isn,tn1);
            double a2 = a(isn,tn2);

            double alpha = -a0*(1.0-2.0*lambda)*h;
            double betta = 1.0 - 2.0*alpha;

            double alpha1 = lambda*a1*h;
            double alpha2 = lambda*a2*h;

            ka[n-1] = alpha;
            kb[n-1] = betta;
            kc[n-1] = alpha;
            kd[n-1] = (alpha1*(u[m-1][n+1] - 2.0*u[m-1][n] + u[m-1][n-1])) + 2.0*u[m-1][n]
                    + (alpha2*(u[m-2][n+1] - 2.0*u[m-2][n] + u[m-2][n-1])) - u[m-2][n]
                    + (ht*ht)*f(isn, tn);
        }

        ka[0]   = 0.0;
        kc[N-2] = 0.0;

        /* border conditions */
        u[m][0] = boundary(lsn, tn);
        u[m][N] = boundary(rsn, tn);

        kd[0]   += a(lsn,tn) * (1.0-2.0*lambda) * h * u[m][0];
        kd[N-2] += a(rsn,tn) * (1.0-2.0*lambda) * h * u[m][N];

        (*algorithm)(ka, kb, kc, kd, rx, N-1);

        for (unsigned int n=1; n<=N-1; n++) u[m][n] = rx[n-1];
    }

    free(ka);
    free(kb);
    free(kc);
    free(kd);
    free(rx);
}

HyperbolicIBVP::~HyperbolicIBVP() {}

void ConjugateCC1IHyperbolicIBVP::explicit_calculate_D2V1(DoubleMatrix &u, double a, double alpha) const
{
    const Dimension &dimX = spaceDimension(Dimension::DimensionX);
    const Dimension &dimY = spaceDimension(Dimension::DimensionX);
    const Dimension &time = timeDimension();

    const unsigned int N = static_cast<unsigned int>( dimX.size() );
    const unsigned int M = static_cast<unsigned int>( dimY.size() );
    const unsigned int L = static_cast<unsigned int>( time.size() );

    const double hx = dimX.step();
    const double hy = dimY.step();
    const double ht = time.step();

    const double alpha_ht_05 = alpha*ht*0.5;
    const double inv__alpha_ht = 1.0/(1.0 + alpha_ht_05);
    const double aa_htht__hxhx = ((a*a*ht*ht)/(hx*hx));
    const double aa_htht__hyhy = ((a*a*ht*ht)/(hy*hy));

    DoubleMatrix u00(M+1, N+1);
    DoubleMatrix u10(M+1, N+1);
    DoubleMatrix u20(M+1, N+1);

    /**********************************************************************/
    initial_calculate(u00, u10, N, hx, M, hy, ht, a, alpha);
    /**********************************************************************/

    for (unsigned int ln=L-2; ln!=(unsigned int)0-1; ln++)
    {
        TimeNodePDE tn20; tn20.i = ln-0; tn20.t = tn20.i*ht;
        TimeNodePDE tn10; tn10.i = ln+1; tn10.t = tn10.i*ht;

        /**************************** border conditions ****************************/
        SpaceNodePDE sn0, sn1;
        sn0.i = static_cast<int>(0); sn0.x = 0*hx;
        sn1.i = static_cast<int>(N); sn1.x = N*hx;
        for (unsigned int m=0; m<=M; m++)
        {
            sn0.j = static_cast<int>(m); sn0.y = m*hy; u20[m][0] = boundary(sn0, tn20);
            sn1.j = static_cast<int>(m); sn1.y = m*hy; u20[m][N] = boundary(sn1, tn20);
        }

        sn0.j = static_cast<int>(0); sn0.y = 0*hy;
        sn1.j = static_cast<int>(M); sn1.y = M*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn0.i = static_cast<int>(n); sn0.x = n*hx; u20[0][n] = boundary(sn0, tn20);
            sn1.i = static_cast<int>(n); sn1.x = n*hx; u20[M][n] = boundary(sn1, tn20);
        }
        /**************************** border conditions ****************************/

        SpaceNodePDE sn;
        for (unsigned int m=1; m<=M-1; m++)
        {
            sn.j = static_cast<int>(m); sn.y = m*hy;
            for (unsigned int n=1; n<=N-1; n++)
            {
                sn.i = static_cast<int>(n); sn.x = n*hx;
                u20[m][n] = inv__alpha_ht * (aa_htht__hxhx*(u10[m][n-1]-2.0*u10[m][n]+u10[m][n+1])
                                           + aa_htht__hyhy*(u10[m-1][n]-2.0*u10[m][n]+u10[m+1][n])
                                           + alpha_ht_05*u00[m][n] + 2.0*u10[m][n] - u00[m][n] + ht*ht*f(sn,tn10));
            }
        }
        layerInfo(u20, ln);

        for (unsigned int m=0; m<=M; m++)
        {
            for (unsigned int n=0; n<=N; n++)
            {
                u00[m][n] = u10[m][n];
                u10[m][n] = u20[m][n];
            }
        }
    }

    for (unsigned int m=0; m<=M; m++)
    {
        for (unsigned int n=0; n<=N; n++)
        {
            u[m][n] = u20[m][n];
        }
    }
}

void ConjugateCC1IHyperbolicIBVP::implicit_calculate_D2V1(DoubleMatrix &, double a, double alpha) const
{
    const Dimension &dimX = spaceDimension(Dimension::DimensionX);
    const Dimension &dimY = spaceDimension(Dimension::DimensionX);
    const Dimension &time = timeDimension();

    const unsigned int N = static_cast<unsigned int>( dimX.size() );
    const unsigned int M = static_cast<unsigned int>( dimY.size() );
    const unsigned int L = static_cast<unsigned int>( time.size() );

    const double hx = dimX.step();
    const double hy = dimY.step();
    const double ht = time.step();

    const double m_aa_htht__hxhx = -(a*a)*((ht*ht)/(hx*hx));
    const double b_aa_htht__hxhx = +(2.0 + 2.0*(a*a)*((ht*ht)/(hx*hx)) + alpha*ht);
    const double p_aa_htht__hyhy = +(a*a)*((ht*ht)/(hy*hy));

    const double m_aa_htht__hyhy = -(a*a)*((ht*ht)/(hy*hy));
    const double b_aa_htht__hyhy = +(2.0 + 2.0*(a*a)*((ht*ht)/(hy*hy)) + alpha*ht);
    const double p_aa_htht__hxhx = +(a*a)*((ht*ht)/(hx*hx));

    const double ht_ht = ht*ht;
    const double alpha_ht = alpha*ht;

    double *ax = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *bx = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *cx = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *dx = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *rx = static_cast<double*>(malloc(sizeof(double)*(N-1)));

    for (unsigned int n=0; n<=N-2; n++)
    {
        ax[n] = m_aa_htht__hxhx;
        bx[n] = b_aa_htht__hxhx;
        cx[n] = m_aa_htht__hxhx;
    }
    ax[0] = 0.0; cx[N-2] = 0.0;

    double *ay = static_cast<double*>(malloc(sizeof(double)*(M-1)));
    double *by = static_cast<double*>(malloc(sizeof(double)*(M-1)));
    double *cy = static_cast<double*>(malloc(sizeof(double)*(M-1)));
    double *dy = static_cast<double*>(malloc(sizeof(double)*(M-1)));
    double *ry = static_cast<double*>(malloc(sizeof(double)*(M-1)));

    for (unsigned int m=0; m<=M-2; m++)
    {
        ay[m] = m_aa_htht__hyhy;
        by[m] = b_aa_htht__hyhy;
        cy[m] = m_aa_htht__hyhy;
    }
    ay[0] = 0.0; cy[M-2] = 0.0;

    DoubleMatrix u00(M+1, N+1);
    DoubleMatrix u10(M+1, N+1);
    DoubleMatrix u15(M+1, N+1);
    DoubleMatrix u20(M+1, N+1);

    /***********************************************************************************************/
    //initial_calculate(u00, u10, N, hx, M, hy, ht, a, alpha, tn);
    /***********************************************************************************************/

    SpaceNodePDE sn;
    for (unsigned int ln=1; ln<=L-1; ln++)
    {
        TimeNodePDE tn10; tn10.i = ln; tn10.t = ln*ht;
        TimeNodePDE tn15; tn15.i = ln; tn15.t = ln*ht+0.5*ht;
        TimeNodePDE tn20; tn20.i = ln; tn20.t = ln*ht+1.0*ht;

        /**************************************************** border conditions ***************************************************/

        SpaceNodePDE sn0;
        SpaceNodePDE sn1;

        sn0.i = 0; sn0.x = 0.0; sn1.i = static_cast<int>(N); sn1.x = hx*N;
        for (unsigned int m=0; m<=M; m++)
        {
            sn0.j = static_cast<int>(m); sn0.y = m*hy; u15[m][0] = boundary(sn0, tn15); u20[m][0] = boundary(sn0, tn20);
            sn1.j = static_cast<int>(m); sn1.y = m*hy; u15[m][N] = boundary(sn1, tn15); u20[m][N] = boundary(sn1, tn20);
        }

        sn0.j = 0; sn0.y = 0.0; sn1.j = static_cast<int>(M); sn1.y = hy*M;
        for (unsigned int n=0; n<=N; n++)
        {
            sn0.i = static_cast<int>(n); sn0.x = n*hx; u15[0][n] = boundary(sn0, tn15); u20[0][n] = boundary(sn0, tn20);
            sn1.i = static_cast<int>(n); sn1.x = n*hx; u15[M][n] = boundary(sn1, tn15); u20[M][n] = boundary(sn1, tn20);
        }

        /**************************************************** border conditions ***************************************************/

        /**************************************************** x direction apprx ***************************************************/

        for (unsigned int m=1; m<=M-1; m++)
        {
            sn.j = static_cast<int>(m); sn.y = m*hy;

            for (unsigned int n=1; n<=N-1; n++)
            {
                sn.i = static_cast<int>(n); sn.x = n*hx;

                dx[n-1] = 0.0;
                dx[n-1] += p_aa_htht__hyhy*(u10[m-1][n] - 2.0*u10[m][n] + u10[m+1][n]);
                dx[n-1] += 2.0*u10[m][n] + (u10[m][n] - u00[m][n]);
                dx[n-1] += alpha_ht*(u10[m][n] - 0.5*(u10[m][n]-u00[m][n]));
                dx[n-1] += ht_ht*f(sn, tn10);
            }
            dx[0]   -= m_aa_htht__hxhx * u15[m][0];
            dx[N-2] -= m_aa_htht__hxhx * u15[m][N];
            tomasAlgorithm(ax, bx, cx, dx, rx, N-1);
            for (unsigned int n=1; n<=N-1; n++) u15[m][n] = rx[n-1];
        }

        /**************************************************** x direction apprx ***************************************************/

        /**************************************************** y direction apprx ***************************************************/

        for (unsigned int n=1; n<=N-1; n++)
        {
            sn.i = static_cast<int>(n); sn.x = n*hx;

            for (unsigned int m=1; m<=M-1; m++)
            {
                sn.j = static_cast<int>(m); sn.y = m*hy;

                dy[m-1] = 0.0;
                dy[m-1] += p_aa_htht__hxhx*(u15[m][n-1] - 2.0*u15[m][n] + u15[m][n+1]);
                dy[m-1] += 2.0*u15[m][n] + (u10[m][n]-u00[m][n]);
                dy[m-1] += alpha_ht*(u15[m][n] - 0.5*(u10[m][n]-u00[m][n]));
                dy[m-1] += ht_ht*f(sn, tn10);
            }
            dy[0]   -= m_aa_htht__hyhy * u20[0][n];
            dy[M-2] -= m_aa_htht__hyhy * u20[M][n];
            tomasAlgorithm(ay, by, cy, dy, ry, M-1);
            for (unsigned int m=1; m<=M-1; m++) u20[m][n] = ry[m-1];
        }

        /**************************************************** y direction apprx ***************************************************/

        for (unsigned int m=0; m<=M; m++)
        {
            for (unsigned int n=0; n<=N; n++)
            {
                u00[m][n] = u10[m][n];
                u10[m][n] = u20[m][n];
            }
        }

        layerInfo(u20, ln+1);
    }

    u00.clear();
    u10.clear();
    u15.clear();
    u20.clear();

    free(ry);
    free(dy);
    free(cy);
    free(by);
    free(ay);

    free(rx);
    free(dx);
    free(cx);
    free(bx);
    free(ax);
}

void ConjugateCC1IHyperbolicIBVP::initial_calculate(DoubleMatrix &u00, DoubleMatrix &u10, unsigned int N, double hx, unsigned int M, double hy, double ht, double a, double alpha) const
{
    //const Dimension &dimX = spaceDimension(Dimension::DimensionX);
    //const Dimension &dimY = spaceDimension(Dimension::DimensionX);
    const Dimension &time = timeDimension();

    //const unsigned int N = static_cast<unsigned int>( dimX.size() );
    //const unsigned int M = static_cast<unsigned int>( dimY.size() );
    const unsigned int L = static_cast<unsigned int>( time.size() );

    //const double hx = dimX.step();
    //const double hy = dimY.step();
    //const double ht = time.step();

    const double aa__hxhx = (a*a)/(hx*hx);
    const double aa__hyhy = (a*a)/(hy*hy);
    const double htht_05 = 0.5*ht*ht;

    /***********************************************************************************************/

    SpaceNodePDE sn;
    for (unsigned int m=0; m<=M; m++)
    {
        sn.j = static_cast<int>(m); sn.y = sn.j*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = static_cast<int>(n); sn.x = sn.i*hx;
            u00[m][n] = initial1(sn);
        }
    }
    layerInfo(u00, L-0);

    /***********************************************************************************************/

    TimeNodePDE tn00; tn00.i = L-0; tn00.t = tn00.i*ht;
    TimeNodePDE tn10; tn10.i = L-1; tn10.t = tn10.i*ht;
    SpaceNodePDE sn0;
    SpaceNodePDE sn1;
    sn0.i = static_cast<int>(0); sn0.x = 0.00;
    sn1.i = static_cast<int>(N); sn1.x = N*hx;
    for (unsigned int m=0; m<=M; m++)
    {
        sn0.j = static_cast<int>(m); sn0.y = m*hy; u10[m][0] = boundary(sn0, tn10);
        sn1.j = static_cast<int>(m); sn1.y = m*hy; u10[m][N] = boundary(sn1, tn10);
    }
    sn0.j = static_cast<int>(0); sn0.y = 0.00;
    sn1.j = static_cast<int>(M); sn1.y = M*hy;
    for (unsigned int n=0; n<=N; n++)
    {
        sn0.i = static_cast<int>(n); sn0.x = n*hx; u10[0][n] = boundary(sn0, tn10);
        sn1.i = static_cast<int>(n); sn1.x = n*hx; u10[M][n] = boundary(sn1, tn10);
    }

    for (unsigned int m=1; m<=M-1; m++)
    {
        sn.j = static_cast<int>(m); sn.y = m*hy;
        for (unsigned int n=1; n<=N-1; n++)
        {
            sn.i = static_cast<int>(n); sn.x = n*hx;
            u10[m][n] = u00[m][n] - ht*initial2(sn);
            u10[m][n] += htht_05*(aa__hxhx*(u00[m][n-1]-2.0*u00[m][n]+u00[m][n+1])+aa__hyhy*(u00[m-1][n]-2.0*u00[m][n]+u00[m+1][n])+alpha*initial2(sn)+f(sn,tn00));
        }
    }
    layerInfo(u10, L-1);

    /***********************************************************************************************/
}
