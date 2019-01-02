#include "hibvp.h"

IHyperbolicIBVP::~IHyperbolicIBVP() {}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

CcIHyperbolicIBVP::~CcIHyperbolicIBVP() {}

void CcIHyperbolicIBVP::explicit_calculate_D1V1(DoubleVector &u, double a) const
{
    const Dimension &dimX = spaceDimension(Dimension::DimensionX);
    const Dimension &time = timeDimension();

    const unsigned int N = static_cast<unsigned int>( dimX.size() );
    const unsigned int M = static_cast<unsigned int>( time.size() );

    const double hx = dimX.step();
    const double ht = time.step();

    const double aa_htht__hxhx = ((a*a*ht*ht)/(hx*hx));
    const double ht_ht = ht*ht;

    DoubleVector u00(N+1);
    DoubleVector u10(N+1);
    DoubleVector u20(N+1);

    /**********************************************************************/
    initial_calculate(u00, u10, N, hx, ht, a);
    /**********************************************************************/

    SpaceNodePDE sn;
    for (unsigned int m=2; m<=M; m++)
    {
        TimeNodePDE tn20; tn20.i = m;   tn20.t = tn20.i*ht;
        TimeNodePDE tn10; tn10.i = m-1; tn10.t = tn10.i*ht;

        /**************************** border conditions ****************************/
        sn.i = static_cast<int>(0); sn.x = sn.i*hx; u20[0] = boundary(sn, tn20);
        sn.i = static_cast<int>(N); sn.x = sn.i*hx; u20[N] = boundary(sn, tn20);
        /**************************** border conditions ****************************/
        for (unsigned int n=1; n<=N-1; n++)
        {
            sn.i = static_cast<int>(n); sn.x = n*hx;
            u20[n] = aa_htht__hxhx*(u10[n-1]-2.0*u10[n]+u10[n+1]) + (2.0*u10[n] - u00[n]) + ht_ht*f(sn,tn10);
        }

        layerInfo(u20, m);

        for (unsigned int n=0; n<=N; n++)
        {
            u00[n] = u10[n];
            u10[n] = u20[n];
        }
    }

    u.clear(); u.resize(N+1); for (unsigned int n=0; n<=N; n++) { u[n] = u20[n]; }

    u20.clear();
    u10.clear();
    u00.clear();
}

void CcIHyperbolicIBVP::implicit_calculate_D1V1(DoubleVector &u, double a) const
{
    const Dimension &dimx = spaceDimension(Dimension::DimensionX);
    const Dimension &time = timeDimension();

    const unsigned int N = static_cast<const unsigned int>(dimx.size());
    const unsigned int M = static_cast<const unsigned int>(time.size());

    const double hx = dimx.step();
    const double ht = time.step();

    const double alpha = -(a*a)*((ht*ht)/(hx*hx));
    const double betta = +(1.0 + 2.0*(a*a)*((ht*ht)/(hx*hx)));
    //
    //
    const double ht_ht = ht*ht;

    DoubleVector u00(N+1);
    DoubleVector u10(N+1);
    DoubleVector u20(N+1);

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

    /**********************************************************************/
    initial_calculate(u00, u10, N, hx, ht, a);
    /**********************************************************************/

    SpaceNodePDE sn;
    TimeNodePDE tn;
    for (unsigned int m=2; m<=M; m++)
    {
        tn.i = m; tn.t = tn.i*ht;

        /**************************** border conditions ****************************/
        sn.i = static_cast<int>(0); sn.x = sn.i*hx; u20[0] = boundary(sn, tn);
        sn.i = static_cast<int>(N); sn.x = sn.i*hx; u20[N] = boundary(sn, tn);
        /**************************** border conditions ****************************/

        for (unsigned int n=1; n<=N-1; n++)
        {
            sn.i = static_cast<int>(n); sn.x = n*hx;
            dd[n-1] = 2.0*u10[n] - u00[n] + ht_ht*f(sn, tn);
        }
        dd[0]   -= alpha*u20[0];
        dd[N-2] -= alpha*u20[N];

        tomasAlgorithm(da, db, dc, dd, rx, N-1);
        for (unsigned int n=1; n<=N-1; n++) u20[n] = rx[n-1];

        layerInfo(u20, m);

        for (unsigned int n=0; n<=N; n++)
        {
            u00[n] = u10[n];
            u10[n] = u20[n];
        }
    }

    u.clear(); u.resize(N+1); for (unsigned int n=0; n<=N; n++) { u[n] = u20[n]; }

    u20.clear();
    u10.clear();
    u00.clear();
}

void CcIHyperbolicIBVP::implicit_calculate_D1V2(DoubleVector &u, double a, double lambda) const
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

    DoubleVector u00(N+1);
    DoubleVector u10(N+1);
    DoubleVector u20(N+1);

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

    /**********************************************************************/
    initial_calculate(u00, u10, N, hx, ht, a);
    /**********************************************************************/

    SpaceNodePDE sn;
    TimeNodePDE tn00, tn10, tn20;
    for (unsigned int m=2; m<=M; m++)
    {
        tn20.i = m+0; tn20.t = tn20.i*ht;
        tn10.i = m-1; tn10.t = tn10.i*ht;
        tn00.i = m-2; tn00.t = tn00.i*ht;

        /**************************** border conditions ****************************/
        sn.i = static_cast<int>(0); sn.x = sn.i*hx; u20[0] = boundary(sn, tn20);
        sn.i = static_cast<int>(N); sn.x = sn.i*hx; u20[N] = boundary(sn, tn20);
        /**************************** border conditions ****************************/

        for (unsigned int n=1; n<=N-1; n++)
        {
            sn.i = static_cast<int>(n); sn.x = n*hx;
            dd[n-1] = gamma*(u10[n-1]-2.0*u10[n]+u10[n+1])
                    + theta*(u00[n-1]-2.0*u00[n]+u00[n+1])
                    + 2.0*u10[n] - u00[n];
            //dd[n-1] += ht_ht*(lambda*f(sn, tn20) + (1.0-2.0*lambda)*f(sn, tn10) + lambda*f(sn, tn00));
            dd[n-1] += ht_ht*f(sn, tn10);
        }
        dd[0]   -= alpha*u20[0];
        dd[N-2] -= alpha*u20[N];

        tomasAlgorithm(da, db, dc, dd, rx, N-1);
        for (unsigned int n=1; n<=N-1; n++) u20[n] = rx[n-1];

        layerInfo(u20, m);

        for (unsigned int n=0; n<=N; n++)
        {
            u00[n] = u10[n];
            u10[n] = u20[n];
        }
    }

    u.clear(); u.resize(N+1); for (unsigned int n=0; n<=N; n++) { u[n] = u20[n]; }

    free(rx);
    free(dd);
    free(dc);
    free(db);
    free(da);

    u20.clear();
    u10.clear();
    u00.clear();
}

void CcIHyperbolicIBVP::explicit_calculate_D2V1(DoubleMatrix &u, double a) const
{
    const Dimension &dimX = spaceDimension(Dimension::DimensionX);
    const Dimension &dimY = spaceDimension(Dimension::DimensionY);
    const Dimension &time = timeDimension();

    const unsigned int N = static_cast<unsigned int>( dimX.size() );
    const unsigned int M = static_cast<unsigned int>( dimY.size() );
    const unsigned int L = static_cast<unsigned int>( time.size() );

    const double hx = dimX.step();
    const double hy = dimY.step();
    const double ht = time.step();

    const double aa_htht__hxhx = ((a*a*ht*ht)/(hx*hx));
    const double aa_htht__hyhy = ((a*a*ht*ht)/(hy*hy));

    DoubleMatrix u00(M+1, N+1);
    DoubleMatrix u10(M+1, N+1);
    DoubleMatrix u20(M+1, N+1);

    /**********************************************************************/
    initial_calculate(u00, u10, N, hx, M, hy, ht, a);
    /**********************************************************************/

    for (unsigned int ln=2; ln<=L; ln++)
    {
        TimeNodePDE tn20; tn20.i = ln;   tn20.t = tn20.i*ht;
        TimeNodePDE tn10; tn10.i = ln-1; tn10.t = tn10.i*ht;

        /**************************** border conditions ****************************/
        border2_calculate(u20, N, hx, M, hy, tn20);
        /**************************** border conditions ****************************/

        SpaceNodePDE sn;
        for (unsigned int m=1; m<=M-1; m++)
        {
            sn.j = static_cast<int>(m); sn.y = m*hy;
            for (unsigned int n=1; n<=N-1; n++)
            {
                sn.i = static_cast<int>(n); sn.x = n*hx;
                u20[m][n] = (aa_htht__hxhx*(u10[m][n-1]-2.0*u10[m][n]+u10[m][n+1])
                        + aa_htht__hyhy*(u10[m-1][n]-2.0*u10[m][n]+u10[m+1][n])
                        + 2.0*u10[m][n] - u00[m][n] + ht*ht*f(sn,tn10));
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

    u.clear(); u.resize(M+1, N+1); for (unsigned int m=0; m<=M; m++) { for (unsigned int n=0; n<=N; n++) { u[m][n] = u20[m][n]; } }
}

void CcIHyperbolicIBVP::implicit_calculate_D2V1(DoubleMatrix &u, double a) const
{
    const Dimension &dimX = spaceDimension(Dimension::DimensionX);
    const Dimension &dimY = spaceDimension(Dimension::DimensionY);
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
    //

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
    initial_calculate(u00, u10, N, hx, M, hy, ht, a);
    /***********************************************************************************************/

    SpaceNodePDE sn;
    for (unsigned int ln=1; ln<=L-1; ln++)
    {
        TimeNodePDE tn10; tn10.i = ln;   tn10.t = tn10.i*ht;
        TimeNodePDE tn15; tn15.i = ln;   tn15.t = tn15.i*ht+0.5*ht;
        TimeNodePDE tn20; tn20.i = ln+1; tn20.t = tn20.i*ht;
        /**************************************************** border conditions ***************************************************/
        border1_calculate(u15, u20, N, hx, M, hy, tn15, tn20);
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
                dx[n-1] += 2.0*u10[m][n] + (u10[m][n]-u00[m][n]);
                //
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
                //
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
    u.clear(); u.resize(M+1, N+1); for (unsigned int m=0; m<=M; m++) { for (unsigned int n=0; n<=N; n++) { u[m][n] = u20[m][n]; } }

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

void CcIHyperbolicIBVP::implicit_calculate_D2V2(DoubleMatrix &, double, double) const {}

void CcIHyperbolicIBVP::implicit_calculate_D2V3(DoubleMatrix &u, double a) const
{
    const Dimension &dimX = spaceDimension(Dimension::DimensionX);
    const Dimension &dimY = spaceDimension(Dimension::DimensionY);
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
    //

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
    initial_calculate(u00, u10, N, hx, M, hy, ht, a);
    /***********************************************************************************************/

    SpaceNodePDE sn;
    for (unsigned int ln=1; ln<=L-1; ln++)
    {
        TimeNodePDE tn10; tn10.i = ln;   tn10.t = tn10.i*ht;
        TimeNodePDE tn15; tn15.i = ln;   tn15.t = tn15.i*ht+0.5*ht;
        TimeNodePDE tn20; tn20.i = ln+1; tn20.t = tn20.i*ht;
        /**************************************************** border conditions ***************************************************/
        border1_calculate(u15, u20, N, hx, M, hy, tn15, tn20);
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
                //
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
                //
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
    u.clear(); u.resize(M+1, N+1); for (unsigned int m=0; m<=M; m++) { for (unsigned int n=0; n<=N; n++) { u[m][n] = u20[m][n]; } }

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

void CcIHyperbolicIBVP::initial_calculate(DoubleVector &u00, DoubleVector &u10, unsigned int N, double hx, double ht, double a) const
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
    }
    layerInfo(u10, 1);
    /***********************************************************************************************/
}

void CcIHyperbolicIBVP::initial_calculate(DoubleMatrix &u00, DoubleMatrix &u10, unsigned int N, double hx, unsigned int M, double hy, double ht, double a) const
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
            u10[m][n] += htht_05*(aa__hxhx*(u00[m][n-1]-2.0*u00[m][n]+u00[m][n+1])+aa__hyhy*(u00[m-1][n]-2.0*u00[m][n]+u00[m+1][n])+f(sn,tn00));
        }
    }
    layerInfo(u10, 1);

    /***********************************************************************************************/
}

void CcIHyperbolicIBVP::border1_calculate(DoubleMatrix &u15, DoubleMatrix &u20, unsigned int N, double hx, unsigned int M, double hy, const TimeNodePDE &tn15, const TimeNodePDE &tn20) const
{
    SpaceNodePDE sn0, sn1;

    sn0.i = static_cast<int>(0); sn0.x = 0*hx;
    sn1.i = static_cast<int>(N); sn1.x = N*hx;
    for (unsigned int m=0; m<=M; m++)
    {
        sn0.j = static_cast<int>(m); sn0.y = m*hy; u15[m][0] = boundary(sn0, tn15); u20[m][0] = boundary(sn0, tn20);
        sn1.j = static_cast<int>(m); sn1.y = m*hy; u15[m][N] = boundary(sn1, tn15); u20[m][N] = boundary(sn1, tn20);
    }

    sn0.j = static_cast<int>(0); sn0.y = 0*hy;
    sn1.j = static_cast<int>(M); sn1.y = M*hy;
    for (unsigned int n=0; n<=N; n++)
    {
        sn0.i = static_cast<int>(n); sn0.x = n*hx; u15[0][n] = boundary(sn0, tn15); u20[0][n] = boundary(sn0, tn20);
        sn1.i = static_cast<int>(n); sn1.x = n*hx; u15[M][n] = boundary(sn1, tn15); u20[M][n] = boundary(sn1, tn20);
    }
}

void CcIHyperbolicIBVP::border2_calculate(DoubleMatrix &u20, unsigned int N, double hx, unsigned int M, double hy, const TimeNodePDE &tn20) const
{
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
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

CdIHyperbolicIBVP::~CdIHyperbolicIBVP() {}

void CdIHyperbolicIBVP::explicit_calculate_D1V1(DoubleVector &, double, double) const {}

void CdIHyperbolicIBVP::implicit_calculate_D1V1(DoubleVector &, double, double) const {}

void CdIHyperbolicIBVP::implicit_calculate_D1V2(DoubleVector &, double, double, double) const {}

void CdIHyperbolicIBVP::explicit_calculate_D2V1(DoubleMatrix &u, double a, double alpha) const
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
        TimeNodePDE tn20; tn20.i = ln;   tn20.t = tn20.i*ht;
        TimeNodePDE tn10; tn10.i = ln-1; tn10.t = tn10.i*ht;

        /**************************** border conditions ****************************/
        border2_calculate(u20, N, hx, M, hy, tn20);
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

    u.clear(); u.resize(M+1, N+1); for (unsigned int m=0; m<=M; m++) { for (unsigned int n=0; n<=N; n++) { u[m][n] = u20[m][n]; } }
}

void CdIHyperbolicIBVP::implicit_calculate_D2V1(DoubleMatrix &u, double a, double alpha) const
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
        TimeNodePDE tn10; tn10.i = ln;   tn10.t = tn10.i*ht;
        TimeNodePDE tn15; tn15.i = ln;   tn15.t = tn15.i*ht+0.5*ht;
        TimeNodePDE tn20; tn20.i = ln+1; tn20.t = tn20.i*ht;
        /**************************************************** border conditions ***************************************************/
        border1_calculate(u15, u20, N, hx, M, hy, tn15, tn20);
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
    u.clear(); u.resize(M+1, N+1); for (unsigned int m=0; m<=M; m++) { for (unsigned int n=0; n<=N; n++) { u[m][n] = u20[m][n]; } }

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

void CdIHyperbolicIBVP::implicit_calculate_D2V3(DoubleMatrix &u, double a, double alpha) const
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
        TimeNodePDE tn10; tn10.i = ln;   tn10.t = tn10.i*ht;
        TimeNodePDE tn15; tn15.i = ln;   tn15.t = tn15.i*ht+0.5*ht;
        TimeNodePDE tn20; tn20.i = ln+1; tn20.t = tn20.i*ht;
        /**************************************************** border conditions ***************************************************/
        border1_calculate(u15, u20, N, hx, M, hy, tn15, tn20);
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
                dy[m-1] += alpha_ht*(u15[m][n] - 0.5*(u10[m][n]-u00[m][n]));
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
    u.clear(); u.resize(M+1, N+1); for (unsigned int m=0; m<=M; m++) { for (unsigned int n=0; n<=N; n++) { u[m][n] = u20[m][n]; } }

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

void CdIHyperbolicIBVP::implicit_calculate_D2V2(DoubleMatrix &, double, double, double) const {}

void CdIHyperbolicIBVP::initial_calculate(DoubleVector &u00, DoubleVector &u10, unsigned int N, double hx, double ht, double a, double alpha) const
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

void CdIHyperbolicIBVP::initial_calculate(DoubleMatrix &u00, DoubleMatrix &u10, unsigned int N, double hx, unsigned int M, double hy, double ht, double a, double alpha) const
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

void CdIHyperbolicIBVP::border1_calculate(DoubleMatrix &u15, DoubleMatrix &u20, unsigned int N, double hx, unsigned int M, double hy, const TimeNodePDE &tn15, const TimeNodePDE &tn20) const
{
    SpaceNodePDE sn0, sn1;

    sn0.i = static_cast<int>(0); sn0.x = 0*hx;
    sn1.i = static_cast<int>(N); sn1.x = N*hx;
    for (unsigned int m=0; m<=M; m++)
    {
        sn0.j = static_cast<int>(m); sn0.y = m*hy; u15[m][0] = boundary(sn0, tn15); u20[m][0] = boundary(sn0, tn20);
        sn1.j = static_cast<int>(m); sn1.y = m*hy; u15[m][N] = boundary(sn1, tn15); u20[m][N] = boundary(sn1, tn20);
    }

    sn0.j = static_cast<int>(0); sn0.y = 0*hy;
    sn1.j = static_cast<int>(M); sn1.y = M*hy;
    for (unsigned int n=0; n<=N; n++)
    {
        sn0.i = static_cast<int>(n); sn0.x = n*hx; u15[0][n] = boundary(sn0, tn15); u20[0][n] = boundary(sn0, tn20);
        sn1.i = static_cast<int>(n); sn1.x = n*hx; u15[M][n] = boundary(sn1, tn15); u20[M][n] = boundary(sn1, tn20);
    }
}

void CdIHyperbolicIBVP::border2_calculate(DoubleMatrix &u20, unsigned int N, double hx, unsigned int M, double hy, const TimeNodePDE &tn20) const
{
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
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

ConjugateCdIHyperbolicIBVP::~ConjugateCdIHyperbolicIBVP() {}

void ConjugateCdIHyperbolicIBVP::explicit_calculate_D2V1(DoubleMatrix &p, double a, double alpha) const
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

    DoubleMatrix p00(M+1, N+1);
    DoubleMatrix p10(M+1, N+1);
    DoubleMatrix p20(M+1, N+1);

    /**********************************************************************/
    initial_calculate(p00, p10, N, hx, M, hy, ht, a, alpha);
    /**********************************************************************/

    const unsigned int size_ln = static_cast<unsigned int>(0)-1;
    for (unsigned int ln=L-2; ln != size_ln; ln--)
    {
        TimeNodePDE tn20; tn20.i = ln;   tn20.t = tn20.i*ht;
        TimeNodePDE tn10; tn10.i = ln+1; tn10.t = tn10.i*ht;

        /**************************** border conditions ****************************/
        border2_calculate(p20, N, hx, M, hy, tn20);
        /**************************** border conditions ****************************/

        SpaceNodePDE sn;
        for (unsigned int m=1; m<=M-1; m++)
        {
            sn.j = static_cast<int>(m); sn.y = m*hy;
            for (unsigned int n=1; n<=N-1; n++)
            {
                sn.i = static_cast<int>(n); sn.x = n*hx;
                p20[m][n] = inv__alpha_ht * (aa_htht__hxhx*(p10[m][n-1]-2.0*p10[m][n]+p10[m][n+1])
                        + aa_htht__hyhy*(p10[m-1][n]-2.0*p10[m][n]+p10[m+1][n])
                        + alpha_ht_05*p00[m][n] + 2.0*p10[m][n] - p00[m][n] + ht*ht*f(sn,tn10));
            }
        }
        layerInfo(p20, ln);

        for (unsigned int m=0; m<=M; m++)
        {
            for (unsigned int n=0; n<=N; n++)
            {
                p00[m][n] = p10[m][n];
                p10[m][n] = p20[m][n];
            }
        }
    }

    p.clear();
    p.resize(M+1, N+1);
    for (unsigned int m=0; m<=M; m++)
    {
        for (unsigned int n=0; n<=N; n++)
        {
            p[m][n] = p20[m][n];
        }
    }
}

void ConjugateCdIHyperbolicIBVP::implicit_calculate_D2V1(DoubleMatrix &p, double a, double alpha) const
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

    DoubleMatrix p00(M+1, N+1);
    DoubleMatrix p10(M+1, N+1);
    DoubleMatrix p15(M+1, N+1);
    DoubleMatrix p20(M+1, N+1);

    /***********************************************************************************************/
    initial_calculate(p00, p10, N, hx, M, hy, ht, a, alpha);
    /***********************************************************************************************/

    SpaceNodePDE sn;
    const unsigned int size_ln = static_cast<unsigned int>(0)-1;
    for (unsigned int ln=L-1; ln!=size_ln; ln--)
    {
        TimeNodePDE tn10; tn10.i = ln;   tn10.t = tn10.i*ht;
        TimeNodePDE tn15; tn15.i = ln;   tn15.t = tn15.i*ht-0.5*ht;
        TimeNodePDE tn20; tn20.i = ln-1; tn20.t = tn20.i*ht;
        /**************************************************** border conditions ***************************************************/
        border1_calculate(p15, p20, N, hx, M, hy, tn15, tn20);
        /**************************************************** border conditions ***************************************************/
        /**************************************************** x direction apprx ***************************************************/
        for (unsigned int m=1; m<=M-1; m++)
        {
            sn.j = static_cast<int>(m); sn.y = m*hy;
            for (unsigned int n=1; n<=N-1; n++)
            {
                sn.i = static_cast<int>(n); sn.x = n*hx;
                dx[n-1] = 0.0;
                dx[n-1] += p_aa_htht__hyhy*(p10[m-1][n] - 2.0*p10[m][n] + p10[m+1][n]);
                dx[n-1] += 2.0*p10[m][n] + (p10[m][n] - p00[m][n]);
                dx[n-1] += 2.0*alpha_ht*p10[m][n];
                dx[n-1] += ht_ht*f(sn, tn15);
            }
            dx[0]   -= m_aa_htht__hxhx * p15[m][0];
            dx[N-2] -= m_aa_htht__hxhx * p15[m][N];
            tomasAlgorithm(ax, bx, cx, dx, rx, N-1);
            for (unsigned int n=1; n<=N-1; n++) p15[m][n] = rx[n-1];
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
                dy[m-1] += p_aa_htht__hxhx*(p15[m][n-1] - 2.0*p15[m][n] + p15[m][n+1]);
                dy[m-1] += 2.0*p15[m][n] + (p10[m][n]-p00[m][n]);
                dy[m-1] += 2.0*alpha_ht*p15[m][n];
                dy[m-1] += ht_ht*f(sn, tn20);
            }
            dy[0]   -= m_aa_htht__hyhy * p20[0][n];
            dy[M-2] -= m_aa_htht__hyhy * p20[M][n];
            tomasAlgorithm(ay, by, cy, dy, ry, M-1);
            for (unsigned int m=1; m<=M-1; m++) p20[m][n] = ry[m-1];
        }
        /**************************************************** y direction apprx ***************************************************/
        for (unsigned int m=0; m<=M; m++)
        {
            for (unsigned int n=0; n<=N; n++)
            {
                p00[m][n] = p10[m][n];
                p10[m][n] = p20[m][n];
            }
        }
        layerInfo(p20, ln-1);
    }

    p.clear();
    p.resize(M+1, N+1);
    for (unsigned int m=0; m<=M; m++)
    {
        for (unsigned int n=0; n<=N; n++)
        {
            p[m][n] = p20[m][n];
        }
    }

    p00.clear();
    p10.clear();
    p15.clear();
    p20.clear();

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

void ConjugateCdIHyperbolicIBVP::implicit_calculate_D2V3(DoubleMatrix &p, double a, double alpha) const
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

    DoubleMatrix p00(M+1, N+1);
    DoubleMatrix p10(M+1, N+1);
    DoubleMatrix p15(M+1, N+1);
    DoubleMatrix p20(M+1, N+1);

    /***********************************************************************************************/
    initial_calculate(p00, p10, N, hx, M, hy, ht, a, alpha);
    /***********************************************************************************************/

    SpaceNodePDE sn;
    const unsigned int size_ln = static_cast<unsigned int>(0)-1;
    for (unsigned int ln=L-1; ln!=size_ln; ln--)
    {
        TimeNodePDE tn10; tn10.i = ln;   tn10.t = tn10.i*ht;
        TimeNodePDE tn15; tn15.i = ln;   tn15.t = tn15.i*ht-0.5*ht;
        TimeNodePDE tn20; tn20.i = ln-1; tn20.t = tn20.i*ht;
        /**************************************************** border conditions ***************************************************/
        border1_calculate(p15, p20, N, hx, M, hy, tn15, tn20);
        /**************************************************** border conditions ***************************************************/
        /**************************************************** x direction apprx ***************************************************/
        for (unsigned int m=1; m<=M-1; m++)
        {
            sn.j = static_cast<int>(m); sn.y = m*hy;
            for (unsigned int n=1; n<=N-1; n++)
            {
                sn.i = static_cast<int>(n); sn.x = n*hx;
                dx[n-1] = 0.0;
                dx[n-1] += p_aa_htht__hyhy*(p10[m-1][n] - 2.0*p10[m][n] + p10[m+1][n]);
                dx[n-1] += 2.0*p10[m][n] + (p10[m][n] - p00[m][n]);
                dx[n-1] += alpha_ht*(p10[m][n] - 0.5*(p10[m][n]-p00[m][n]));
                dx[n-1] += ht_ht*f(sn, tn10);
            }
            dx[0]   -= m_aa_htht__hxhx * p15[m][0];
            dx[N-2] -= m_aa_htht__hxhx * p15[m][N];
            tomasAlgorithm(ax, bx, cx, dx, rx, N-1);
            for (unsigned int n=1; n<=N-1; n++) p15[m][n] = rx[n-1];
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
                dy[m-1] += p_aa_htht__hxhx*(p15[m][n-1] - 2.0*p15[m][n] + p15[m][n+1]);
                dy[m-1] += 2.0*p15[m][n] + (p10[m][n]-p00[m][n]);
                dy[m-1] += alpha_ht*(p15[m][n] - 0.5*(p10[m][n]-p00[m][n]));
                dy[m-1] += ht_ht*f(sn, tn10);
            }
            dy[0]   -= m_aa_htht__hyhy * p20[0][n];
            dy[M-2] -= m_aa_htht__hyhy * p20[M][n];
            tomasAlgorithm(ay, by, cy, dy, ry, M-1);
            for (unsigned int m=1; m<=M-1; m++) p20[m][n] = ry[m-1];
        }
        /**************************************************** y direction apprx ***************************************************/
        for (unsigned int m=0; m<=M; m++)
        {
            for (unsigned int n=0; n<=N; n++)
            {
                p00[m][n] = p10[m][n];
                p10[m][n] = p20[m][n];
            }
        }

        layerInfo(p20, ln-1);
    }

    p.clear();
    p.resize(M+1, N+1);
    for (unsigned int m=0; m<=M; m++)
    {
        for (unsigned int n=0; n<=N; n++)
        {
            p[m][n] = p20[m][n];
        }
    }

    p00.clear();
    p10.clear();
    p15.clear();
    p20.clear();

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

void ConjugateCdIHyperbolicIBVP::initial_calculate(DoubleMatrix &p00, DoubleMatrix &p10, unsigned int N, double hx, unsigned int M, double hy, double ht, double a, double alpha) const
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
            p00[m][n] = initial1(sn);
        }
    }
    layerInfo(p00, L);

    /***********************************************************************************************/

    TimeNodePDE tn00; tn00.i = L-0; tn00.t = tn00.i*ht;
    TimeNodePDE tn10; tn10.i = L-1; tn10.t = tn10.i*ht;

    border2_calculate(p10, N, hx, M, hy, tn10);

    for (unsigned int m=1; m<=M-1; m++)
    {
        sn.j = static_cast<int>(m); sn.y = m*hy;
        for (unsigned int n=1; n<=N-1; n++)
        {
            sn.i = static_cast<int>(n); sn.x = n*hx;
            p10[m][n] = p00[m][n] - ht*initial2(sn);
            p10[m][n] += htht_05*(aa__hxhx*(p00[m][n-1]-2.0*p00[m][n]+p00[m][n+1])+aa__hyhy*(p00[m-1][n]-2.0*p00[m][n]+p00[m+1][n])+alpha*initial2(sn)+f(sn,tn00));
        }
    }
    layerInfo(p10, L-1);

    /***********************************************************************************************/
}

void ConjugateCdIHyperbolicIBVP::border2_calculate(DoubleMatrix &p, unsigned int N, double hx, unsigned int M, double hy, const TimeNodePDE &tn) const
{
    SpaceNodePDE sn0, sn1;

    sn0.i = static_cast<int>(0); sn0.x = 0*hx;
    sn1.i = static_cast<int>(N); sn1.x = N*hx;
    for (unsigned int m=0; m<=M; m++)
    {
        sn0.j = static_cast<int>(m); sn0.y = m*hy; p[m][0] = boundary(sn0, tn);
        sn1.j = static_cast<int>(m); sn1.y = m*hy; p[m][N] = boundary(sn1, tn);
    }

    sn0.j = static_cast<int>(0); sn0.y = 0*hy;
    sn1.j = static_cast<int>(M); sn1.y = M*hy;
    for (unsigned int n=0; n<=N; n++)
    {
        sn0.i = static_cast<int>(n); sn0.x = n*hx; p[0][n] = boundary(sn0, tn);
        sn1.i = static_cast<int>(n); sn1.x = n*hx; p[M][n] = boundary(sn1, tn);
    }
}

void ConjugateCdIHyperbolicIBVP::border1_calculate(DoubleMatrix &p15, DoubleMatrix &p20, unsigned int N, double hx, unsigned int M, double hy, const TimeNodePDE &tn15, const TimeNodePDE &tn20) const
{
    SpaceNodePDE sn0, sn1;

    sn0.i = static_cast<int>(0); sn0.x = 0*hx;
    sn1.i = static_cast<int>(N); sn1.x = N*hx;
    for (unsigned int m=0; m<=M; m++)
    {
        sn0.j = static_cast<int>(m); sn0.y = m*hy; p15[m][0] = boundary(sn0, tn15); p20[m][0] = boundary(sn0, tn20);
        sn1.j = static_cast<int>(m); sn1.y = m*hy; p15[m][N] = boundary(sn1, tn15); p20[m][N] = boundary(sn1, tn20);
    }

    sn0.j = static_cast<int>(0); sn0.y = 0*hy;
    sn1.j = static_cast<int>(M); sn1.y = M*hy;
    for (unsigned int n=0; n<=N; n++)
    {
        sn0.i = static_cast<int>(n); sn0.x = n*hx; p15[0][n] = boundary(sn0, tn15); p20[0][n] = boundary(sn0, tn20);
        sn1.i = static_cast<int>(n); sn1.x = n*hx; p15[M][n] = boundary(sn1, tn15); p20[M][n] = boundary(sn1, tn20);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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


