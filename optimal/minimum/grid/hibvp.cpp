#include "hibvp.h"

//--------------------------------------------------------------------------------------------------------------//

IWaveEquationIBVP::IWaveEquationIBVP(double waveSpeed, double waveDissipation) :
    _waveSpeed(waveSpeed), _waveDissipation(waveDissipation) {}

IWaveEquationIBVP::~IWaveEquationIBVP() {}

auto IWaveEquationIBVP::waveSpeed() const -> double { return _waveSpeed; }

auto IWaveEquationIBVP::setWaveSpeed(double waveSpeed) -> void { _waveSpeed = waveSpeed; }

auto IWaveEquationIBVP::waveDissipation() const -> double { return _waveDissipation; }

auto IWaveEquationIBVP::setWaveDissipation(double waveDissipation) -> void { _waveDissipation = waveDissipation; }

auto IWaveEquationIBVP::lambda() const -> double { return 0.25; }

/********************/

void IWaveEquationIBVP::explicit_calculate_D1V1(DoubleVector &, double, double) const {}

void IWaveEquationIBVP::explicit_calculate_D1V1_initial(DoubleVector &u00, DoubleVector &u10, unsigned int N, double hx, double ht, double a, double alpha) const
{
    const double aa__hxhx = (a*a)/(hx*hx);
    const double htht_05 = 0.5*ht*ht;
    TimeNodePDE tn00; tn00.i = 0; tn00.t = 0.0;
    TimeNodePDE tn10; tn10.i = 1; tn10.t = ht;

    /***********************************************************************************************/

    SpaceNodePDE sn; sn.j = 0; sn.y = 0.0;
    for (unsigned int n=0; n<=N; n++)
    {
        sn.i = static_cast<int>(n); sn.x = sn.i*hx;
        u00[n] = initial(sn, InitialCondition::InitialValue);
    }
    layerInfo(u00, tn00);

    /***********************************************************************************************/

    sn.i = static_cast<int>(0); sn.x = 0.0;   u10[0] = boundary(sn, tn10);
    sn.i = static_cast<int>(N); sn.x = N*hx;  u10[N] = boundary(sn, tn10);
    for (unsigned int n=1; n<=N-1; n++)
    {
        sn.i = static_cast<int>(n); sn.x = n*hx;
        u10[n] = u00[n] + ht*initial(sn, InitialCondition::FirstDerivative);
        u10[n] += htht_05*(aa__hxhx*(u00[n-1]-2.0*u00[n]+u00[n+1])+f(sn,tn00));
        u10[n] += htht_05*(-alpha*initial(sn,InitialCondition::FirstDerivative));
    }
    layerInfo(u10, tn10);

    /***********************************************************************************************/
}

void IWaveEquationIBVP::explicit_calculate_D1V1_border(DoubleVector &u, unsigned int N, double hx, double ht, const TimeNodePDE &tn) const
{}

void IWaveEquationIBVP::implicit_calculate_D1V1_initial(DoubleVector &u00, DoubleVector &u10, unsigned int N, double hx, double ht, double a, double alpha) const
{}

void IWaveEquationIBVP::implicit_calculate_D1V1_border(DoubleVector &u20, unsigned int N, double hx, const TimeNodePDE &tn20) const
{
    SpaceNodePDE sn0, sn1;
    sn0.i = static_cast<int>(0); sn0.x = 0*hx;
    sn1.i = static_cast<int>(N); sn1.x = N*hx;
    u20[0] = boundary(sn0, tn20); u20[N] = boundary(sn1, tn20);
}

/********************/

void IWaveEquationIBVP::implicit_calculate_D1V1() const
{
    const Dimension &dimX = spaceDimension(Dimension::DimensionX);
    const Dimension &time = timeDimension();

    const unsigned int N = static_cast<unsigned int>( dimX.size() );
    const unsigned int L = static_cast<unsigned int>( time.size() );

    const double hx = dimX.step();
    const double ht = time.step();

    const double alpha       = waveDissipation();
    const double a           = waveSpeed();
    const double lmbd        = lambda();
    const double m1_2lambda  = 1.0-2.0*lmbd;
    const double ht_ht       = ht*ht;
    const double alpha_ht_05 = 0.5*alpha*ht;

    const double m_aa_htht__hxhx_lambda = -(a*a)*((ht*ht)/(hx*hx))*lmbd;
    const double b_aa_htht__hxhx = +(1.0 + 2.0*(a*a)*((ht*ht)/(hx*hx))*lmbd + alpha_ht_05);
    const double p_aa_htht__hxhx_1m2lambda = +(a*a)*((ht*ht)/(hx*hx))*(1.0-2.0*lmbd);
    const double p_aa_htht__hxhx_lambda = +(a*a)*((ht*ht)/(hx*hx))*lmbd;

    double *ax = static_cast<double*>(malloc(sizeof(double)*(N+1)));
    double *bx = static_cast<double*>(malloc(sizeof(double)*(N+1)));
    double *cx = static_cast<double*>(malloc(sizeof(double)*(N+1)));
    double *dx = static_cast<double*>(malloc(sizeof(double)*(N+1)));
    double *rx = static_cast<double*>(malloc(sizeof(double)*(N+1)));

    for (unsigned int n=0; n<=N; n++)
    {
        ax[n] = m_aa_htht__hxhx_lambda;
        bx[n] = b_aa_htht__hxhx;
        cx[n] = m_aa_htht__hxhx_lambda;
    }
    ax[0] = 0.0; cx[N] = 0.0;

    const double aa__hxhx = (a*a)/(hx*hx);
    const double htht_05 = 0.5*ht*ht;

    DoubleVector u00(N+1);
    DoubleVector u10(N+1);
    DoubleVector u20(N+1);

    /***********************************************************************************************/

    TimeNodePDE tn00; tn00.i = 0; tn00.t = 0.5*tn00.i*ht;
    TimeNodePDE tn05; tn05.i = 1; tn05.t = 0.5*tn05.i*ht;
    TimeNodePDE tn10; tn10.i = 2; tn10.t = 0.5*tn10.i*ht;

    SpaceNodePDE sn;
    for (unsigned int n=0; n<=N; n++)
    {
        sn.i = static_cast<int>(n); sn.x = sn.i*hx;
        u00[n] = initial(sn, InitialCondition::InitialValue);
    }
    layerInfo(u00, tn00);

    /***********************************************************************************************/

    implicit_calculate_D1V1_border(u10, N, hx, tn10);

    for (unsigned int n=1; n<=N-1; n++)
    {
        sn.i = static_cast<int>(n); sn.x = sn.i*hx;
        double firstDerivative = initial(sn, InitialCondition::FirstDerivative);
        double u0tt = aa__hxhx*(u00[n-1]-2.0*u00[n]+u00[n+1]) + f(sn,tn00) - alpha*firstDerivative;
        u10[n] = u00[n] + firstDerivative*ht;
        u10[n] += u0tt * htht_05;
    }
    std::cout << "OK" << std::endl;
    layerInfo(u10, tn05);

    /***********************************************************************************************/

    //SpaceNodePDE sn;
    for (unsigned int ln=2; ln<=L; ln++)
    {
        TimeNodePDE tn10; tn10.i = ln-1; tn10.t = tn10.i*ht;
        TimeNodePDE tn20; tn20.i = ln;   tn20.t = tn20.i*ht;
        /**************************************************** border conditions ***************************************************/
        //implicit_calculate_D1V1_border(u20, N, hx, tn20);

        SpaceNodePDE sn0, sn1;
        sn0.i = static_cast<int>(0); sn0.x = 0*hx;
        sn1.i = static_cast<int>(N); sn1.x = N*hx;

        BoundaryConditionPDE cn0, cn1;
        double b0, bN;
        //b0 = boundary(sn0, tn20, cn0);
        //bN = boundary(sn1, tn20, cn1);
        unsigned int s,e;

        if (cn0.type == BoundaryCondition::Dirichlet)
        {
            s = 1;
            u20[0] = b0;
            ax[s] = 0.0;
            bx[s] = b_aa_htht__hxhx;
            cx[s] = m_aa_htht__hxhx_lambda;
        }
        else
        {
            s = 0;
            ax[s] = 0.0;
            bx[s] = 1.0 + (2.0*lmbd*(a*a)*(ht*ht))/(hx*hx) + alpha_ht_05
                    - 2.0*a*a*ht*ht*(lmbd/hx)*(cn0.coefficientValue/cn0.coefficientDerivative);
            cx[s] = -(2.0*lmbd*(a*a)*(ht*ht))/(hx*hx);
        }

        if (cn1.type == BoundaryCondition::Dirichlet)
        {
            e = N-1;
            u20[N] = bN;
            ax[e] = m_aa_htht__hxhx_lambda;
            bx[e] = b_aa_htht__hxhx;
            cx[e] = 0.0;
        }
        else
        {
            e = N;
            ax[e] = -(2.0*lmbd*(a*a)*(ht*ht))/(hx*hx);
            bx[e] = 1.0 + (2.0*lmbd*(a*a)*(ht*ht))/(hx*hx) + alpha_ht_05
                    - 2.0*a*a*ht*ht*(lmbd/hx)*(cn1.coefficientValue/cn1.coefficientDerivative);
            cx[e] = 0.0;
        }

        /**************************************************** border conditions ***************************************************/
        /**************************************************** x direction apprx ***************************************************/
        for (unsigned int n=1; n<=N-1; n++)
        {
            sn.i = static_cast<int>(n); sn.x = n*hx;
            dx[n] = 0.0;
            dx[n] += (u10[n-1] - 2.0*u10[n] + u10[n+1])*p_aa_htht__hxhx_1m2lambda;
            dx[n] += (u00[n-1] - 2.0*u00[n] + u00[n+1])*p_aa_htht__hxhx_lambda;
            dx[n] += 2.0*u10[n] - u00[n] + alpha_ht_05*u00[n];
            dx[n] += ht_ht*f(sn, tn10);
            //dx[n-1] += ht_ht*(lmbd*f(sn, tn05) + m1_2lambda*f(sn, tn10) + lmbd*f(sn, tn15));
        }

        if (cn0.type == BoundaryCondition::Dirichlet)
        {
            dx[1]   -= u20[0]*m_aa_htht__hxhx_lambda;
        }
        else
        {
            dx[0] = 0.0;
            sn.i = static_cast<int>(0); sn.x = 0*hx;
            dx[0] += ((u10[1] - u10[0])/hx+(cn0.coefficientValue/cn0.coefficientDerivative)*u10[0])*((2.0*a*a*ht*ht*(1.0-2.0*lmbd))/(hx));
            dx[0] += ((u00[1] - u00[0])/hx+(cn0.coefficientValue/cn0.coefficientDerivative)*u00[0])*((2.0*a*a*ht*ht*lmbd)/(hx));
            dx[0] += 2.0*u10[0] - u00[0] + alpha_ht_05*u00[0];
            dx[0] += -(2.0*a*a*ht*ht*b0)/(cn0.coefficientDerivative*hx);
            dx[0] += ht_ht*f(sn, tn10);
        }

        if (cn1.type == BoundaryCondition::Dirichlet)
        {
            dx[N-1] -= u20[N]*m_aa_htht__hxhx_lambda;
        }
        else
        {
            dx[N] = 0.0;
            sn.i = static_cast<int>(N); sn.x = N*hx;
            dx[N] += ((u10[N-1] - u10[N])/hx+(cn1.coefficientValue/cn1.coefficientDerivative)*u10[N])*((2.0*a*a*ht*ht*(1.0-2.0*lmbd))/(hx));
            dx[N] += ((u00[N-1] - u00[N])/hx+(cn1.coefficientValue/cn1.coefficientDerivative)*u00[N])*((2.0*a*a*ht*ht*lmbd)/(hx));
            dx[N] += 2.0*u10[N] - u00[N] + alpha_ht_05*u00[N];
            dx[N] += -(2.0*a*a*ht*ht*bN)/(cn1.coefficientDerivative*hx);
            dx[N] += ht_ht*f(sn, tn10);
        }

        //std::cout << s << " " << e << std::endl;
        tomasAlgorithm(ax+s, bx+s, cx+s, dx+s, rx+s, e-s+1);
        for (unsigned int n=1; n<=N-1; n++) u20[n] = rx[n];
        layerInfo(u20, tn20);
        /**************************************************** x direction apprx ***************************************************/
        u00 = u10;
        u10 = u20;
    }
}

/********************/

void IWaveEquationIBVP::explicit_calculate_D2V1() const
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

    const double a = waveSpeed();
    const double alpha = waveDissipation();

    const double alpha_ht_05 = alpha*ht*0.5;
    const double inv__1malpha_ht05 = 1.0/(1.0 + alpha_ht_05);
    const double aa_htht__hxhx = ((a*a*ht*ht)/(hx*hx));
    const double aa_htht__hyhy = ((a*a*ht*ht)/(hy*hy));

    DoubleMatrix u00(M+1, N+1);
    DoubleMatrix u10(M+1, N+1);
    DoubleMatrix u20(M+1, N+1);

    /**********************************************************************/
    explicit_calculate_D2V1_initial(u00, u10, N, hx, M, hy, ht, a, alpha);
    /**********************************************************************/

    for (unsigned int ln=2; ln<=L; ln++)
    {
        TimeNodePDE tn20; tn20.i = ln;   tn20.t = tn20.i*ht;
        TimeNodePDE tn10; tn10.i = ln-1; tn10.t = tn10.i*ht;

        /**************************** border conditions ****************************/
        explicit_calculate_D2V1_border(u20, N, hx, M, hy, tn20);
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
                        + 2.0*u10[m][n] - u00[m][n]+ alpha_ht_05*u00[m][n]
                        + ht*ht*f(sn,tn10))*inv__1malpha_ht05;
            }
        }
        layerInfo(u20, tn20);

        for (unsigned int m=0; m<=M; m++)
        {
            for (unsigned int n=0; n<=N; n++)
            {
                u00[m][n] = u10[m][n];
                u10[m][n] = u20[m][n];
            }
        }
    }
}

void IWaveEquationIBVP::explicit_calculate_D2V1_initial(DoubleMatrix &u00, DoubleMatrix &u10, unsigned int N, double hx, unsigned int M, double hy, double ht, double a, double alpha) const
{
    const double aa__hxhx = (a*a)/(hx*hx);
    const double aa__hyhy = (a*a)/(hy*hy);
    const double htht_05 = 0.5*ht*ht;

    TimeNodePDE tn00; tn00.i = 0; tn00.t = 0.0;
    TimeNodePDE tn10; tn10.i = 1; tn10.t = ht;

    /***********************************************************************************************/

    SpaceNodePDE sn;
    for (unsigned int m=0; m<=M; m++)
    {
        sn.j = static_cast<int>(m); sn.y = sn.j*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = static_cast<int>(n); sn.x = sn.i*hx;
            u00[m][n] = initial(sn, InitialCondition::InitialValue);
        }
    }
    layerInfo(u00, tn00);

    /***********************************************************************************************/

    explicit_calculate_D2V1_border(u10, N, hx, M, hy, tn10);

    for (unsigned int m=1; m<=M-1; m++)
    {
        sn.j = static_cast<int>(m); sn.y = m*hy;
        for (unsigned int n=1; n<=N-1; n++)
        {
            sn.i = static_cast<int>(n); sn.x = n*hx;
            u10[m][n] = u00[m][n] + ht*initial(sn, InitialCondition::FirstDerivative);
            u10[m][n] += htht_05*(aa__hxhx*(u00[m][n-1]-2.0*u00[m][n]+u00[m][n+1])+aa__hyhy*(u00[m-1][n]-2.0*u00[m][n]+u00[m+1][n])+f(sn,tn00)-alpha*initial(sn, InitialCondition::FirstDerivative));
        }
    }
    layerInfo(u10, tn10);

    /***********************************************************************************************/
}

void IWaveEquationIBVP::explicit_calculate_D2V1_border(DoubleMatrix &u, unsigned int N, double hx, unsigned int M, double hy, const TimeNodePDE &tn) const
{
    SpaceNodePDE sn0, sn1;

    sn0.i = static_cast<int>(0); sn0.x = 0*hx;
    sn1.i = static_cast<int>(N); sn1.x = N*hx;
    for (unsigned int m=0; m<=M; m++)
    {
        sn0.j = static_cast<int>(m); sn0.y = m*hy; u[m][0] = boundary(sn0, tn);
        sn1.j = static_cast<int>(m); sn1.y = m*hy; u[m][N] = boundary(sn1, tn);
    }

    sn0.j = static_cast<int>(0); sn0.y = 0*hy;
    sn1.j = static_cast<int>(M); sn1.y = M*hy;
    for (unsigned int n=0; n<=N; n++)
    {
        sn0.i = static_cast<int>(n); sn0.x = n*hx; u[0][n] = boundary(sn0, tn);
        sn1.i = static_cast<int>(n); sn1.x = n*hx; u[M][n] = boundary(sn1, tn);
    }
}

/********************/

void IWaveEquationIBVP::implicit_calculate_D2V1() const
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

    const double _lambda = lambda();
    const double a = waveSpeed();
    const double alpha = waveDissipation();
    const double ht_ht_025 = ht*ht*0.25;
    const double alpha_ht_025 = alpha*ht*0.25;

    const double m_aa_htht__hxhx_025_lambda = -(0.25*a*a)*((ht*ht)/(hx*hx))*_lambda;
    const double b_aa_htht__hxhx = +(1.0 + 0.5*(a*a)*((ht*ht)/(hx*hx))*_lambda + alpha_ht_025);
    const double p_aa_htht__hyhy_025 = +(0.25*a*a)*((ht*ht)/(hy*hy));
    const double p_aa_htht__hxhx_025_1m2lambda = +(0.25*a*a)*((ht*ht)/(hx*hx))*(1.0-2.0*_lambda);
    const double p_aa_htht__hxhx_025_lambda = +(0.25*a*a)*((ht*ht)/(hx*hx))*_lambda;

    const double m_aa_htht__hyhy_025_lambda = -(0.25*a*a)*((ht*ht)/(hy*hy))*_lambda;
    const double b_aa_htht__hyhy = +(1.0 + 0.5*(a*a)*((ht*ht)/(hy*hy))*_lambda + alpha_ht_025);
    const double p_aa_htht__hxhx_025 = +(0.25*a*a)*((ht*ht)/(hx*hx));
    const double p_aa_htht__hyhy_025_1m2lambda = +(0.25*a*a)*((ht*ht)/(hy*hy))*(1.0-2.0*_lambda);
    const double p_aa_htht__hyhy_025_lambda = +(0.25*a*a)*((ht*ht)/(hy*hy))*_lambda;

    double *ax = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *bx = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *cx = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *dx = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *rx = static_cast<double*>(malloc(sizeof(double)*(N-1)));

    for (unsigned int n=0; n<=N-2; n++)
    {
        ax[n] = m_aa_htht__hxhx_025_lambda;
        bx[n] = b_aa_htht__hxhx;
        cx[n] = m_aa_htht__hxhx_025_lambda;
    }
    ax[0] = 0.0; cx[N-2] = 0.0;

    double *ay = static_cast<double*>(malloc(sizeof(double)*(M-1)));
    double *by = static_cast<double*>(malloc(sizeof(double)*(M-1)));
    double *cy = static_cast<double*>(malloc(sizeof(double)*(M-1)));
    double *dy = static_cast<double*>(malloc(sizeof(double)*(M-1)));
    double *ry = static_cast<double*>(malloc(sizeof(double)*(M-1)));

    for (unsigned int m=0; m<=M-2; m++)
    {
        ay[m] = m_aa_htht__hyhy_025_lambda;
        by[m] = b_aa_htht__hyhy;
        cy[m] = m_aa_htht__hyhy_025_lambda;
    }
    ay[0] = 0.0; cy[M-2] = 0.0;

    DoubleMatrix u00(M+1, N+1);
    DoubleMatrix u05(M+1, N+1);
    DoubleMatrix u10(M+1, N+1);
    DoubleMatrix u15(M+1, N+1);
    DoubleMatrix u20(M+1, N+1);

    /***********************************************************************************************/
    implicit_calculate_D2V1_initial(u00, u05, u10, N, hx, M, hy, ht, a, alpha);
    /***********************************************************************************************/

    SpaceNodePDE sn;
    for (unsigned int ln=2; ln<=L; ln++)
    {
        TimeNodePDE tn05; tn05.i = 2*ln-3; tn05.t = 0.5*tn05.i*ht;
        TimeNodePDE tn10; tn10.i = 2*ln-2; tn10.t = 0.5*tn10.i*ht;
        TimeNodePDE tn15; tn15.i = 2*ln-1; tn15.t = 0.5*tn15.i*ht;
        TimeNodePDE tn20; tn20.i = 2*ln-0; tn20.t = 0.5*tn20.i*ht;
        /**************************************************** border conditions ***************************************************/
        implicit_calculate_D2V1_border(u15, u20, N, hx, M, hy, tn15, tn20);
        /**************************************************** border conditions ***************************************************/
        /**************************************************** x direction apprx ***************************************************/
        for (unsigned int m=1; m<=M-1; m++)
        {
            sn.j = static_cast<int>(m); sn.y = m*hy;
            for (unsigned int n=1; n<=N-1; n++)
            {
                sn.i = static_cast<int>(n); sn.x = n*hx;
                dx[n-1] = 0.0;
                dx[n-1] += (u10[m][n-1] - 2.0*u10[m][n] + u10[m][n+1])*p_aa_htht__hxhx_025_1m2lambda;
                dx[n-1] += (u05[m][n-1] - 2.0*u05[m][n] + u05[m][n+1])*p_aa_htht__hxhx_025_lambda;
                dx[n-1] += (u10[m-1][n] - 2.0*u10[m][n] + u10[m+1][n])*p_aa_htht__hyhy_025;
                dx[n-1] += 2.0*u10[m][n] - u05[m][n] + alpha_ht_025*u05[m][n];
                dx[n-1] += ht_ht_025*f(sn, tn10);
                //dx[n-1] += ht_ht_025*(lambda*f(sn, tn05)+(1.0-2.0*lambda)*f(sn, tn10)+lambda*f(sn, tn15));
            }
            dx[0]   -= u15[m][0]*m_aa_htht__hxhx_025_lambda;
            dx[N-2] -= u15[m][N]*m_aa_htht__hxhx_025_lambda;
            tomasAlgorithm(ax, bx, cx, dx, rx, N-1);
            for (unsigned int n=1; n<=N-1; n++) u15[m][n] = rx[n-1];
        }
        layerInfo(u15, tn15);
        /**************************************************** x direction apprx ***************************************************/
        /**************************************************** y direction apprx ***************************************************/
        for (unsigned int n=1; n<=N-1; n++)
        {
            sn.i = static_cast<int>(n); sn.x = n*hx;
            for (unsigned int m=1; m<=M-1; m++)
            {
                sn.j = static_cast<int>(m); sn.y = m*hy;
                dy[m-1] = 0.0;
                dy[m-1] += (u15[m-1][n] - 2.0*u15[m][n] + u15[m+1][n])*p_aa_htht__hyhy_025_1m2lambda;
                dy[m-1] += (u10[m-1][n] - 2.0*u10[m][n] + u10[m+1][n])*p_aa_htht__hyhy_025_lambda;
                dy[m-1] += (u15[m][n-1] - 2.0*u15[m][n] + u15[m][n+1])*p_aa_htht__hxhx_025;
                dy[m-1] += 2.0*u15[m][n] - u10[m][n] + alpha_ht_025*u10[m][n];
                dy[m-1] += ht_ht_025*f(sn, tn15);
                //dy[m-1] += ht_ht_025*(lambda*f(sn, tn10)+(1.0-2.0*lambda)*f(sn, tn15)+lambda*f(sn, tn20));
            }
            dy[0]   -= u20[0][n]*m_aa_htht__hyhy_025_lambda;
            dy[M-2] -= u20[M][n]*m_aa_htht__hyhy_025_lambda;
            tomasAlgorithm(ay, by, cy, dy, ry, M-1);
            for (unsigned int m=1; m<=M-1; m++) u20[m][n] = ry[m-1];
        }
        layerInfo(u20, tn20);
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

    u00.clear();
    u05.clear();
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

void IWaveEquationIBVP::implicit_calculate_D2V1_initial(DoubleMatrix &u00, DoubleMatrix &u05, DoubleMatrix &u10, unsigned int N, double hx, unsigned int M, double hy, double ht, double a, double alpha) const
{
    const double aa__hxhx = (a*a)/(hx*hx);
    const double aa__hyhy = (a*a)/(hy*hy);
    const double ht_05 = ht*0.5;
    const double htht_05 = 0.5*ht*ht;
    const double htht_05_025 = htht_05*0.25;

    TimeNodePDE tn00; tn00.i = 0; tn00.t = 0.5*tn00.i*ht;
    TimeNodePDE tn05; tn05.i = 1; tn05.t = 0.5*tn05.i*ht;
    TimeNodePDE tn10; tn10.i = 2; tn10.t = 0.5*tn10.i*ht;

    /***********************************************************************************************/

    SpaceNodePDE sn;
    for (unsigned int m=0; m<=M; m++)
    {
        sn.j = static_cast<int>(m); sn.y = sn.j*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = static_cast<int>(n); sn.x = sn.i*hx;
            u00[m][n] = initial(sn, InitialCondition::InitialValue);
        }
    }
    layerInfo(u00, tn00);

    /***********************************************************************************************/

    implicit_calculate_D2V1_border(u05, u10, N, hx, M, hy, tn05, tn10);

    for (unsigned int m=1; m<=M-1; m++)
    {
        sn.j = static_cast<int>(m); sn.y = sn.j*hy;
        for (unsigned int n=1; n<=N-1; n++)
        {
            sn.i = static_cast<int>(n); sn.x = sn.i*hx;
            double firstDerivative = initial(sn, InitialCondition::FirstDerivative);
            double secndDerivative = 0.0;
            secndDerivative += aa__hxhx*(u00[m][n-1]-2.0*u00[m][n]+u00[m][n+1]);
            secndDerivative += aa__hyhy*(u00[m-1][n]-2.0*u00[m][n]+u00[m+1][n]);
            secndDerivative += f(sn,tn00);
            secndDerivative -=  alpha*firstDerivative;

            u05[m][n] = u00[m][n] + firstDerivative*ht_05 + secndDerivative*htht_05_025;
            //u10[m][n] = u00[m][n] + firstDerivative*ht    + secndDerivative*htht_05;
        }
    }
    layerInfo(u05, tn05);

    for (unsigned int m=1; m<=M-1; m++)
    {
        sn.j = static_cast<int>(m); sn.y = sn.j*hy;
        for (unsigned int n=1; n<=N-1; n++)
        {
            sn.i = static_cast<int>(n); sn.x = sn.i*hx;
            u10[m][n]  = 2.0*u05[m][n]-u00[m][n];
            u10[m][n] += ((a*a*ht*ht*0.25)/(hx*hx))*(u05[m][n-1]-2.0*u05[m][n]+u05[m][n+1]);
            u10[m][n] += ((a*a*ht*ht*0.25)/(hy*hy))*(u05[m-1][n]-2.0*u05[m][n]+u05[m+1][n]);
            u10[m][n] += (ht*ht*0.25)*f(sn,tn05);
            u10[m][n] += (ht*0.25)*alpha*u00[m][n];
            u10[m][n] /= (1.0+0.25*ht*alpha);
        }
    }

    layerInfo(u10, tn10);
}

void IWaveEquationIBVP::implicit_calculate_D2V1_border(DoubleMatrix &u05, DoubleMatrix &u10, unsigned int N, double hx, unsigned int M, double hy, const TimeNodePDE &tn05, const TimeNodePDE &tn10) const
{
    SpaceNodePDE sn0, sn1;

    sn0.i = static_cast<int>(0); sn0.x = 0*hx;
    sn1.i = static_cast<int>(N); sn1.x = N*hx;
    for (unsigned int m=0; m<=M; m++)
    {
        sn0.j = static_cast<int>(m); sn0.y = m*hy; u05[m][0] = boundary(sn0, tn05); u10[m][0] = boundary(sn0, tn10);
        sn1.j = static_cast<int>(m); sn1.y = m*hy; u05[m][N] = boundary(sn1, tn05); u10[m][N] = boundary(sn1, tn10);
    }

    sn0.j = static_cast<int>(0); sn0.y = 0*hy;
    sn1.j = static_cast<int>(M); sn1.y = M*hy;
    for (unsigned int n=0; n<=N; n++)
    {
        sn0.i = static_cast<int>(n); sn0.x = n*hx; u05[0][n] = boundary(sn0, tn05); u10[0][n] = boundary(sn0, tn10);
        sn1.i = static_cast<int>(n); sn1.x = n*hx; u05[M][n] = boundary(sn1, tn05); u10[M][n] = boundary(sn1, tn10);
    }
}

/********************/

ConjugateCdIHyperbolicIBVP::~ConjugateCdIHyperbolicIBVP() {}

void ConjugateCdIHyperbolicIBVP::explicit_calculate_D2V1(DoubleMatrix &p, double a, double alpha) const
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

    const double alpha_ht_05 = alpha*ht*0.5;
    const double inv__1malpha_ht05 = 1.0/(1.0 - alpha_ht_05);
    const double aa_htht__hxhx = ((a*a*ht*ht)/(hx*hx));
    const double aa_htht__hyhy = ((a*a*ht*ht)/(hy*hy));

    DoubleMatrix p00(M+1, N+1);
    DoubleMatrix p10(M+1, N+1);
    DoubleMatrix p20(M+1, N+1);

    /**********************************************************************/
    explicit_calculate_D2V1_initial(p00, p10, N, hx, M, hy, ht, a, alpha, L);
    /**********************************************************************/

    const unsigned int size_ln = static_cast<unsigned int>(0)-1;
    for (unsigned int ln=L-2; ln != size_ln; ln--)
    {
        TimeNodePDE tn20; tn20.i = ln;   tn20.t = tn20.i*ht;
        TimeNodePDE tn10; tn10.i = ln+1; tn10.t = tn10.i*ht;

        /**************************** border conditions ****************************/
        explicit_calculate_D2V1_border(p20, N, hx, M, hy, tn20);
        /**************************** border conditions ****************************/

        SpaceNodePDE sn;
        for (unsigned int m=1; m<=M-1; m++)
        {
            sn.j = static_cast<int>(m); sn.y = m*hy;
            for (unsigned int n=1; n<=N-1; n++)
            {
                sn.i = static_cast<int>(n); sn.x = n*hx;
                p20[m][n] = (aa_htht__hxhx*(p10[m][n-1]-2.0*p10[m][n]+p10[m][n+1])
                        + aa_htht__hyhy*(p10[m-1][n]-2.0*p10[m][n]+p10[m+1][n])
                        + 2.0*p10[m][n] - p00[m][n] - alpha_ht_05*p00[m][n]
                        + ht*ht*f(sn,tn10)) * inv__1malpha_ht05;
            }
        }
        layerInfo(p20, tn20);

        for (unsigned int m=0; m<=M; m++)
        {
            for (unsigned int n=0; n<=N; n++)
            {
                p00[m][n] = p10[m][n];
                p10[m][n] = p20[m][n];
            }
        }
    }

    p.clear(); p.resize(M+1, N+1); for (unsigned int m=0; m<=M; m++) { for (unsigned int n=0; n<=N; n++) { p[m][n] = p20[m][n]; } }
}

void ConjugateCdIHyperbolicIBVP::explicit_calculate_D2V1_initial(DoubleMatrix &p00, DoubleMatrix &p10, unsigned int N, double hx, unsigned int M, double hy, double ht, double a, double alpha, unsigned int L) const
{
    const double aa__hxhx = (a*a)/(hx*hx);
    const double aa__hyhy = (a*a)/(hy*hy);
    const double htht_05 = 0.5*ht*ht;

    TimeNodePDE tn00; tn00.i = L;   tn00.t = tn00.i*ht;
    TimeNodePDE tn10; tn10.i = L-1; tn10.t = tn10.i*ht;

    /***********************************************************************************************/

    SpaceNodePDE sn;
    for (unsigned int m=0; m<=M; m++)
    {
        sn.j = static_cast<int>(m); sn.y = sn.j*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = static_cast<int>(n); sn.x = sn.i*hx;
            p00[m][n] = initial(sn, InitialCondition::InitialValue);
        }
    }
    layerInfo(p00, tn00);

    /***********************************************************************************************/

    explicit_calculate_D2V1_border(p10, N, hx, M, hy, tn10);

    for (unsigned int m=1; m<=M-1; m++)
    {
        sn.j = static_cast<int>(m); sn.y = m*hy;
        for (unsigned int n=1; n<=N-1; n++)
        {
            sn.i = static_cast<int>(n); sn.x = n*hx;
            p10[m][n] = p00[m][n] + ht*initial(sn, InitialCondition::FirstDerivative);
            p10[m][n] += htht_05*(aa__hxhx*(p00[m][n-1]-2.0*p00[m][n]+p00[m][n+1])+aa__hyhy*(p00[m-1][n]-2.0*p00[m][n]+p00[m+1][n])+f(sn,tn00)+alpha*initial(sn, InitialCondition::FirstDerivative));
        }
    }
    layerInfo(p10, tn10);

    /***********************************************************************************************/
}

void ConjugateCdIHyperbolicIBVP::explicit_calculate_D2V1_border(DoubleMatrix &p, unsigned int N, double hx, unsigned int M, double hy, const TimeNodePDE &tn) const
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

void ConjugateCdIHyperbolicIBVP::implicit_calculate_D2V1(DoubleMatrix &p, double a, double alpha, double lambda) const
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

    const double ht_ht_025 = ht*ht*0.25;
    const double alpha_ht_025 = alpha*ht*0.25;

    const double m_aa_htht__hxhx_025_lambda = -(0.25*a*a)*((ht*ht)/(hx*hx))*lambda;
    const double b_aa_htht__hxhx = +(1.0 + 0.5*(a*a)*((ht*ht)/(hx*hx))*lambda + alpha_ht_025);
    const double p_aa_htht__hyhy_025 = +(0.25*a*a)*((ht*ht)/(hy*hy));
    const double p_aa_htht__hxhx_025_1m2lambda = +(0.25*a*a)*((ht*ht)/(hx*hx))*(1.0-2.0*lambda);
    const double p_aa_htht__hxhx_025_lambda = +(0.25*a*a)*((ht*ht)/(hx*hx))*lambda;

    const double m_aa_htht__hyhy_025_lambda = -(0.25*a*a)*((ht*ht)/(hy*hy))*lambda;
    const double b_aa_htht__hyhy = +(1.0 + 0.5*(a*a)*((ht*ht)/(hy*hy))*lambda + alpha_ht_025);
    const double p_aa_htht__hxhx_025 = +(0.25*a*a)*((ht*ht)/(hx*hx));
    const double p_aa_htht__hyhy_025_1m2lambda = +(0.25*a*a)*((ht*ht)/(hy*hy))*(1.0-2.0*lambda);
    const double p_aa_htht__hyhy_025_lambda = +(0.25*a*a)*((ht*ht)/(hy*hy))*lambda;

    double *ax = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *bx = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *cx = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *dx = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *rx = static_cast<double*>(malloc(sizeof(double)*(N-1)));

    for (unsigned int n=0; n<=N-2; n++)
    {
        ax[n] = m_aa_htht__hxhx_025_lambda;
        bx[n] = b_aa_htht__hxhx;
        cx[n] = m_aa_htht__hxhx_025_lambda;
    }
    ax[0] = 0.0; cx[N-2] = 0.0;

    double *ay = static_cast<double*>(malloc(sizeof(double)*(M-1)));
    double *by = static_cast<double*>(malloc(sizeof(double)*(M-1)));
    double *cy = static_cast<double*>(malloc(sizeof(double)*(M-1)));
    double *dy = static_cast<double*>(malloc(sizeof(double)*(M-1)));
    double *ry = static_cast<double*>(malloc(sizeof(double)*(M-1)));

    for (unsigned int m=0; m<=M-2; m++)
    {
        ay[m] = m_aa_htht__hyhy_025_lambda;
        by[m] = b_aa_htht__hyhy;
        cy[m] = m_aa_htht__hyhy_025_lambda;
    }
    ay[0] = 0.0; cy[M-2] = 0.0;

    DoubleMatrix p00(M+1, N+1);
    DoubleMatrix p05(M+1, N+1);
    DoubleMatrix p10(M+1, N+1);
    DoubleMatrix p15(M+1, N+1);
    DoubleMatrix p20(M+1, N+1);

    /***********************************************************************************************/
    implicit_calculate_D2V1_initial(p00, p05, p10, N, hx, M, hy, ht, a, alpha, L);
    /***********************************************************************************************/

    SpaceNodePDE sn;
    const unsigned int size_ln = static_cast<unsigned int>(0)-1;
    for (unsigned int ln=L-2; ln !=size_ln; ln--)
    {
        TimeNodePDE tn05; tn05.i = 2*ln+3; tn05.t = 0.5*tn05.i*ht;//tn05.i*ht+0.5*ht;
        TimeNodePDE tn10; tn10.i = 2*ln+2; tn10.t = 0.5*tn10.i*ht;//tn10.i*ht;
        TimeNodePDE tn15; tn15.i = 2*ln+1; tn15.t = 0.5*tn15.i*ht;//tn15.i*ht-0.5*ht;
        TimeNodePDE tn20; tn20.i = 2*ln;   tn20.t = 0.5*tn20.i*ht;//tn20.i*ht;
        /**************************************************** border conditions ***************************************************/
        implicit_calculate_D2V1_border(p15, p20, N, hx, M, hy, tn15, tn20);
        /**************************************************** border conditions ***************************************************/
        /**************************************************** x direction apprx ***************************************************/
        for (unsigned int m=1; m<=M-1; m++)
        {
            sn.j = static_cast<int>(m); sn.y = m*hy;
            for (unsigned int n=1; n<=N-1; n++)
            {
                sn.i = static_cast<int>(n); sn.x = n*hx;
                dx[n-1] = 0.0;
                dx[n-1] += (p10[m][n-1] - 2.0*p10[m][n] + p10[m][n+1])*p_aa_htht__hxhx_025_1m2lambda;
                dx[n-1] += (p05[m][n-1] - 2.0*p05[m][n] + p05[m][n+1])*p_aa_htht__hxhx_025_lambda;
                dx[n-1] += (p10[m-1][n] - 2.0*p10[m][n] + p10[m+1][n])*p_aa_htht__hyhy_025;
                dx[n-1] += 2.0*p10[m][n] - p05[m][n] + alpha_ht_025*p05[m][n];
                //dx[n-1] += ht_ht_025*f(sn, tn10);
                dx[n-1] += ht_ht_025*(lambda*f(sn, tn05)+(1.0-2.0*lambda)*f(sn, tn10)+lambda*f(sn, tn15));
            }
            dx[0]   -= p15[m][0]*m_aa_htht__hxhx_025_lambda;
            dx[N-2] -= p15[m][N]*m_aa_htht__hxhx_025_lambda;
            tomasAlgorithm(ax, bx, cx, dx, rx, N-1);
            for (unsigned int n=1; n<=N-1; n++) p15[m][n] = rx[n-1];
        }
        layerInfo(p15, tn15);
        /**************************************************** x direction apprx ***************************************************/
        /**************************************************** y direction apprx ***************************************************/
        for (unsigned int n=1; n<=N-1; n++)
        {
            sn.i = static_cast<int>(n); sn.x = n*hx;
            for (unsigned int m=1; m<=M-1; m++)
            {
                sn.j = static_cast<int>(m); sn.y = m*hy;
                dy[m-1] = 0.0;
                dy[m-1] += (p15[m-1][n] - 2.0*p15[m][n] + p15[m+1][n])*p_aa_htht__hyhy_025_1m2lambda;
                dy[m-1] += (p10[m-1][n] - 2.0*p10[m][n] + p10[m+1][n])*p_aa_htht__hyhy_025_lambda;
                dy[m-1] += (p15[m][n-1] - 2.0*p15[m][n] + p15[m][n+1])*p_aa_htht__hxhx_025;
                dy[m-1] += 2.0*p15[m][n] - p10[m][n] + alpha_ht_025*p10[m][n];
                //dy[m-1] += ht_ht_025*f(sn, tn15);
                dy[m-1] += ht_ht_025*(lambda*f(sn, tn10)+(1.0-2.0*lambda)*f(sn, tn15)+lambda*f(sn, tn20));
            }
            dy[0]   -= p20[0][n]*m_aa_htht__hyhy_025_lambda;
            dy[M-2] -= p20[M][n]*m_aa_htht__hyhy_025_lambda;
            tomasAlgorithm(ay, by, cy, dy, ry, M-1);
            for (unsigned int m=1; m<=M-1; m++) p20[m][n] = ry[m-1];
        }
        layerInfo(p20, tn20);
        /**************************************************** y direction apprx ***************************************************/
        //p00 = p10;
        //p05 = p15;
        //p10 = p20;
        for (unsigned int m=0; m<=M; m++)
        {
            for (unsigned int n=0; n<=N; n++)
            {
                p00[m][n] = p10[m][n];
                p05[m][n] = p15[m][n];
                p10[m][n] = p20[m][n];
            }
        }
    }
    p.clear(); p.resize(M+1, N+1); for (unsigned int m=0; m<=M; m++) { for (unsigned int n=0; n<=N; n++) { p[m][n] = p20[m][n]; } }

    p00.clear();
    p05.clear();
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

void ConjugateCdIHyperbolicIBVP::implicit_calculate_D2V1_initial(DoubleMatrix &p00, DoubleMatrix &p05, DoubleMatrix &p10, unsigned int N, double hx, unsigned int M, double hy, double ht, double a, double alpha, unsigned int L) const
{
    const double aa__hxhx = (a*a)/(hx*hx);
    const double aa__hyhy = (a*a)/(hy*hy);
    const double ht_05 = 0.5*ht;
    const double htht_05 = 0.5*ht*ht;
    const double htht_05_025 = 0.5*ht*ht*0.25;

    TimeNodePDE tn00; tn00.i = 2*L;   tn00.t = 0.5*tn00.i*ht;
    TimeNodePDE tn05; tn05.i = 2*L-1; tn05.t = 0.5*tn05.i*ht;
    TimeNodePDE tn10; tn10.i = 2*L-2; tn10.t = 0.5*tn10.i*ht;

    /***********************************************************************************************/

    SpaceNodePDE sn;
    for (unsigned int m=0; m<=M; m++)
    {
        sn.j = static_cast<int>(m); sn.y = sn.j*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = static_cast<int>(n); sn.x = sn.i*hx;
            p00[m][n] = initial(sn, InitialCondition::InitialValue);
        }
    }
    layerInfo(p00, tn00);

    /***********************************************************************************************/

    implicit_calculate_D2V1_border(p05, p10, N, hx, M, hy, tn05, tn10);

    for (unsigned int m=1; m<=M-1; m++)
    {
        sn.j = static_cast<int>(m); sn.y = sn.j*hy;
        for (unsigned int n=1; n<=N-1; n++)
        {
            sn.i = static_cast<int>(n); sn.x = sn.i*hx;
            double firstDerivative = initial(sn, InitialCondition::FirstDerivative);
            double secndDericative = 0.0;
            secndDericative += aa__hxhx*(p00[m][n-1]-2.0*p00[m][n]+p00[m][n+1]);
            secndDericative += aa__hyhy*(p00[m-1][n]-2.0*p00[m][n]+p00[m+1][n]);
            secndDericative += f(sn,tn00);
            secndDericative += alpha*firstDerivative;

            p05[m][n] = p00[m][n] - firstDerivative*ht_05 + secndDericative*htht_05_025;
            //p10[m][n] = p00[m][n] - firstDerivative*ht    + secndDericative*htht_05;
        }
    }
    layerInfo(p05, tn05);

    for (unsigned int m=1; m<=M-1; m++)
    {
        sn.j = static_cast<int>(m); sn.y = sn.j*hy;
        for (unsigned int n=1; n<=N-1; n++)
        {
            sn.i = static_cast<int>(n); sn.x = sn.i*hx;
            p10[m][n]  = 2.0*p05[m][n]-p00[m][n];
            p10[m][n] += ((a*a*ht*ht*0.25)/(hx*hx))*(p05[m][n-1]-2.0*p05[m][n]+p05[m][n+1]);
            p10[m][n] += ((a*a*ht*ht*0.25)/(hy*hy))*(p05[m-1][n]-2.0*p05[m][n]+p05[m+1][n]);
            p10[m][n] += (ht*ht*0.25)*f(sn,tn05);
            p10[m][n] -= (ht*0.25)*alpha*p00[m][n];
            p10[m][n] /= (1.0-0.25*ht*alpha);
        }
    }

    layerInfo(p10, tn10);
}

void ConjugateCdIHyperbolicIBVP::implicit_calculate_D2V1_border(DoubleMatrix &p05, DoubleMatrix &p10, unsigned int N, double hx, unsigned int M, double hy, const TimeNodePDE &tn05, const TimeNodePDE &tn10) const
{
    SpaceNodePDE sn0, sn1;

    sn0.i = static_cast<int>(0); sn0.x = 0*hx;
    sn1.i = static_cast<int>(N); sn1.x = N*hx;
    for (unsigned int m=0; m<=M; m++)
    {
        sn0.j = static_cast<int>(m); sn0.y = m*hy; p05[m][0] = boundary(sn0, tn05); p10[m][0] = boundary(sn0, tn10);
        sn1.j = static_cast<int>(m); sn1.y = m*hy; p05[m][N] = boundary(sn1, tn05); p10[m][N] = boundary(sn1, tn10);
    }

    sn0.j = static_cast<int>(0); sn0.y = 0*hy;
    sn1.j = static_cast<int>(M); sn1.y = M*hy;
    for (unsigned int n=0; n<=N; n++)
    {
        sn0.i = static_cast<int>(n); sn0.x = n*hx; p05[0][n] = boundary(sn0, tn05); p10[0][n] = boundary(sn0, tn10);
        sn1.i = static_cast<int>(n); sn1.x = n*hx; p05[M][n] = boundary(sn1, tn05); p10[M][n] = boundary(sn1, tn10);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
