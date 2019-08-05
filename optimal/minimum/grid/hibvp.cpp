#include "hibvp.h"

//--------------------------------------------------------------------------------------------------------------//

IWaveEquationIBVP::IWaveEquationIBVP(double waveSpeed, double waveDissipation) :
    _waveSpeed(waveSpeed), _waveDissipation(waveDissipation) {}

IWaveEquationIBVP::IWaveEquationIBVP(const IWaveEquationIBVP& ibvp) :
    _waveSpeed(ibvp._waveSpeed), _waveDissipation(ibvp._waveDissipation)
{
    this->_timeDimension = ibvp._timeDimension;
    this->_spaceDimensionX = ibvp._spaceDimensionX;
    this->_spaceDimensionY = ibvp._spaceDimensionY;
    this->_spaceDimensionZ = ibvp._spaceDimensionZ;
}

IWaveEquationIBVP & IWaveEquationIBVP::operator =(const IWaveEquationIBVP &other)
{
    if (this == &other) { return *this; }

    this->_timeDimension = other._timeDimension;
    this->_spaceDimensionX = other._spaceDimensionX;
    this->_spaceDimensionY = other._spaceDimensionY;
    this->_spaceDimensionZ = other._spaceDimensionZ;
    this->_waveSpeed = other._waveSpeed;
    this->_waveDissipation = other._waveDissipation;
    return *this;
}

IWaveEquationIBVP::~IWaveEquationIBVP() {}

double IWaveEquationIBVP::waveSpeed() const { return _waveSpeed; }

double IWaveEquationIBVP::waveDissipation() const { return _waveDissipation; }

void IWaveEquationIBVP::setWaveSpeed(double waveSpeed) { _waveSpeed = waveSpeed; }

void IWaveEquationIBVP::setWaveDissipation(double waveDissipation) { _waveDissipation = waveDissipation; }

void IWaveEquationIBVP::explicit_calculate_D1V1() const
{
    const Dimension &dimX = spaceDimensionX();
    const Dimension &time = timeDimension();

    const unsigned int N = static_cast<unsigned int>( dimX.size() );
    const unsigned int L = static_cast<unsigned int>( time.size() );

    const double hx = dimX.step();
    const double ht = time.step();

    const double a               = waveSpeed();
    const double wd              = waveDissipation();
    const double ht_ht           = ht*ht;
    const double wvdsp_ht_05     = 0.5*wd*ht;
    const double p_aa_htht__hx   = ((a*a)*(ht*ht))/hx;
    const double p_aa_htht__hxhx = ((a*a)*(ht*ht))/(hx*hx);

    const double aa__hxhx = (a*a)/(hx*hx);
    const double htht_05 = 0.5*ht*ht;

    DoubleVector u00(N+1);
    DoubleVector u10(N+1);
    DoubleVector u20(N+1);

    /***********************************************************************************************/
    /***********************************************************************************************/

    TimeNodePDE tn00; tn00.i = 0; tn00.t = tn00.i*ht;
    TimeNodePDE tn10; tn10.i = 1; tn10.t = tn10.i*ht;

    SpaceNodePDE sn;
    for (unsigned int n=0; n<=N; n++)
    {
        sn.i = static_cast<int>(n); sn.x = sn.i*hx;
        u00[n] = initial(sn, InitialCondition::InitialValue);
    }
    layerInfo(u00, tn00);

    /***********************************************************************************************/

    for (unsigned int n=0; n<=N; n++)
    {
        sn.i = static_cast<int>(n); sn.x = sn.i*hx;
        double firstDerivative = initial(sn, InitialCondition::FirstDerivative);
        u10[n] = u00[n] + ht*firstDerivative;
        double uxx = 0.0;
        if (n==0)
        {
            uxx = aa__hxhx*(+2.0*u00[0]   -5.0*u00[1]   +4.0*u00[2]   -1.0*u00[3]);
        }
        else if (n==N)
        {
            uxx = aa__hxhx*(-1.0*u00[N-3] +4.0*u00[N-2] -5.0*u00[N-1] +2.0*u00[N]);
        }
        else
        {
            uxx = aa__hxhx*(u00[n-1]-2.0*u00[n]+u00[n+1]);
        }
        u10[n] += (uxx + f(sn,tn00) - wd*firstDerivative) * htht_05;
    }
    layerInfo(u10, tn10);

    /***********************************************************************************************/
    /***********************************************************************************************/

    for (unsigned int ln=2; ln<=L; ln++)
    {
        TimeNodePDE tn00; tn00.i = ln-2; tn00.t = tn00.i*ht;
        TimeNodePDE tn10; tn10.i = ln-1; tn10.t = tn10.i*ht;
        TimeNodePDE tn20; tn20.i = ln;   tn20.t = tn20.i*ht;

        BoundaryConditionPDE condition; double alpha, beta, gamma, value;

        for (unsigned int n=1; n<=N-1; n++)
        {
            sn.i = static_cast<int>(n); sn.x = n*hx;

            u20[n]  = 2.0*u10[n] - u00[n] + wvdsp_ht_05*u00[n] + ht_ht*f(sn, tn10);
            u20[n] += (u10[n-1] - 2.0*u10[n] + u10[n+1])*p_aa_htht__hxhx;
            u20[n] /= (1.0+wvdsp_ht_05);
        }

        sn.i = static_cast<int>(0); sn.x = 0*hx;
        value = boundary(sn, tn10, condition);
        alpha = condition.alpha();
        beta  = condition.beta();
        gamma = condition.gamma();
        if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
        {
            u20[0] = (gamma/alpha)*value;
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
        {
            u20[0]  = beta *(2.0*u10[0] - u00[0] + wvdsp_ht_05*u00[0]);
            u20[0] += beta *(2.0*p_aa_htht__hxhx*(u10[1]-u10[0]));
            u20[0] += beta *(ht_ht*f(sn,tn10));
            u20[0] += gamma*(-2.0*p_aa_htht__hx)*value;
            u20[0] /= beta*(1.0+wvdsp_ht_05);
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Robin)
        {
            u20[0]  = beta *(2.0*u10[0] - u00[0] + wvdsp_ht_05*u00[0]);
            u20[0] += beta *(2.0*p_aa_htht__hxhx*(u10[1]-u10[0]));
            u20[0] += beta *(ht_ht*f(sn,tn10));
            u20[0] += gamma*(-2.0*p_aa_htht__hx)*value;
            u20[0] += alpha*(+2.0*p_aa_htht__hx)*u10[0];
            u20[0] /= beta*(1.0+wvdsp_ht_05);
        }

        sn.i = static_cast<int>(N); sn.x = N*hx;
        value = boundary(sn, tn10, condition);
        alpha = condition.alpha();
        beta  = condition.beta();
        gamma = condition.gamma();
        if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
        {
            u20[N] = (gamma/alpha)*value;
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
        {
            u20[N]  = beta *(2.0*u10[N] - u00[N] + wvdsp_ht_05*u00[N]);
            u20[N] += beta *(2.0*p_aa_htht__hxhx*(u10[N-1]-u10[N]));
            u20[N] += beta *(ht_ht*f(sn,tn10));
            u20[N] += gamma*(+2.0*p_aa_htht__hx)*value;
            u20[N] /= beta*(1.0+wvdsp_ht_05);
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Robin)
        {
            u20[N]  = beta *(2.0*u10[N] - u00[N] + wvdsp_ht_05*u00[N]);
            u20[N] += beta *(2.0*p_aa_htht__hxhx*(u10[N-1]-u10[N]));
            u20[N] += beta *(ht_ht*f(sn,tn10));
            u20[N] += gamma*(+2.0*p_aa_htht__hx)*value;
            u20[N] += alpha*(-2.0*p_aa_htht__hx)*u10[N];
            u20[N] /= beta*(1.0+wvdsp_ht_05);
        }
        layerInfo(u20, tn20);

        for (unsigned int i=0; i<=N; i++)
        {
            u00[i] = u10[i];
            u10[i] = u20[i];
        }
    }
}

void IWaveEquationIBVP::implicit_calculate_D1V1() const
{
    const Dimension &dimX = spaceDimensionX();
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
    const double p_aa_htht__hx   = ((a*a)*(ht*ht))/hx;
    const double p_aa_htht__hxhx = ((a*a)*(ht*ht))/(hx*hx);

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
    /***********************************************************************************************/

    TimeNodePDE tn00; tn00.i = 0; tn00.t = tn00.i*ht;
    TimeNodePDE tn10; tn10.i = 1; tn10.t = tn10.i*ht;

    SpaceNodePDE sn;
    for (unsigned int n=0; n<=N; n++)
    {
        sn.i = static_cast<int>(n); sn.x = sn.i*hx;
        u00[n] = initial(sn, InitialCondition::InitialValue);
    }
    layerInfo(u00, tn00);

    /***********************************************************************************************/

    for (unsigned int n=0; n<=N; n++)
    {
        sn.i = static_cast<int>(n); sn.x = sn.i*hx;
        double firstDerivative = initial(sn, InitialCondition::FirstDerivative);
        u10[n] = u00[n] + firstDerivative*ht;
        double secndDerivative = 0.0;
        if (n==0) { secndDerivative = aa__hxhx*(+2.0*u00[0]   -5.0*u00[1]   +4.0*u00[2]   -1.0*u00[3]); }
        else if (n==N) { secndDerivative = aa__hxhx*(-1.0*u00[N-3] +4.0*u00[N-2] -5.0*u00[N-1] +2.0*u00[N]); }
        else { secndDerivative = aa__hxhx*(u00[n-1]-2.0*u00[n]+u00[n+1]); }
        u10[n] += (secndDerivative + f(sn,tn00) - alpha*firstDerivative) * htht_05;
    }
    layerInfo(u10, tn10);

    /***********************************************************************************************/
    /***********************************************************************************************/

    for (unsigned int ln=2; ln<=L; ln++)
    {
        TimeNodePDE tn00; tn00.i = ln-2; tn00.t = tn00.i*ht;
        TimeNodePDE tn10; tn10.i = ln-1; tn10.t = tn10.i*ht;
        TimeNodePDE tn20; tn20.i = ln;   tn20.t = tn20.i*ht;
        unsigned int s=0, e=N;
        BoundaryConditionPDE condition; double alpha, beta, gamma, value;

        for (unsigned int n=1; n<=N-1; n++)
        {
            sn.i = static_cast<int>(n); sn.x = n*hx;
            dx[n] = 0.0;
            dx[n] += (u10[n-1] - 2.0*u10[n] + u10[n+1])*p_aa_htht__hxhx_1m2lambda;
            dx[n] += (u00[n-1] - 2.0*u00[n] + u00[n+1])*p_aa_htht__hxhx_lambda;
            dx[n] += 2.0*u10[n] - u00[n] + alpha_ht_05*u00[n];
            dx[n] += ht_ht*f(sn, tn10);
            //dx[n] += ht_ht*(lmbd*f(sn, tn00) + m1_2lambda*f(sn, tn10) + lmbd*f(sn, tn20));
        }

        sn.i = static_cast<int>(0); sn.x = 0*hx;
        value = boundary(sn, tn20, condition);
        alpha = condition.alpha();
        beta  = condition.beta();
        gamma = condition.gamma();
        if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
        {
            s = 1;

            u20[0] = (gamma/alpha)*value;
            dx[1] -= u20[0]*m_aa_htht__hxhx_lambda;
            ax[1] = ax[0] = bx[0] = cx[0] = dx[0] = rx[0] = 0.0;
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
        {
            s = 0;

            ax[0]  = 0.0;
            bx[0]  = beta *(+1.0 + 2.0*p_aa_htht__hxhx*lmbd + alpha_ht_05);
            cx[0]  = beta *(-2.0*p_aa_htht__hxhx*lmbd);

            dx[0]  = beta *(2.0*u10[0] - u00[0] + alpha_ht_05*u00[0]);
            dx[0] += beta *(p_aa_htht__hxhx*(2.0*u10[0]-5.0*u10[1]+4.0*u10[2]-u10[3])*m1_2lambda);
            dx[0] += beta *(p_aa_htht__hxhx*(2.0*u00[0]-5.0*u00[1]+4.0*u00[2]-u00[3])*lmbd);
            dx[0] += beta *(ht_ht*f(sn,tn10));
            dx[0] += gamma*(-2.0*p_aa_htht__hx*lmbd)*value;
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Robin)
        {
            s = 0;

            ax[0]  = 0.0;
            bx[0]  = beta *(+1.0 + 2.0*p_aa_htht__hxhx*lmbd + alpha_ht_05);
            bx[0] += alpha*(-2.0*p_aa_htht__hx*lmbd);
            cx[0]  = beta *(-2.0*p_aa_htht__hxhx*lmbd);

            dx[0]  = beta *(2.0*u10[0] - u00[0] + alpha_ht_05*u00[0]);
            dx[0] += beta *(p_aa_htht__hxhx*(2.0*u10[0]-5.0*u10[1]+4.0*u10[2]-u10[3])*m1_2lambda);
            dx[0] += beta *(p_aa_htht__hxhx*(2.0*u00[0]-5.0*u00[1]+4.0*u00[2]-u00[3])*lmbd);
            dx[0] += beta *(ht_ht*f(sn,tn10));
            dx[0] += gamma*(-2.0*p_aa_htht__hx*lmbd)*value;
        }

        sn.i = static_cast<int>(N); sn.x = N*hx;
        value = boundary(sn, tn20, condition);
        alpha = condition.alpha();
        beta  = condition.beta();
        gamma = condition.gamma();
        if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
        {
            e = N-1;
            u20[N] = (gamma/alpha)*value;
            dx[N-1] -= u20[N]*m_aa_htht__hxhx_lambda;
            cx[N-1] = ax[N] = bx[N] = cx[N] = dx[N] = rx[N] = 0.0;
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
        {
            e = N;

            ax[N]  = beta *(-2.0*p_aa_htht__hxhx*lmbd);
            bx[N]  = beta *(+1.0 + 2.0*p_aa_htht__hxhx*lmbd + alpha_ht_05);
            cx[N]  = 0.0;

            dx[N]  = beta *(2.0*u10[N] - u00[N] + alpha_ht_05*u00[N]);
            dx[N] += beta *(p_aa_htht__hxhx*(-u10[N-3]+4.0*u10[N-2]-5.0*u10[N-1]+2.0*u10[N])*m1_2lambda);
            dx[N] += beta *(p_aa_htht__hxhx*(-u00[N-3]+4.0*u00[N-2]-5.0*u00[N-1]+2.0*u00[N])*lmbd);
            dx[N] += beta *(ht_ht*f(sn,tn10));
            dx[N] += gamma*(+2.0*p_aa_htht__hx*lmbd)*value;
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Robin)
        {
            e = N;

            ax[N]  = beta *(-2.0*p_aa_htht__hxhx*lmbd);
            bx[N]  = beta *(+1.0 + 2.0*p_aa_htht__hxhx*lmbd + alpha_ht_05);
            bx[N] += alpha*(+2.0*p_aa_htht__hx*lmbd);
            cx[N]  = 0.0;

            dx[N]  = beta *(2.0*u10[N] - u00[N] + alpha_ht_05*u00[N]);
            dx[N] += beta *(p_aa_htht__hxhx*(-u10[N-3]+4.0*u10[N-2]-5.0*u10[N-1]+2.0*u10[N])*m1_2lambda);
            dx[N] += beta *(p_aa_htht__hxhx*(-u00[N-3]+4.0*u00[N-2]-5.0*u00[N-1]+2.0*u00[N])*lmbd);
            dx[N] += beta *(ht_ht*f(sn,tn10));
            dx[N] += gamma*(+2.0*p_aa_htht__hx*lmbd)*value;
        }

        tomasAlgorithm(ax+s, bx+s, cx+s, dx+s, rx+s, e-s+1);
        for (unsigned int n=s; n<=e; n++) u20[n] = rx[n];
        layerInfo(u20, tn20);

        for (unsigned int i=0; i<=N; i++)
        {
            u00[i] = u10[i];
            u10[i] = u20[i];
        }
    }
}

void IWaveEquationIBVP::explicit_calculate_D2V1() const
{
    const unsigned int N = static_cast<unsigned int>( spaceDimensionX().size() );
    const unsigned int M = static_cast<unsigned int>( spaceDimensionY().size() );
    const unsigned int L = static_cast<unsigned int>( timeDimension().size() );

    const double hx = spaceDimensionX().step();
    const double hy = spaceDimensionY().step();
    const double ht = timeDimension().step();

    const double a = waveSpeed();
    const double alpha = waveDissipation();
    const double alpha_ht_05 = alpha*ht*0.5;
    const double ht_ht_05 = ht*ht*0.5;
    const double inv__1malpha_ht05 = 1.0/(1.0 + alpha_ht_05);
    const double aa_htht__hxhx = ((a*a*ht*ht)/(hx*hx));
    const double aa_htht__hyhy = ((a*a*ht*ht)/(hy*hy));
    const double aa__hxhx = (a*a)/(hx*hx);
    const double aa__hyhy = (a*a)/(hy*hy);

    DoubleMatrix u00(M+1, N+1);
    DoubleMatrix u10(M+1, N+1);
    DoubleMatrix u20(M+1, N+1);

    /***********************************************************************************************/
    /***********************************************************************************************/

    TimeNodePDE tn00; tn00.i = 0; tn00.t = tn00.i*ht;
    TimeNodePDE tn10; tn10.i = 1; tn10.t = tn10.i*ht;

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

    for (unsigned int m=0; m<=M; m++)
    {
        sn.j = static_cast<int>(m); sn.y = m*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = static_cast<int>(n); sn.x = n*hx;
            double firstDerivative = initial(sn, InitialCondition::FirstDerivative);

            double secndDerivative = 0.0;
            if (m==0) { secndDerivative += aa__hyhy*(+2.0*u00[0][n]-5.0*u00[1][n]+4.0*u00[2][n]-1.0*u00[3][n]); }
            else
                if (m==M) { secndDerivative += aa__hyhy*(-1.0*u00[M-3][n]+4.0*u00[M-2][n]-5.0*u00[M-1][n]+2.0*u00[M][n]); }
                else { secndDerivative += aa__hyhy*(u00[m-1][n]-2.0*u00[m][n]+u00[m+1][n]); }

            if (n==0) { secndDerivative += aa__hxhx*(+2.0*u00[m][0]-5.0*u00[m][1]+4.0*u00[m][2]-1.0*u00[m][3]); }
            else
                if (n==N) { secndDerivative += aa__hxhx*(-1.0*u00[m][N-3]+4.0*u00[m][N-2]-5.0*u00[m][N-1]+2.0*u00[m][N]); }
                else { secndDerivative += aa__hxhx*(u00[m][n-1]-2.0*u00[m][n]+u00[m][n+1]); }

            secndDerivative += f(sn,tn00);
            secndDerivative -=  alpha*firstDerivative;

            u10[m][n] = u00[m][n] + firstDerivative*ht + secndDerivative*ht_ht_05;
        }
    }
    layerInfo(u10, tn10);
    /***********************************************************************************************/

    for (unsigned int ln=2; ln<=L; ln++)
    {
        TimeNodePDE tn00; tn00.i = ln-2; tn00.t = tn00.i*ht;
        TimeNodePDE tn10; tn10.i = ln-1; tn10.t = tn10.i*ht;
        TimeNodePDE tn20; tn20.i = ln-0; tn20.t = tn20.i*ht;

        /**************************************************** border conditions ***************************************************/

        SpaceNodePDE sn0, sn1;
        BoundaryConditionPDE condition;

        sn0.i = static_cast<int>(0); sn0.x = 0*hx;
        sn1.i = static_cast<int>(N); sn1.x = N*hx;
        for (unsigned int m=0; m<=M; m++)
        {
            sn0.j = static_cast<int>(m); sn0.y = m*hy; u20[m][0] = boundary(sn0, tn20, condition);
            sn1.j = static_cast<int>(m); sn1.y = m*hy; u20[m][N] = boundary(sn1, tn20, condition);
        }

        sn0.j = static_cast<int>(0); sn0.y = 0*hy;
        sn1.j = static_cast<int>(M); sn1.y = M*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn0.i = static_cast<int>(n); sn0.x = n*hx; u20[0][n] = boundary(sn0, tn20, condition);
            sn1.i = static_cast<int>(n); sn1.x = n*hx; u20[M][n] = boundary(sn1, tn20, condition);
        }

        /**************************************************** border conditions ***************************************************/

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

    u00.clear();
    u10.clear();
    u20.clear();
}

void IWaveEquationIBVP::implicit_calculate_D2V1() const
{
    const unsigned int N = static_cast<unsigned int>( _spaceDimensionX.size() );
    const unsigned int M = static_cast<unsigned int>( _spaceDimensionY.size() );
    const unsigned int L = static_cast<unsigned int>( _timeDimension.size() );

    const double hx = _spaceDimensionX.step();
    const double hy = _spaceDimensionY.step();
    const double ht = _timeDimension.step();

    const double _lambda = lambda();
    const double a = waveSpeed();
    const double alpha = waveDissipation();
    const double ht_050 = ht*0.5;
    const double ht_ht_025 = ht*ht*0.25;
    const double ht_ht_025_05 = ht_ht_025*0.50;
    const double alpha_ht_025 = alpha*ht*0.25;
    const double aa__hxhx = (a*a)/(hx*hx);
    const double aa__hyhy = (a*a)/(hy*hy);

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
    /***********************************************************************************************/

    TimeNodePDE tn00; tn00.i = 0; tn00.t = 0.5*tn00.i*ht;
    TimeNodePDE tn05; tn05.i = 1; tn05.t = 0.5*tn05.i*ht;
    TimeNodePDE tn10; tn10.i = 2; tn10.t = 0.5*tn10.i*ht;

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

    for (unsigned int m=0; m<=M; m++)
    {
        sn.j = static_cast<int>(m); sn.y = sn.j*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = static_cast<int>(n); sn.x = sn.i*hx;
            double firstDerivative = initial(sn, InitialCondition::FirstDerivative);

            double secndDerivative = 0.0;
            if (m==0) { secndDerivative += aa__hyhy*(+2.0*u00[0][n]-5.0*u00[1][n]+4.0*u00[2][n]-1.0*u00[3][n]); }
            else
                if (m==M) { secndDerivative += aa__hyhy*(-1.0*u00[M-3][n]+4.0*u00[M-2][n]-5.0*u00[M-1][n]+2.0*u00[M][n]); }
                else { secndDerivative += aa__hyhy*(u00[m-1][n]-2.0*u00[m][n]+u00[m+1][n]); }

            if (n==0) { secndDerivative += aa__hxhx*(+2.0*u00[m][0]-5.0*u00[m][1]+4.0*u00[m][2]-1.0*u00[m][3]); }
            else
                if (n==N) { secndDerivative += aa__hxhx*(-1.0*u00[m][N-3]+4.0*u00[m][N-2]-5.0*u00[m][N-1]+2.0*u00[m][N]); }
                else { secndDerivative += aa__hxhx*(u00[m][n-1]-2.0*u00[m][n]+u00[m][n+1]); }

            secndDerivative += f(sn,tn00);
            secndDerivative -=  alpha*firstDerivative;

            u05[m][n] = u00[m][n] + firstDerivative*ht_050 + secndDerivative*ht_ht_025_05;
        }
    }
    layerInfo(u05, tn05);

    for (unsigned int m=0; m<=M; m++)
    {
        sn.j = static_cast<int>(m); sn.y = sn.j*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = static_cast<int>(n); sn.x = sn.i*hx;
            u10[m][n]  = 2.0*u05[m][n]-u00[m][n];

            double secndDerivative = 0.0;
            if (m==0) { secndDerivative += aa__hyhy*(+2.0*u05[0][n]-5.0*u05[1][n]+4.0*u05[2][n]-1.0*u05[3][n]); } else
                if (m==M) { secndDerivative += aa__hyhy*(-1.0*u05[M-3][n]+4.0*u05[M-2][n]-5.0*u05[M-1][n]+2.0*u05[M][n]); }
                else { secndDerivative += aa__hyhy*(u05[m-1][n]-2.0*u05[m][n]+u05[m+1][n]); }

            if (n==0) { secndDerivative += aa__hxhx*(+2.0*u05[m][0]-5.0*u05[m][1]+4.0*u05[m][2]-1.0*u05[m][3]); } else
                if (n==N) { secndDerivative += aa__hxhx*(-1.0*u05[m][N-3]+4.0*u05[m][N-2]-5.0*u05[m][N-1]+2.0*u05[m][N]); }
                else { secndDerivative += aa__hxhx*(u05[m][n-1]-2.0*u05[m][n]+u05[m][n+1]); }

            secndDerivative += f(sn,tn05);
            u10[m][n] += ht_ht_025*secndDerivative;
            u10[m][n] += alpha_ht_025*u00[m][n];

            u10[m][n] /= (1.0+alpha_ht_025);
        }
    }
    layerInfo(u10, tn10);
    /***********************************************************************************************/

    for (unsigned int ln=2; ln<=L; ln++)
    {
        TimeNodePDE tn05; tn05.i = 2*ln-3; tn05.t = 0.5*tn05.i*ht;
        TimeNodePDE tn10; tn10.i = 2*ln-2; tn10.t = 0.5*tn10.i*ht;
        TimeNodePDE tn15; tn15.i = 2*ln-1; tn15.t = 0.5*tn15.i*ht;
        TimeNodePDE tn20; tn20.i = 2*ln-0; tn20.t = 0.5*tn20.i*ht;

        /**************************************************** border conditions ***************************************************/

        SpaceNodePDE sn0, sn1;
        BoundaryConditionPDE condition;

        sn0.i = static_cast<int>(0); sn0.x = 0*hx;
        sn1.i = static_cast<int>(N); sn1.x = N*hx;
        for (unsigned int m=0; m<=M; m++)
        {
            sn0.j = static_cast<int>(m); sn0.y = m*hy; u15[m][0] = boundary(sn0, tn15, condition); u20[m][0] = boundary(sn0, tn20, condition);
            sn1.j = static_cast<int>(m); sn1.y = m*hy; u15[m][N] = boundary(sn1, tn15, condition); u20[m][N] = boundary(sn1, tn20, condition);
        }

        sn0.j = static_cast<int>(0); sn0.y = 0*hy;
        sn1.j = static_cast<int>(M); sn1.y = M*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn0.i = static_cast<int>(n); sn0.x = n*hx; u15[0][n] = boundary(sn0, tn15, condition); u20[0][n] = boundary(sn0, tn20, condition);
            sn1.i = static_cast<int>(n); sn1.x = n*hx; u15[M][n] = boundary(sn1, tn15, condition); u20[M][n] = boundary(sn1, tn20, condition);
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
                dx[n-1] += (u10[m][n-1] - 2.0*u10[m][n] + u10[m][n+1])*p_aa_htht__hxhx_025_1m2lambda;
                dx[n-1] += (u05[m][n-1] - 2.0*u05[m][n] + u05[m][n+1])*p_aa_htht__hxhx_025_lambda;
                dx[n-1] += (u10[m-1][n] - 2.0*u10[m][n] + u10[m+1][n])*p_aa_htht__hyhy_025;
                dx[n-1] += 2.0*u10[m][n] - u05[m][n] + alpha_ht_025*u05[m][n];
                dx[n-1] += ht_ht_025*f(sn, tn10);
                //dx[n-1] += ht_ht_025*(_lambda*f(sn, tn05)+(1.0-2.0*_lambda)*f(sn, tn10)+_lambda*f(sn, tn15));
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
                //dy[m-1] += ht_ht_025*(_lambda*f(sn, tn10)+(1.0-2.0*_lambda)*f(sn, tn15)+_lambda*f(sn, tn20));
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

double IWaveEquationIBVP::lambda() const { return 0.25; }

//--------------------------------------------------------------------------------------------------------------//

IFinalWaveEquationIBVP::IFinalWaveEquationIBVP(double waveSpeed, double waveDissipation)
    : _waveSpeed(waveSpeed), _waveDissipation(waveDissipation) {}

IFinalWaveEquationIBVP::IFinalWaveEquationIBVP(const IFinalWaveEquationIBVP &ibvp) :
    _waveSpeed(ibvp._waveSpeed), _waveDissipation(ibvp._waveDissipation) {}

IFinalWaveEquationIBVP & IFinalWaveEquationIBVP::operator =(const IFinalWaveEquationIBVP &other)
{
    this->_waveSpeed = other._waveSpeed;
    this->_waveDissipation = other._waveDissipation;
    return *this;
}

IFinalWaveEquationIBVP::~IFinalWaveEquationIBVP() {}

double IFinalWaveEquationIBVP::waveSpeed() const  { return _waveSpeed; }

double IFinalWaveEquationIBVP::waveDissipation() const { return _waveDissipation; }

void IFinalWaveEquationIBVP::setWaveSpeed(double waveSpeed) { _waveSpeed = waveSpeed; }

void IFinalWaveEquationIBVP::setWaveDissipation(double waveDissipation) { _waveDissipation = waveDissipation; }

void IFinalWaveEquationIBVP::explicit_calculate_D1V1() const {}

void IFinalWaveEquationIBVP::implicit_calculate_D1V1() const {}

void IFinalWaveEquationIBVP::explicit_calculate_D2V1() const
{
    const unsigned int N = static_cast<unsigned int>( _spaceDimensionX.size() );
    const unsigned int M = static_cast<unsigned int>( _spaceDimensionY.size() );
    const unsigned int L = static_cast<unsigned int>( _timeDimension.size() );

    const double hx = _spaceDimensionX.step();
    const double hy = _spaceDimensionY.step();
    const double ht = _timeDimension.step();

    const double a = waveSpeed();
    const double alpha = waveDissipation();
    const double alpha_ht_05 = alpha*ht*0.5;
    const double ht_ht_05 = ht*ht*0.5;
    const double inv__1malpha_ht05 = 1.0/(1.0 + alpha_ht_05);
    const double aa_htht__hxhx = ((a*a*ht*ht)/(hx*hx));
    const double aa_htht__hyhy = ((a*a*ht*ht)/(hy*hy));
    const double aa__hxhx = (a*a)/(hx*hx);
    const double aa__hyhy = (a*a)/(hy*hy);


    DoubleMatrix p00(M+1, N+1);
    DoubleMatrix p10(M+1, N+1);
    DoubleMatrix p20(M+1, N+1);

    /***********************************************************************************************/
    /***********************************************************************************************/

    TimeNodePDE tn00; tn00.i = L-0; tn00.t = tn00.i*ht;
    TimeNodePDE tn10; tn10.i = L-1; tn10.t = tn10.i*ht;

    SpaceNodePDE sn;
    for (unsigned int m=0; m<=M; m++)
    {
        sn.j = static_cast<int>(m); sn.y = sn.j*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = static_cast<int>(n); sn.x = sn.i*hx;
            p00[m][n] = final(sn, FinalCondition::FinalValue);
        }
    }
    layerInfo(p00, tn00);

    /***********************************************************************************************/

    for (unsigned int m=0; m<=M; m++)
    {
        sn.j = static_cast<int>(m); sn.y = m*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = static_cast<int>(n); sn.x = n*hx;
            double firstDerivative = final(sn, FinalCondition::FinalFirstDerivative);

            double secndDerivative = 0.0;
            if (m==0) { secndDerivative += aa__hyhy*(+2.0*p00[0][n]-5.0*p00[1][n]+4.0*p00[2][n]-1.0*p00[3][n]); }
            else
                if (m==M) { secndDerivative += aa__hyhy*(-1.0*p00[M-3][n]+4.0*p00[M-2][n]-5.0*p00[M-1][n]+2.0*p00[M][n]); }
                else { secndDerivative += aa__hyhy*(p00[m-1][n]-2.0*p00[m][n]+p00[m+1][n]); }

            if (n==0) { secndDerivative += aa__hxhx*(+2.0*p00[m][0]-5.0*p00[m][1]+4.0*p00[m][2]-1.0*p00[m][3]); }
            else
                if (n==N) { secndDerivative += aa__hxhx*(-1.0*p00[m][N-3]+4.0*p00[m][N-2]-5.0*p00[m][N-1]+2.0*p00[m][N]); }
                else { secndDerivative += aa__hxhx*(p00[m][n-1]-2.0*p00[m][n]+p00[m][n+1]); }

            secndDerivative += f(sn,tn00);
            secndDerivative +=  alpha*firstDerivative;

            p10[m][n] = p00[m][n] - firstDerivative*ht + secndDerivative*ht_ht_05;
        }
    }
    layerInfo(p10, tn10);
    /***********************************************************************************************/
    const unsigned int size_ln = static_cast<unsigned int>(0)-1;
    for (unsigned int ln=L-2; ln != size_ln; ln--)
    {
        TimeNodePDE tn00; tn00.i = ln+2; tn00.t = tn00.i*ht;
        TimeNodePDE tn10; tn10.i = ln+1; tn10.t = tn10.i*ht;
        TimeNodePDE tn20; tn20.i = ln+0; tn20.t = tn20.i*ht;

        /**************************************************** border conditions ***************************************************/

        SpaceNodePDE sn0, sn1;
        BoundaryConditionPDE condition;

        sn0.i = static_cast<int>(0); sn0.x = 0*hx;
        sn1.i = static_cast<int>(N); sn1.x = N*hx;
        for (unsigned int m=0; m<=M; m++)
        {
            sn0.j = static_cast<int>(m); sn0.y = m*hy; p20[m][0] = boundary(sn0, tn20, condition);
            sn1.j = static_cast<int>(m); sn1.y = m*hy; p20[m][N] = boundary(sn1, tn20, condition);
        }

        sn0.j = static_cast<int>(0); sn0.y = 0*hy;
        sn1.j = static_cast<int>(M); sn1.y = M*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn0.i = static_cast<int>(n); sn0.x = n*hx; p20[0][n] = boundary(sn0, tn20, condition);
            sn1.i = static_cast<int>(n); sn1.x = n*hx; p20[M][n] = boundary(sn1, tn20, condition);
        }

        /**************************************************** border conditions ***************************************************/

        SpaceNodePDE sn;
        for (unsigned int m=1; m<=M-1; m++)
        {
            sn.j = static_cast<int>(m); sn.y = m*hy;
            for (unsigned int n=1; n<=N-1; n++)
            {
                sn.i = static_cast<int>(n); sn.x = n*hx;
                p20[m][n] = (aa_htht__hxhx*(p10[m][n-1]-2.0*p10[m][n]+p10[m][n+1])
                        + aa_htht__hyhy*(p10[m-1][n]-2.0*p10[m][n]+p10[m+1][n])
                        + 2.0*p10[m][n] - p00[m][n] + alpha_ht_05*p00[m][n]
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

    p00.clear();
    p10.clear();
    p20.clear();
}

void IFinalWaveEquationIBVP::implicit_calculate_D2V1() const
{
    const unsigned int N = static_cast<unsigned int>( _spaceDimensionX.size() );
    const unsigned int M = static_cast<unsigned int>( _spaceDimensionY.size() );
    const unsigned int L = static_cast<unsigned int>( _timeDimension.size() );

    const double hx = _spaceDimensionX.step();
    const double hy = _spaceDimensionY.step();
    const double ht = _timeDimension.step();

    const double _lambda = lambda();
    const double a = waveSpeed();
    const double alpha = waveDissipation();
    const double ht_050 = ht*0.50;
    const double ht_ht_025 = ht*ht*0.25;
    const double ht_ht_025_05 = ht_ht_025*0.50;
    const double alpha_ht_025 = alpha*ht*0.25;
    const double aa__hxhx = (a*a)/(hx*hx);
    const double aa__hyhy = (a*a)/(hy*hy);

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

    DoubleMatrix p00(M+1, N+1);
    DoubleMatrix p05(M+1, N+1);
    DoubleMatrix p10(M+1, N+1);
    DoubleMatrix p15(M+1, N+1);
    DoubleMatrix p20(M+1, N+1);

    /***********************************************************************************************/
    /***********************************************************************************************/

    TimeNodePDE tn00; tn00.i = 2*L-0; tn00.t = 0.5*tn00.i*ht;
    TimeNodePDE tn05; tn05.i = 2*L-1; tn05.t = 0.5*tn05.i*ht;
    TimeNodePDE tn10; tn10.i = 2*L-2; tn10.t = 0.5*tn10.i*ht;

    SpaceNodePDE sn;
    for (unsigned int m=0; m<=M; m++)
    {
        sn.j = static_cast<int>(m); sn.y = sn.j*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = static_cast<int>(n); sn.x = sn.i*hx;
            p00[m][n] = final(sn, FinalCondition::FinalValue);
        }
    }
    layerInfo(p00, tn00);

    /***********************************************************************************************/

    for (unsigned int m=0; m<=M; m++)
    {
        sn.j = static_cast<int>(m); sn.y = sn.j*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = static_cast<int>(n); sn.x = sn.i*hx;
            double firstDerivative = final(sn, FinalCondition::FinalFirstDerivative);

            double secndDerivative = 0.0;
            if (m==0) { secndDerivative += aa__hyhy*(+2.0*p00[0][n]-5.0*p00[1][n]+4.0*p00[2][n]-1.0*p00[3][n]); }
            else
                if (m==M) { secndDerivative += aa__hyhy*(-1.0*p00[M-3][n]+4.0*p00[M-2][n]-5.0*p00[M-1][n]+2.0*p00[M][n]); }
                else { secndDerivative += aa__hyhy*(p00[m-1][n]-2.0*p00[m][n]+p00[m+1][n]); }

            if (n==0) { secndDerivative += aa__hxhx*(+2.0*p00[m][0]-5.0*p00[m][1]+4.0*p00[m][2]-1.0*p00[m][3]); }
            else
                if (n==N) { secndDerivative += aa__hxhx*(-1.0*p00[m][N-3]+4.0*p00[m][N-2]-5.0*p00[m][N-1]+2.0*p00[m][N]); }
                else { secndDerivative += aa__hxhx*(p00[m][n-1]-2.0*p00[m][n]+p00[m][n+1]); }

            secndDerivative += f(sn,tn00);
            secndDerivative += alpha*firstDerivative;

            p05[m][n] = p00[m][n] - firstDerivative*ht_050 + secndDerivative*ht_ht_025_05;
            //p10[m][n] = p00[m][n] - firstDerivative*ht    + secndDericative*ht_ht_05;
        }
    }
    layerInfo(p05, tn05);

    for (unsigned int m=0; m<=M; m++)
    {
        sn.j = static_cast<int>(m); sn.y = sn.j*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = static_cast<int>(n); sn.x = sn.i*hx;
            p10[m][n]  = 2.0*p05[m][n]-p00[m][n];

            double secndDerivative = 0.0;
            if (m==0) { secndDerivative += aa__hyhy*(+2.0*p05[0][n]-5.0*p05[1][n]+4.0*p05[2][n]-1.0*p05[3][n]); } else
                if (m==M) { secndDerivative += aa__hyhy*(-1.0*p05[M-3][n]+4.0*p05[M-2][n]-5.0*p05[M-1][n]+2.0*p05[M][n]); }
                else { secndDerivative += aa__hyhy*(p05[m-1][n]-2.0*p05[m][n]+p05[m+1][n]); }

            if (n==0) { secndDerivative += aa__hxhx*(+2.0*p05[m][0]-5.0*p05[m][1]+4.0*p05[m][2]-1.0*p05[m][3]); } else
                if (n==N) { secndDerivative += aa__hxhx*(-1.0*p05[m][N-3]+4.0*p05[m][N-2]-5.0*p05[m][N-1]+2.0*p05[m][N]); }
                else { secndDerivative += aa__hxhx*(p05[m][n-1]-2.0*p05[m][n]+p05[m][n+1]); }

            secndDerivative += f(sn,tn05);
            p10[m][n] += ht_ht_025*secndDerivative;
            p10[m][n] += alpha_ht_025*p00[m][n];

            p10[m][n] /= (1.0+alpha_ht_025);
        }
    }
    layerInfo(p10, tn10);
    /***********************************************************************************************/
    const unsigned int size_ln = static_cast<unsigned int>(0)-1;
    for (unsigned int ln=L-2; ln !=size_ln; ln--)
    {
        TimeNodePDE tn05; tn05.i = 2*ln+3; tn05.t = 0.5*tn05.i*ht;
        TimeNodePDE tn10; tn10.i = 2*ln+2; tn10.t = 0.5*tn10.i*ht;
        TimeNodePDE tn15; tn15.i = 2*ln+1; tn15.t = 0.5*tn15.i*ht;
        TimeNodePDE tn20; tn20.i = 2*ln+0; tn20.t = 0.5*tn20.i*ht;

        /**************************************************** border conditions ***************************************************/

        SpaceNodePDE sn0, sn1;
        BoundaryConditionPDE condition;

        sn0.i = static_cast<int>(0); sn0.x = 0*hx;
        sn1.i = static_cast<int>(N); sn1.x = N*hx;
        for (unsigned int m=0; m<=M; m++)
        {
            sn0.j = static_cast<int>(m); sn0.y = m*hy; p15[m][0] = boundary(sn0, tn15, condition); p20[m][0] = boundary(sn0, tn20, condition);
            sn1.j = static_cast<int>(m); sn1.y = m*hy; p15[m][N] = boundary(sn1, tn15, condition); p20[m][N] = boundary(sn1, tn20, condition);
        }

        sn0.j = static_cast<int>(0); sn0.y = 0*hy;
        sn1.j = static_cast<int>(M); sn1.y = M*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn0.i = static_cast<int>(n); sn0.x = n*hx; p15[0][n] = boundary(sn0, tn15, condition); p20[0][n] = boundary(sn0, tn20, condition);
            sn1.i = static_cast<int>(n); sn1.x = n*hx; p15[M][n] = boundary(sn1, tn15, condition); p20[M][n] = boundary(sn1, tn20, condition);
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
                dx[n-1] += (p10[m][n-1] - 2.0*p10[m][n] + p10[m][n+1])*p_aa_htht__hxhx_025_1m2lambda;
                dx[n-1] += (p05[m][n-1] - 2.0*p05[m][n] + p05[m][n+1])*p_aa_htht__hxhx_025_lambda;
                dx[n-1] += (p10[m-1][n] - 2.0*p10[m][n] + p10[m+1][n])*p_aa_htht__hyhy_025;
                dx[n-1] += 2.0*p10[m][n] - p05[m][n] + alpha_ht_025*p05[m][n];
                dx[n-1] += ht_ht_025*f(sn, tn10);
                //dx[n-1] += ht_ht_025*(_lambda*f(sn, tn05)+(1.0-2.0*_lambda)*f(sn, tn10)+_lambda*f(sn, tn15));
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
                dy[m-1] += ht_ht_025*f(sn, tn15);
                //dy[m-1] += ht_ht_025*(_lambda*f(sn, tn10)+(1.0-2.0*_lambda)*f(sn, tn15)+_lambda*f(sn, tn20));
            }
            dy[0]   -= p20[0][n]*m_aa_htht__hyhy_025_lambda;
            dy[M-2] -= p20[M][n]*m_aa_htht__hyhy_025_lambda;
            tomasAlgorithm(ay, by, cy, dy, ry, M-1);
            for (unsigned int m=1; m<=M-1; m++) p20[m][n] = ry[m-1];
        }
        layerInfo(p20, tn20);
        /**************************************************** y direction apprx ***************************************************/

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

double IFinalWaveEquationIBVP::lambda() const { return 0.25; }

//--------------------------------------------------------------------------------------------------------------//
