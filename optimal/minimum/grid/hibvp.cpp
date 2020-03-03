#include "hibvp.h"

void IHyperbolicIBVP::layerInfo(const DoubleVector &, const TimeNodePDE &) const {}

void IHyperbolicIBVP::layerInfo(const DoubleMatrix &, const TimeNodePDE &) const {}

//--------------------------------------------------------------------------------------------------------------//

void IHyperbolicFBVP::layerInfo(const DoubleVector &, const TimeNodePDE &) const {}

void IHyperbolicFBVP::layerInfo(const DoubleMatrix &, const TimeNodePDE &) const {}

//--------------------------------------------------------------------------------------------------------------//

IWaveEquationIBVP::IWaveEquationIBVP(double waveSpeed, double waveDissipation, double unknownC, double unknownD) : IHyperbolicIBVP(),
    _waveSpeed(waveSpeed), _waveDissipation(waveDissipation), _unknownB(unknownC), _restoration(unknownD) {}

IWaveEquationIBVP::IWaveEquationIBVP(const IWaveEquationIBVP& other) : IHyperbolicIBVP(other),
    _waveSpeed(other._waveSpeed), _waveDissipation(other._waveDissipation), _unknownB(other._unknownB), _restoration(other._restoration) {}

IWaveEquationIBVP& IWaveEquationIBVP::operator=(const IWaveEquationIBVP &other)
{
    if (this == &other) { return *this; }

    this->_waveSpeed = other._waveSpeed;
    this->_waveDissipation = other._waveDissipation;
    this->_unknownB = other._unknownB;
    this->_restoration = other._restoration;
    return *this;
}

IWaveEquationIBVP::~IWaveEquationIBVP() {}

double IWaveEquationIBVP::waveSpeed() const { return _waveSpeed; }

void IWaveEquationIBVP::setWaveSpeed(double waveSpeed) { this->_waveSpeed = waveSpeed; }

double IWaveEquationIBVP::waveDissipation() const { return _waveDissipation; }

void IWaveEquationIBVP::setWaveDissipation(double waveDissipation) { this->_waveDissipation = waveDissipation; }

void IWaveEquationIBVP::setUnknownB(double unknownB) { this->_unknownB = unknownB; }

double IWaveEquationIBVP::unknownB() const { return _unknownB; }

void IWaveEquationIBVP::setRestoration(double restoration) { this->_restoration = restoration; }

double IWaveEquationIBVP::restoration() const { return _restoration; }

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
        double firstDerivative = initial(sn, InitialCondition::InitialFirstDerivative);
        u10[n] = u00[n] + ht*firstDerivative;
        double secndDerivative = 0.0;
        if (n==0)      { secndDerivative = aa__hxhx*(+2.0*u00[0] - 5.0*u00[1] + 4.0*u00[2] - 1.0*u00[3]); }
        else if (n==N) { secndDerivative = aa__hxhx*(-1.0*u00[N-3] + 4.0*u00[N-2] - 5.0*u00[N-1] + 2.0*u00[N]); }
        else           { secndDerivative = aa__hxhx*(u00[n-1] - 2.0*u00[n] + u00[n+1]); }
        u10[n] += (secndDerivative + f(sn,tn00) - wd*firstDerivative) * htht_05;
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
    const unsigned int N = static_cast<unsigned int>(spaceDimensionX().size()) - 1;
    const unsigned int L = static_cast<unsigned int>(timeDimension().size()) - 1;

    const int xmin = spaceDimensionX().min();
    const int xmax = spaceDimensionX().max();
    //const int tmin = timeDimension().min();
    //const int tmax = timeDimension().max();

    const double hx = spaceDimensionX().step();
    const double ht = timeDimension().step();

    const double a = waveSpeed();
    const double b = unknownB();
    const double c = restoration();
    const double d = waveDissipation();
    const double w = weight();
    const double w1 = w;
    const double w2 = 1.0 - 2.0*w;
    const double w3 = w;

    // equation parameters
    const double k11 = -(a*a)*((ht*ht)/(hx*hx))*w1 + ((b*ht*ht)/(2.0*hx))*w1;
    const double k12 = +1.0 + 2.0*(a*a)*((ht*ht)/(hx*hx))*w1 - c*ht*ht*w1 + 0.5*d*ht;
    const double k13 = -(a*a)*((ht*ht)/(hx*hx))*w1 - ((b*ht*ht)/(2.0*hx))*w1;
    const double k21 = +(a*a)*((ht*ht)/(hx*hx))*w2 - ((b*ht*ht)/(2.0*hx))*w2;
    const double k22 = +2.0 - 2.0*(a*a)*((ht*ht)/(hx*hx))*w2 + c*ht*ht*w2;
    const double k23 = +(a*a)*((ht*ht)/(hx*hx))*w2 + ((b*ht*ht)/(2.0*hx))*w2;
    const double k31 = +(a*a)*((ht*ht)/(hx*hx))*w3 - ((b*ht*ht)/(2.0*hx))*w3;
    const double k32 = -1.0 - 2.0*(a*a)*((ht*ht)/(hx*hx))*w3 + c*ht*ht*w3 + 0.5*d*ht;
    const double k33 = +(a*a)*((ht*ht)/(hx*hx))*w3 + ((b*ht*ht)/(2.0*hx))*w3;

    // left border condition parameters
    const double b11 = +1.0 + 2.0*((a*a*ht*ht)/(hx*hx))*w1 - c*ht*ht*w1 + 0.5*d*ht;
    const double b12 = -2.0*((a*a*ht*ht)/hx)*w1 + b*ht*ht*w1;
    const double b13 = -2.0*((a*a*ht*ht)/(hx*hx))*w1;
    const double b14 = +2.0 - 2.0*((a*a*ht*ht)/(hx*hx))*w2 + c*ht*ht*w2;
    const double b15 = +2.0*((a*a*ht*ht)/hx)*w2 - b*ht*ht*w2;
    const double b16 = +2.0*((a*a*ht*ht)/(hx*hx))*w2;
    const double b17 = -1.0 - 2.0*((a*a*ht*ht)/(hx*hx))*w3 + c*ht*ht*w3 + 0.5*d*ht;
    const double b18 = +2.0*((a*a*ht*ht)/hx)*w3 - b*ht*ht*w3;
    const double b19 = +2.0*((a*a*ht*ht)/(hx*hx))*w3;
    const double b20 = -2.0*((a*a*ht*ht)/hx) + b*ht*ht;

    // right border condition parameters
    const double b21 = -2.0*((a*a*ht*ht)/(hx*hx))*w1;
    const double b22 = +2.0*((a*a*ht*ht)/hx)*w1 + b*ht*ht*w1;
    const double b23 = +1.0 + 2.0*((a*a*ht*ht)/(hx*hx))*w1 - c*ht*ht*w1 + 0.5*d*ht;
    const double b24 = +2.0*((a*a*ht*ht)/(hx*hx))*w2;
    const double b25 = -2.0*((a*a*ht*ht)/hx)*w2 - b*ht*ht*w2;
    const double b26 = +2.0 - 2.0*((a*a*ht*ht)/(hx*hx))*w2 + c*ht*ht*w2;
    const double b27 = +2.0*((a*a*ht*ht)/(hx*hx))*w3;
    const double b28 = -2.0*((a*a*ht*ht)/hx)*w3 - b*ht*ht*w3;
    const double b29 = -1.0 - 2.0*((a*a*ht*ht)/(hx*hx))*w3 + c*ht*ht*w3 + 0.5*d*ht;
    const double b30 = +2.0*((a*a*ht*ht)/hx) + b*ht*ht;

    // initial condition parameters
    const double aa__hxhx = (a*a)/(hx*hx);
    const double htht_05  = 0.5*ht*ht;
    const double ht_ht    = ht*ht;

    double *ax = static_cast<double*>(malloc(sizeof(double)*(N+1)));
    double *bx = static_cast<double*>(malloc(sizeof(double)*(N+1)));
    double *cx = static_cast<double*>(malloc(sizeof(double)*(N+1)));
    double *dx = static_cast<double*>(malloc(sizeof(double)*(N+1)));
    double *rx = static_cast<double*>(malloc(sizeof(double)*(N+1)));

    for (unsigned int n=0; n<=N; n++)
    {
        ax[n] = k11;
        bx[n] = k12;
        cx[n] = k13;
    }
    ax[0] = 0.0; cx[N] = 0.0;

    DoubleVector u00(N+1);
    DoubleVector u10(N+1);
    DoubleVector u20(N+1);

    /***********************************************************************************************/
    /***********************************************************************************************/

    TimeNodePDE tn00; tn00.i = 0; tn00.t = tn00.i*ht;
    TimeNodePDE tn10; tn10.i = 1; tn10.t = tn10.i*ht;

    SpaceNodePDE sn;
    unsigned int i = 0;
    for (int n=xmin; n<=xmax; n++, i++)
    {
        sn.i = n; sn.x = n*hx;
        u00[i] = initial(sn, InitialCondition::InitialValue);
    }
    layerInfo(u00, tn00);

    /***********************************************************************************************/

    i = 0;
    for (int n=xmin; n<=xmax; n++, i++)
    {
        sn.i = n; sn.x = n*hx;
        double firstDerivativeT = initial(sn, InitialCondition::InitialFirstDerivative);
        u10[i] = u00[i] + firstDerivativeT*ht;
        double secndDerivativeX = 0.0;
        double firstDerivativeX = 0.0;
        if (n==xmin)
        {
            secndDerivativeX = aa__hxhx*(+2.0*u00[0]-5.0*u00[1]+4.0*u00[2]-1.0*u00[3]);
            firstDerivativeX = (-3.0*u00[0]+4.0*u00[1]-u00[2])/(2.0*hx);
        }
        else if (n==xmax)
        {
            secndDerivativeX = aa__hxhx*(-1.0*u00[N-3]+4.0*u00[N-2] -5.0*u00[N-1] +2.0*u00[N]);
            firstDerivativeX = (u00[N-2]-4.0*u00[N-1]+3.0*u00[N])/(2.0*hx);
        }
        else
        {
            secndDerivativeX = aa__hxhx*(u00[i-1]-2.0*u00[i]+u00[i+1]);
            firstDerivativeX = (u00[i+1]-u00[i-1])/(2.0*hx);
        }
        u10[i] += (secndDerivativeX + b*firstDerivativeX + c*u00[i] + f(sn,tn00) - d*firstDerivativeT) * htht_05;
    }
    layerInfo(u10, tn10);

    /***********************************************************************************************/
    /***********************************************************************************************/

    for (unsigned int ln=2; ln<=L; ln++)
    {
        TimeNodePDE tn00; tn00.i = ln-2; tn00.t = tn00.i*ht;
        TimeNodePDE tn10; tn10.i = ln-1; tn10.t = tn10.i*ht;
        TimeNodePDE tn20; tn20.i = ln;   tn20.t = tn20.i*ht;

        unsigned int i = 1;
        for (int n=xmin+1; n<=xmax-1; n++, i++)
        {
            sn.i = n; sn.x = n*hx;
            dx[i] = 0.0;
            dx[i] += k21 * u10[i-1] + k22 * u10[i] + k23 * u10[i+1];
            dx[i] += k31 * u00[i-1] + k32 * u00[i] + k33 * u00[i+1];

            dx[i] += ht_ht*f(sn, tn10);
            //dx[i] += ht_ht*w1*f(sn, tn20) + ht_ht*w2*f(sn, tn10) + ht_ht*w3*f(sn, tn00);
        }

        unsigned int s=0, e=N;
        BoundaryConditionPDE condition; double alpha, beta, gamma, value;

        sn.i = xmin; sn.x = xmin*hx;
        value = boundary(sn, tn20, condition);
        alpha = condition.alpha();
        beta  = condition.beta();
        gamma = condition.gamma();
        if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
        {
            s = 1;

            u20[0] = (gamma/alpha)*value;
            dx[1] -= k11 * u20[0];
            ax[1] = ax[0] = bx[0] = cx[0] = dx[0] = rx[0] = 0.0;
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
        {
            s = 0;

            ax[s]  = 0.0;
            bx[s]  = beta * b11 + alpha * b12;
            cx[s]  = beta * b13;

            dx[s]  = beta * b14 * u10[s] + beta * b16 * u10[s+1];
            dx[s] += beta * b17 * u00[s] + beta * b19 * u00[s+1];

            dx[s] += gamma * b20 * boundary(sn, tn10, condition);
            dx[s] += beta  * ht_ht*f(sn,tn10);

            //dx[s] += gamma * b20 * (w1*value + w2*boundary(sn, tn10, condition) + w3*boundary(sn, tn00, condition));
            //dx[s] += beta * (ht_ht*w1*f(sn, tn20) + ht_ht*w2*f(sn, tn10) + ht_ht*w3*f(sn, tn00));
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Robin)
        {
            s = 0;

            //ax[s] = 0.0;
            //bx[s] = hx*alpha-beta;
            //cx[s] = beta;
            //dx[s] = hx*gamma*value;

            ax[s]  = 0.0;
            bx[s]  = beta * b11 + alpha * b12;
            cx[s]  = beta * b13;

            dx[s]  = beta * b14 * u10[s] + alpha * b15 * u10[s] + beta * b16 * u10[s+1];
            dx[s] += beta * b17 * u00[s] + alpha * b18 * u00[s] + beta * b19 * u00[s+1];

            dx[s] += gamma * b20 * boundary(sn, tn10, condition);
            dx[s] += beta  * ht_ht*f(sn,tn10);

            //dx[s] += gamma * b20 * (w1*value + w2*boundary(sn, tn10, condition) + w3*boundary(sn, tn00, condition));
            //dx[s] += beta * (ht_ht*w1*f(sn, tn20) + ht_ht*w2*f(sn, tn10) + ht_ht*w3*f(sn, tn00));
        }

        sn.i = xmax; sn.x = xmax*hx;
        value = boundary(sn, tn20, condition);
        alpha = condition.alpha();
        beta  = condition.beta();
        gamma = condition.gamma();
        if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
        {
            e = N-1;

            u20[N] = (gamma/alpha)*value;
            dx[N-1] -= k13 * u20[N];
            cx[N-1] = ax[N] = bx[N] = cx[N] = dx[N] = rx[N] = 0.0;
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
        {
            e = N;

            ax[e]  = beta * b21;
            bx[e]  = beta * b23 + alpha * b22;
            cx[e]  = 0.0;

            dx[e]  = beta * b24 * u10[e-1] + beta * b26 * u10[e];
            dx[e] += beta * b27 * u00[e-1] + beta * b29 * u00[e];

            dx[e] += gamma * b30 * boundary(sn, tn10, condition);
            dx[e] += beta  * ht_ht*f(sn,tn10);

            //dx[e] += gamma * b30 * (w1*value + w2*boundary(sn, tn10, condition) + w3*boundary(sn, tn00, condition));
            //dx[e] += beta * (ht_ht*w1*f(sn, tn20) + ht_ht*w2*f(sn, tn10) + ht_ht*w3*f(sn, tn00));
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Robin)
        {
            e = N;

            //ax[e] = -beta;
            //bx[e] = hx*alpha+beta;
            //cx[e] = 0.0;
            //dx[e] = hx*gamma*value;

            ax[e]  = beta * b21;
            bx[e]  = beta * b23 + alpha * b22;
            cx[e]  = 0.0;

            dx[e]  = beta * b24 * u10[e-1] + alpha * b25 * u10[e] + beta * b26 * u10[e];
            dx[e] += beta * b27 * u00[e-1] + alpha * b28 * u00[e] + beta * b29 * u00[e];

            dx[e] += gamma * b30 * boundary(sn, tn10, condition);
            dx[e] += beta  * ht_ht*f(sn,tn10);

            //dx[e] += gamma * b30 * (w1*value + w2*boundary(sn, tn10, condition) + w3*boundary(sn, tn00, condition));
            //dx[e] += beta * (ht_ht*w1*f(sn, tn20) + ht_ht*w2*f(sn, tn10) + ht_ht*w3*f(sn, tn00));
        }

        tomasAlgorithmLeft2Right(ax+s, bx+s, cx+s, dx+s, rx+s, e-s+1);
        for (unsigned int n=s; n<=e; n++) u20[n] = rx[n];
        layerInfo(u20, tn20);

        for (unsigned int i=0; i<=N; i++)
        {
            u00[i] = u10[i];
            u10[i] = u20[i];
        }
    }

    u00.clear();
    u10.clear();
    u20.clear();

    free(rx);
    free(dx);
    free(cx);
    free(bx);
    free(ax);
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
            double firstDerivative = initial(sn, InitialCondition::InitialFirstDerivative);

            double secndDerivative = 0.0;
            if (m==0)      { secndDerivative += aa__hyhy*(+2.0*u00[0][n]-5.0*u00[1][n]+4.0*u00[2][n]-1.0*u00[3][n]); }
            else if (m==M) { secndDerivative += aa__hyhy*(-1.0*u00[M-3][n]+4.0*u00[M-2][n]-5.0*u00[M-1][n]+2.0*u00[M][n]); }
            else           { secndDerivative += aa__hyhy*(u00[m-1][n]-2.0*u00[m][n]+u00[m+1][n]); }

            if (n==0)      { secndDerivative += aa__hxhx*(+2.0*u00[m][0]-5.0*u00[m][1]+4.0*u00[m][2]-1.0*u00[m][3]); }
            else if (n==N) { secndDerivative += aa__hxhx*(-1.0*u00[m][N-3]+4.0*u00[m][N-2]-5.0*u00[m][N-1]+2.0*u00[m][N]); }
            else           { secndDerivative += aa__hxhx*(u00[m][n-1]-2.0*u00[m][n]+u00[m][n+1]); }

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
    const unsigned int N = static_cast<unsigned int>(spaceDimensionX().size()) -1;
    const unsigned int M = static_cast<unsigned int>(spaceDimensionY().size()) -1;
    const unsigned int L = static_cast<unsigned int>(timeDimension().size()) -1;

    const int xmin = spaceDimensionX().min();
    const int xmax = spaceDimensionX().max();
    const int ymin = spaceDimensionY().min();
    const int ymax = spaceDimensionY().max();
    //const int tmin = timeDimension().min();
    //const int tmax = timeDimension().max();

    const double hx = spaceDimensionX().step();
    const double hy = spaceDimensionY().step();
    const double ht = timeDimension().step();

    const double a1 = waveSpeed();
    const double a2 = waveSpeed();
    const double b1 = unknownB();
    const double b2 = unknownB();
    const double c  = restoration();
    const double d  = waveDissipation();

    const double w1 = weight();
    const double w2 = 1.0 - 2.0*w1;
    const double w3 = w1;

    // common parameters
    const double htht05 = ht*ht*0.50;
    const double ht_050 = ht*0.50;
    const double htht25 = ht*ht*0.25;

    // equation parameters
    const double k101 = -(0.25*a1*a1)*((ht*ht)/(hx*hx))*w1 /*+ (0.125*b1)*((ht*ht)/(hx))*w1*/;           // i-1
    const double k102 = +1.0 + (0.5*a1*a1)*((ht*ht)/(hx*hx))*w1 /*- (0.25*c)*(ht*ht)*w1*/ + 0.25*d*ht;  // i
    const double k103 = -(0.25*a1*a1)*((ht*ht)/(hx*hx))*w1 /*- (0.125*b1)*((ht*ht)/(hx))*w1*/;           // i+1

    const double k104 = +(0.25*a1*a1)*((ht*ht)/(hx*hx))*w2 /*- (0.125*b1)*((ht*ht)/(hx))*w2*/;           // i-1
    const double k105 = +2.0 - (0.5*a1*a1)*((ht*ht)/(hx*hx))*w2 /*+ (0.25*c)*(ht*ht)*w2*/;               // i
    const double k106 = +(0.25*a1*a1)*((ht*ht)/(hx*hx))*w2 /*+ (0.125*b1)*((ht*ht)/(hx))*w2*/;           // i+1

    const double k107 = +(0.25*a1*a1)*((ht*ht)/(hx*hx))*w3 /*- (0.125*b1)*((ht*ht)/(hx))*w3*/;           // i-1
    const double k108 = -1.0 - (0.5*a1*a1)*((ht*ht)/(hx*hx))*w3 /*+ (0.25*c)*(ht*ht)*w3*/ + 0.25*d*ht;  // i
    const double k109 = +(0.25*a1*a1)*((ht*ht)/(hx*hx))*w3 /*+ (0.125*b1)*((ht*ht)/(hx))*w3*/;           // i+1

    const double k110 = +(0.25*a2*a2)*((ht*ht)/(hy*hy)) /*- (0.125*b2)*((ht*ht)/(hy))*/;                 // j-1
    const double k111 = -(0.50*a2*a2)*((ht*ht)/(hy*hy));                                             // j
    const double k112 = +(0.25*a2*a2)*((ht*ht)/(hy*hy)) /*+ (0.125*b2)*((ht*ht)/(hy))*/;                 // j+1

    const double k201 = -(0.25*a2*a2)*((ht*ht)/(hy*hy))*w1 /*+ (0.125*b2)*((ht*ht)/(hy))*w1*/;           // j-1
    const double k202 = +1.0 + (0.5*a2*a2)*((ht*ht)/(hy*hy))*w1 /*- (0.25*c)*(ht*ht)*w1*/ + 0.25*d*ht;  // j
    const double k203 = -(0.25*a2*a2)*((ht*ht)/(hy*hy))*w1 /*- (0.125*b2)*((ht*ht)/(hy))*w1*/;           // j+1
    const double k204 = +(0.25*a2*a2)*((ht*ht)/(hy*hy))*w2 /*- (0.125*b2)*((ht*ht)/(hy))*w2*/;           // j-1
    const double k205 = +2.0 - (0.5*a2*a2)*((ht*ht)/(hy*hy))*w2 /*+ (0.25*c)*(ht*ht)*w2*/;               // j
    const double k206 = +(0.25*a2*a2)*((ht*ht)/(hy*hy))*w2 /*+ (0.125*b2)*((ht*ht)/(hy))*w2*/;           // j+1
    const double k207 = +(0.25*a2*a2)*((ht*ht)/(hy*hy))*w3 /*- (0.125*b2)*((ht*ht)/(hy))*w3*/;           // j-1
    const double k208 = -1.0 - (0.5*a2*a2)*((ht*ht)/(hy*hy))*w3 /*+ (0.25*c)*(ht*ht)*w3*/ + 0.25*d*ht;  // j
    const double k209 = +(0.25*a2*a2)*((ht*ht)/(hy*hy))*w3 /*+ (0.125*b2)*((ht*ht)/(hy))*w3*/;           // j+1
    const double k210 = +(0.25*a1*a1)*((ht*ht)/(hx*hx)) /*- (0.125*b2)*((ht*ht)/(hx))*/;                 // i-1
    const double k211 = -(0.50*a1*a1)*((ht*ht)/(hx*hx));                                             // i
    const double k212 = +(0.25*a1*a1)*((ht*ht)/(hx*hx)) /*+ (0.125*b2)*((ht*ht)/(hx))*/;                 // i+1

    double *ax = static_cast<double*>(malloc(sizeof(double)*(N+1)));
    double *bx = static_cast<double*>(malloc(sizeof(double)*(N+1)));
    double *cx = static_cast<double*>(malloc(sizeof(double)*(N+1)));
    double *dx = static_cast<double*>(malloc(sizeof(double)*(N+1)));
    double *rx = static_cast<double*>(malloc(sizeof(double)*(N+1)));

    for (unsigned int n=0; n<=N; n++)
    {
        ax[n] = k101;
        bx[n] = k102;
        cx[n] = k103;
    }
    ax[0] = 0.0; cx[N] = 0.0;

    double *ay = static_cast<double*>(malloc(sizeof(double)*(M+1)));
    double *by = static_cast<double*>(malloc(sizeof(double)*(M+1)));
    double *cy = static_cast<double*>(malloc(sizeof(double)*(M+1)));
    double *dy = static_cast<double*>(malloc(sizeof(double)*(M+1)));
    double *ry = static_cast<double*>(malloc(sizeof(double)*(M+1)));

    for (unsigned int m=0; m<=M; m++)
    {
        ay[m] = k201;
        by[m] = k202;
        cy[m] = k203;
    }
    ay[0] = 0.0; cy[M] = 0.0;

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

    unsigned int i=0, j=0, s=0, e=N;
    int n = 0, m = 0;
    BoundaryConditionPDE condition; double value, alpha, beta, gamma;

    for (m=ymin, sn.j = m, sn.y = m*hy, j=0; m<=ymax; ++m, sn.j=m, sn.y=m*hy, ++j)
    {
        for (n=xmin, sn.i = n, sn.x = n*hx, i=0; n<=xmax; ++n, sn.i=n, sn.x=n*hx, ++i)
        {
            u00[j][i] = initial(sn, InitialCondition::InitialValue);
        }
    }
    layerInfo(u00, tn00);

    /***********************************************************************************************/

    for (m=ymin, sn.j = m, sn.y = m*hy, j=0; m<=ymax; ++m, sn.j=m, sn.y=m*hy, ++j)
    {
        for (n=xmin, sn.i = n, sn.x = n*hx, i=0; n<=xmax; ++n, sn.i=n, sn.x=n*hx, ++i)
        {
            double firstDerivative = initial(sn, InitialCondition::InitialFirstDerivative);

            double secndDerivative = 0.0;
            if (j==0)
            {
                secndDerivative += ((a2*a2)/(hy*hy))*(+2.0*u00[0][i]-5.0*u00[1][i]+4.0*u00[2][i]-1.0*u00[3][i]);
            }
            else if (j==M)
            {
                secndDerivative += ((a2*a2)/(hy*hy))*(-1.0*u00[M-3][i]+4.0*u00[M-2][i]-5.0*u00[M-1][i]+2.0*u00[M][i]);
            }
            else
            {
                secndDerivative += ((a2*a2)/(hy*hy))*(u00[j-1][i]-2.0*u00[j][i]+u00[j+1][i]);
            }

            if (i==0)
            {
                secndDerivative += ((a1*a1)/(hx*hx))*(+2.0*u00[j][0]-5.0*u00[j][1]+4.0*u00[j][2]-1.0*u00[j][3]);
            }
            else if (i==N)
            {
                secndDerivative += ((a1*a1)/(hx*hx))*(-1.0*u00[j][N-3]+4.0*u00[j][N-2]-5.0*u00[j][N-1]+2.0*u00[j][N]);
            }
            else
            {
                secndDerivative += ((a1*a1)/(hx*hx))*(u00[j][i-1]-2.0*u00[j][i]+u00[j][i+1]);
            }

            secndDerivative += f(sn,tn00);
            secndDerivative -=  d*firstDerivative;

            u05[j][i] = u00[j][i] + firstDerivative*ht_050 + secndDerivative*0.5*0.25*ht*ht;

            u10[j][i] = u00[j][i] + firstDerivative*ht + secndDerivative*0.5*ht*ht;
        }
    }

    layerInfo(u05, tn05);

//    for (m=ymin, sn.j = m, sn.y = m*hy, j=0; m<=ymax; ++m, sn.j=m, sn.y=m*hy, ++j)
//    {
//        for (n=xmin, sn.i = n, sn.x = n*hx, i=0; n<=xmax; ++n, sn.i=n, sn.x=n*hx, ++i)
//        {
//            u10[j][i] = 2.0*u05[j][i]-u00[j][i];

//            double secndDerivative = 0.0;
//            if (j==0)
//            {
//                secndDerivative += ((a2*a2)/(hy*hy))*(+2.0*u05[0][i]-5.0*u05[1][i]+4.0*u05[2][i]-1.0*u05[3][i]);
//            }
//            else if (j==M)
//            {
//                secndDerivative += ((a2*a2)/(hy*hy))*(-1.0*u05[M-3][i]+4.0*u05[M-2][i]-5.0*u05[M-1][i]+2.0*u05[M][i]);
//            }
//            else
//            {
//                secndDerivative += ((a2*a2)/(hy*hy))*(u05[j-1][i]-2.0*u05[j][i]+u05[j+1][i]);
//            }

//            if (i==0)
//            {
//                secndDerivative += ((a1*a1)/(hx*hx))*(+2.0*u05[j][0]-5.0*u05[j][1]+4.0*u05[j][2]-1.0*u05[j][3]);
//            }
//            else if (n==xmax)
//            {
//                secndDerivative += ((a1*a1)/(hx*hx))*(-1.0*u05[j][N-3]+4.0*u05[j][N-2]-5.0*u05[j][N-1]+2.0*u05[j][N]);
//            }
//            else
//            {
//                secndDerivative += ((a1*a1)/(hx*hx))*(u05[j][i-1]-2.0*u05[j][i]+u05[j][i+1]);
//            }

//            secndDerivative += f(sn,tn05);
//            u10[j][i] += ht*ht*0.25*secndDerivative;
//            u10[j][i] += d*ht*0.25*u00[m][n];

//            u10[j][i] /= (1.0+d*ht*0.25);
//        }
//    }
    layerInfo(u10, tn10);

    /***********************************************************************************************/

    for (unsigned int ln=2; ln<=L; ln++)
    {
        TimeNodePDE tn05; tn05.i = 2*ln-3; tn05.t = 0.5*tn05.i*ht;
        TimeNodePDE tn10; tn10.i = 2*ln-2; tn10.t = 0.5*tn10.i*ht;
        TimeNodePDE tn15; tn15.i = 2*ln-1; tn15.t = 0.5*tn15.i*ht;
        TimeNodePDE tn20; tn20.i = 2*ln-0; tn20.t = 0.5*tn20.i*ht;

        /**************************************************** x direction apprx ***************************************************/

        for (m=ymin+1, sn.j=m, sn.y=m*hy, j=1; m<=ymax-1; ++m, sn.j=m, sn.y=m*hy, ++j)
        {
            for (n=xmin+1, sn.i=n, sn.x=n*hx, i=1; n<=xmax-1; ++n, sn.i=n, sn.x=n*hx, ++i)
            {
                dx[i] = 0.0;
                dx[i] += k104*u10[j][i-1] + k105*u10[j][i] + k106*u10[j][i+1];
                dx[i] += k107*u05[j][i-1] + k108*u05[j][i] + k109*u05[j][i+1];
                dx[i] += k110*u10[j-1][i] + k111*u10[j][i] + k112*u10[j+1][i];
                dx[i] += htht25 * f(sn, tn10);
            }

            sn.i = xmin; sn.x = xmin*hx;
            value = boundary(sn, tn15, condition);
            alpha = condition.alpha();
            beta  = condition.beta();
            gamma = condition.gamma();

            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                s = 1;
                u15[j][0] = (gamma/alpha)*value;
                dx[1] -= k101 * u15[j][0];
                ax[1] = ax[0] = bx[0] = cx[0] = dx[0] = rx[0] = 0.0;
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {}
            else if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {}

            sn.i = xmax; sn.x = xmax*hx;
            value = boundary(sn, tn15, condition);
            alpha = condition.alpha();
            beta  = condition.beta();
            gamma = condition.gamma();

            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                e = N-1;
                u15[j][N] = (gamma/alpha)*value;
                dx[N-1] -= k103 * u15[j][N];
                cx[N-1] = ax[N] = bx[N] = cx[N] = dx[N] = rx[N] = 0.0;
            }
            if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {}
            if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {}

            tomasAlgorithmLeft2Right(ax+s, bx+s, cx+s, dx+s, rx+s, e-s+1);
            for (unsigned int i=s; i<=e; i++) u15[j][i] = rx[i];
        }

        for (n=xmin, sn.i=n, sn.x=n*hx, i=0; n<=xmax; ++n, sn.i=n, sn.x=n*hx, ++i)
        {
            sn.j = ymin; sn.y = ymin*hy;
            value = boundary(sn, tn15, condition);
            alpha = condition.alpha();
            beta  = condition.beta();
            gamma = condition.gamma();

            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                u15[0][i] = (gamma/alpha)*value;
            }
            if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
                u15[0][i] = (3.5*u15[1][i] - 2.0*u15[2][i] + 0.5*u15[3][i] + hy*(gamma/beta)*value)/(2.0);
                //u05[0][i] = (u05[1][i] + hy*(gamma/beta)*value));
            }
            if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                u15[0][i] = (3.5*u15[1][i] - 2.0*u15[2][i] + 0.5*u15[3][i] + hy*(gamma/beta)*value)/(2.0 + (alpha/beta)*hy);
                //u05[0][i] = (u05[1][i] + hy*(gamma/beta)*value)/(1.0+hy*(alpha/beta));
            }

            sn.j = ymax; sn.y = ymax*hy;
            value = boundary(sn, tn15, condition);
            alpha = condition.alpha();
            beta  = condition.beta();
            gamma = condition.gamma();

            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                u15[M][i] = (gamma/alpha)*value;
            }
            if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
                u15[M][i] = (3.5*u15[M-1][i] - 2.0*u15[M-2][i] + 0.5*u15[M-3][i] + hy*(gamma/beta)*value)/(2.0);
                //u05[M][i] = (u05[M-1][i] + hy*(gamma/beta)*value));
            }
            if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                u15[M][i] = (3.5*u15[M-1][i] - 2.0*u15[M-2][i] + 0.5*u15[M-3][i] + hy*(gamma/beta)*value)/(2.0 + (alpha/beta)*hy);
                //u05[M][i] = (u05[M-1][i] + hy*(gamma/beta)*value)/(1.0 + hy*(alpha/beta));
            }
        }

        layerInfo(u15, tn15);

        /**************************************************** x direction apprx ***************************************************/

        /**************************************************** y direction apprx ***************************************************/

        for (n=xmin+1, sn.i=n, sn.x=n*hx, i=1; n<=xmax-1; ++n, sn.i=n, sn.x=n*hx, ++i)
        {
            for (m=ymin+1, sn.j=m, sn.y=m*hy, j=1; m<=ymax-1; ++m, sn.j=m, sn.y=m*hy, ++j)
            {
                dy[j] = 0.0;
                dy[j] += k204*u15[j-1][i] + k205*u15[j][i] + k206*u15[j+1][i];
                dy[j] += k207*u10[j-1][i] + k208*u10[j][i] + k209*u10[j+1][i];
                dy[j] += k210*u15[j][i-1] + k211*u15[j][i] + k212*u15[j][i+1];
                dy[j] += htht25 * f(sn, tn15);
            }

            sn.j = ymin; sn.y = ymin*hy;
            value = boundary(sn, tn20, condition);
            alpha = condition.alpha();
            beta  = condition.beta();
            gamma = condition.gamma();

            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                s = 1;
                u20[0][i] = (gamma/alpha)*value;
                dy[1] -= k201 * u20[0][i];
                ay[1] = ay[0] = by[0] = cy[0] = dy[0] = ry[0] = 0.0;
            }
            if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {}
            if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {}

            sn.j = ymax; sn.y = ymax*hy;
            value = boundary(sn, tn20, condition);
            alpha = condition.alpha();
            beta  = condition.beta();
            gamma = condition.gamma();

            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                e = M-1;
                u20[M][i] = (gamma/alpha)*value;
                dy[M-1] -= k203 * u20[M][i];
                cy[M-1] = ay[M] = by[M] = cy[M] = dy[M] = ry[M] = 0.0;
            }
            if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {}
            if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {}

            tomasAlgorithmLeft2Right(ay+s, by+s, cy+s, dy+s, ry+s, e-s+1);
            for (unsigned int j=s; j<=e; j++) u20[j][i] = ry[j];
        }

        for (m=ymin, sn.j=m, sn.y=m*hy, j=0; m<=ymax; ++m, sn.j=m, sn.y=m*hy, ++j)
        {
            sn.i = xmin; sn.x = xmin*hx;
            value = boundary(sn, tn20, condition);
            alpha = condition.alpha();
            beta  = condition.beta();
            gamma = condition.gamma();

            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                u20[j][0] = (gamma/alpha)*value;
            }
            if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
                u10[j][0] = (3.5*u10[j][1] - 2.0*u10[j][2] + 0.5*u10[j][3] + hx*(gamma/beta)*value)/(2.0);
                //u10[j][0] = (u10[j][1] + hx*(gamma/beta)*value);
            }
            if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                u10[j][0] = (3.5*u10[j][1] - 2.0*u10[j][2] + 0.5*u10[j][3] + hx*(gamma/beta)*value)/(2.0 + (alpha/beta)*hx);
                //u10[j][0] = (u10[j][1] + hx*(gamma/beta)*value)/(1.0+hx*(alpha/beta));
            }

            sn.i = xmax; sn.x = xmax*hx;
            value = boundary(sn, tn20, condition);
            alpha = condition.alpha();
            beta  = condition.beta();
            gamma = condition.gamma();

            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                u20[j][N] = (gamma/alpha)*value;
            }
            if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
                u10[j][N] = (3.5*u10[j][N-1] - 2.0*u10[j][N-2] + 0.5*u10[j][N-3] + hx*(gamma/beta)*value)/(2.0);
                //u10[j][N] = (u10[j][N-1] + hx*(gamma/beta)*value);
            }
            if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                u10[j][N] = (3.5*u10[j][N-1] - 2.0*u10[j][N-2] + 0.5*u10[j][N-3] + hx*(gamma/beta)*value)/(2.0 + (alpha/beta)*hx);
                //u10[j][N] = (u10[j][N-1] + hx*(gamma/beta)*value)/(1.0 + hx*(alpha/beta));
            }
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

double IWaveEquationIBVP::weight() const { return 0.25; }

//--------------------------------------------------------------------------------------------------------------//

IWaveEquationFBVP::IWaveEquationFBVP(double waveSpeed, double waveDissipation, double unknownB, double restoration) : IHyperbolicFBVP(),
    _waveSpeed(waveSpeed), _waveDissipation(waveDissipation), _unknownB(unknownB), _restoration(restoration) {}

IWaveEquationFBVP::IWaveEquationFBVP(const IWaveEquationFBVP &other) : IHyperbolicFBVP(other),
    _waveSpeed(other._waveSpeed), _waveDissipation(other._waveDissipation), _unknownB(other._unknownB), _restoration(other._restoration) {}

IWaveEquationFBVP& IWaveEquationFBVP::operator=(const IWaveEquationFBVP &other)
{
    if (this == &other) { return *this; }

    this->_waveSpeed = other._waveSpeed;
    this->_waveDissipation = other._waveDissipation;
    this->_unknownB = other._unknownB;
    this->_restoration = other._restoration;
    return *this;
}

IWaveEquationFBVP::~IWaveEquationFBVP() {}

double IWaveEquationFBVP::waveSpeed() const  { return _waveSpeed; }

double IWaveEquationFBVP::waveDissipation() const { return _waveDissipation; }

void IWaveEquationFBVP::setWaveSpeed(double waveSpeed) { _waveSpeed = waveSpeed; }

void IWaveEquationFBVP::setWaveDissipation(double waveDissipation) { _waveDissipation = waveDissipation; }

void IWaveEquationFBVP::setUnknownB(double unknownB) { this->_unknownB = unknownB; }

double IWaveEquationFBVP::unknownB() const { return _unknownB; }

void IWaveEquationFBVP::setRestoration(double restoration) { this->_restoration = restoration; }

double IWaveEquationFBVP::restoration() const { return _restoration; }

void IWaveEquationFBVP::explicit_calculate_D1V1() const {}

void IWaveEquationFBVP::implicit_calculate_D1V1() const
{
    const unsigned int N = static_cast<unsigned int>(spaceDimensionX().size()) - 1;
    const unsigned int L = static_cast<unsigned int>(timeDimension().size()) - 1;

    const int xmin = spaceDimensionX().min();
    const int xmax = spaceDimensionX().max();
    //const int tmin = timeDimension().min();
    //const int tmax = timeDimension().max();

    const double hx = spaceDimensionX().step();
    const double ht = timeDimension().step();

    const double a = waveSpeed();
    const double b = unknownB();
    const double c = restoration();
    const double d = waveDissipation();
    const double w = weight();
    const double w1 = w;
    const double w2 = 1.0 - 2.0*w;
    const double w3 = w;

    // equation parameters
    const double k11 = -(a*a)*((ht*ht)/(hx*hx))*w1 + ((b*ht*ht)/(2.0*hx))*w1;
    const double k12 = +1.0 + 2.0*(a*a)*((ht*ht)/(hx*hx))*w1 - c*ht*ht*w1 - 0.5*d*ht;
    const double k13 = -(a*a)*((ht*ht)/(hx*hx))*w1 - ((b*ht*ht)/(2.0*hx))*w1;
    const double k21 = +(a*a)*((ht*ht)/(hx*hx))*w2 - ((b*ht*ht)/(2.0*hx))*w2;
    const double k22 = +2.0 - 2.0*(a*a)*((ht*ht)/(hx*hx))*w2 + c*ht*ht*w2;
    const double k23 = +(a*a)*((ht*ht)/(hx*hx))*w2 + ((b*ht*ht)/(2.0*hx))*w2;
    const double k31 = +(a*a)*((ht*ht)/(hx*hx))*w3 - ((b*ht*ht)/(2.0*hx))*w3;
    const double k32 = -1.0 - 2.0*(a*a)*((ht*ht)/(hx*hx))*w3 + c*ht*ht*w3 - 0.5*d*ht;
    const double k33 = +(a*a)*((ht*ht)/(hx*hx))*w3 + ((b*ht*ht)/(2.0*hx))*w3;

    // left border condition parameters
    const double b11 = +1.0 + 2.0*((a*a*ht*ht)/(hx*hx))*w1 - c*ht*ht*w1 - 0.5*d*ht;
    const double b12 = -2.0*((a*a*ht*ht)/hx)*w1 + b*ht*ht*w1;
    const double b13 = -2.0*((a*a*ht*ht)/(hx*hx))*w1;
    const double b14 = +2.0 - 2.0*((a*a*ht*ht)/(hx*hx))*w2 + c*ht*ht*w2;
    const double b15 = +2.0*((a*a*ht*ht)/hx)*w2 - b*ht*ht*w2;
    const double b16 = +2.0*((a*a*ht*ht)/(hx*hx))*w2;
    const double b17 = -1.0 - 2.0*((a*a*ht*ht)/(hx*hx))*w3 + c*ht*ht*w3 - 0.5*d*ht;
    const double b18 = +2.0*((a*a*ht*ht)/hx)*w3 - b*ht*ht*w3;
    const double b19 = +2.0*((a*a*ht*ht)/(hx*hx))*w3;
    const double b20 = -2.0*((a*a*ht*ht)/hx) + b*ht*ht;

    // right border condition parameters
    const double b21 = -2.0*((a*a*ht*ht)/(hx*hx))*w1;
    const double b22 = +2.0*((a*a*ht*ht)/hx)*w1 + b*ht*ht*w1;
    const double b23 = +1.0 + 2.0*((a*a*ht*ht)/(hx*hx))*w1 - c*ht*ht*w1 - 0.5*d*ht;
    const double b24 = +2.0*((a*a*ht*ht)/(hx*hx))*w2;
    const double b25 = -2.0*((a*a*ht*ht)/hx)*w2 - b*ht*ht*w2;
    const double b26 = +2.0 - 2.0*((a*a*ht*ht)/(hx*hx))*w2 + c*ht*ht*w2;
    const double b27 = +2.0*((a*a*ht*ht)/(hx*hx))*w3;
    const double b28 = -2.0*((a*a*ht*ht)/hx)*w3 - b*ht*ht*w3;
    const double b29 = -1.0 - 2.0*((a*a*ht*ht)/(hx*hx))*w3 + c*ht*ht*w3 - 0.5*d*ht;
    const double b30 = +2.0*((a*a*ht*ht)/hx) + b*ht*ht;

    // initial condition parameters
    const double aa__hxhx = (a*a)/(hx*hx);
    const double htht_05  = 0.5*ht*ht;
    const double ht_ht    = ht*ht;

    double *ax = static_cast<double*>(malloc(sizeof(double)*(N+1)));
    double *bx = static_cast<double*>(malloc(sizeof(double)*(N+1)));
    double *cx = static_cast<double*>(malloc(sizeof(double)*(N+1)));
    double *dx = static_cast<double*>(malloc(sizeof(double)*(N+1)));
    double *rx = static_cast<double*>(malloc(sizeof(double)*(N+1)));

    for (unsigned int n=0; n<=N; n++)
    {
        ax[n] = k11;
        bx[n] = k12;
        cx[n] = k13;
    }
    ax[0] = 0.0; cx[N] = 0.0;

    DoubleVector u00(N+1);
    DoubleVector u10(N+1);
    DoubleVector u20(N+1);

    /***********************************************************************************************/
    /***********************************************************************************************/

    TimeNodePDE tn00; tn00.i = L;   tn00.t = tn00.i*ht;
    TimeNodePDE tn10; tn10.i = L-1; tn10.t = tn10.i*ht;

    SpaceNodePDE sn;
    unsigned int i = 0;
    for (int n=xmin; n<=xmax; n++, i++)
    {
        sn.i = n; sn.x = n*hx;
        u00[i] = final(sn, FinalCondition::FinalValue);
    }
    layerInfo(u00, tn00);

    /***********************************************************************************************/

    i = 0;
    for (int n=xmin; n<=xmax; n++, i++)
    {
        sn.i = n; sn.x = n*hx;
        double firstDerivativeT = final(sn, FinalCondition::FinalFirstDerivative);
        u10[i] = u00[i] - firstDerivativeT*ht;
        double secndDerivativeX = 0.0;
        double firstDerivativeX = 0.0;
        if (n==xmin)
        {
            secndDerivativeX = aa__hxhx*(+2.0*u00[0]-5.0*u00[1]+4.0*u00[2]-1.0*u00[3]);
            firstDerivativeX = (-3.0*u00[0]+4.0*u00[1]-u00[2])/(2.0*hx);
        }
        else if (n==xmax)
        {
            secndDerivativeX = aa__hxhx*(-1.0*u00[N-3]+4.0*u00[N-2] -5.0*u00[N-1] +2.0*u00[N]);
            firstDerivativeX = (u00[N-2]-4.0*u00[N-1]+3.0*u00[N])/(2.0*hx);
        }
        else
        {
            secndDerivativeX = aa__hxhx*(u00[i-1]-2.0*u00[i]+u00[i+1]);
            firstDerivativeX = (u00[i+1]-u00[i-1])/(2.0*hx);
        }
        u10[i] += (secndDerivativeX + b*firstDerivativeX + c*u00[i] + f(sn,tn00) - d*firstDerivativeT) * htht_05;
    }
    layerInfo(u10, tn10);

    /***********************************************************************************************/
    /***********************************************************************************************/

    for (unsigned int ln=L-2, size_ln = static_cast<unsigned int>(0)-1; ln != size_ln; ln--)
    {
        TimeNodePDE tn00; tn00.i = ln+2; tn00.t = tn00.i*ht;
        TimeNodePDE tn10; tn10.i = ln+1; tn10.t = tn10.i*ht;
        TimeNodePDE tn20; tn20.i = ln;   tn20.t = tn20.i*ht;

        unsigned int i = 1;
        for (int n=xmin+1; n<=xmax-1; n++, i++)
        {
            sn.i = n; sn.x = n*hx;
            dx[i] = 0.0;
            dx[i] += k21 * u10[i-1] + k22 * u10[i] + k23 * u10[i+1];
            dx[i] += k31 * u00[i-1] + k32 * u00[i] + k33 * u00[i+1];

            dx[i] += ht_ht*f(sn, tn10);
            //dx[i] += ht_ht*w1*f(sn, tn20) + ht_ht*w2*f(sn, tn10) + ht_ht*w3*f(sn, tn00);
        }

        unsigned int s=0, e=N;
        BoundaryConditionPDE condition; double alpha, beta, gamma, value;

        sn.i = xmin; sn.x = xmin*hx;
        value = boundary(sn, tn20, condition);
        alpha = condition.alpha();
        beta  = condition.beta();
        gamma = condition.gamma();
        if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
        {
            s = 1;

            u20[0] = (gamma/alpha)*value;
            dx[1] -= k11 * u20[0];
            ax[1] = ax[0] = bx[0] = cx[0] = dx[0] = rx[0] = 0.0;
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
        {
            s = 0;

            ax[s]  = 0.0;
            bx[s]  = beta * b11 + alpha * b12;
            cx[s]  = beta * b13;

            dx[s]  = beta * b14 * u10[s] + beta * b16 * u10[s+1];
            dx[s] += beta * b17 * u00[s] + beta * b19 * u00[s+1];

            dx[s] += gamma * b20 * boundary(sn, tn10, condition);
            dx[s] += beta  * ht_ht*f(sn,tn10);

            //dx[s] += gamma * b20 * (w1*value + w2*boundary(sn, tn10, condition) + w3*boundary(sn, tn00, condition));
            //dx[s] += beta * (ht_ht*w1*f(sn, tn20) + ht_ht*w2*f(sn, tn10) + ht_ht*w3*f(sn, tn00));
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Robin)
        {
            s = 0;

            //ax[s] = 0.0;
            //bx[s] = hx*alpha-beta;
            //cx[s] = beta;
            //dx[s] = hx*gamma*value;

            ax[s]  = 0.0;
            bx[s]  = beta * b11 + alpha * b12;
            cx[s]  = beta * b13;

            dx[s]  = beta * b14 * u10[s] + alpha * b15 * u10[s] + beta * b16 * u10[s+1];
            dx[s] += beta * b17 * u00[s] + alpha * b18 * u00[s] + beta * b19 * u00[s+1];

            dx[s] += gamma * b20 * boundary(sn, tn10, condition);
            dx[s] += beta  * ht_ht*f(sn,tn10);

            //dx[s] += gamma * b20 * (w1*value + w2*boundary(sn, tn10, condition) + w3*boundary(sn, tn00, condition));
            //dx[s] += beta * (ht_ht*w1*f(sn, tn20) + ht_ht*w2*f(sn, tn10) + ht_ht*w3*f(sn, tn00));
        }

        sn.i = xmax; sn.x = xmax*hx;
        value = boundary(sn, tn20, condition);
        alpha = condition.alpha();
        beta  = condition.beta();
        gamma = condition.gamma();
        if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
        {
            e = N-1;

            u20[N] = (gamma/alpha)*value;
            dx[N-1] -= k13 * u20[N];
            cx[N-1] = ax[N] = bx[N] = cx[N] = dx[N] = rx[N] = 0.0;
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
        {
            e = N;

            ax[e]  = beta * b21;
            bx[e]  = beta * b23 + alpha * b22;
            cx[e]  = 0.0;

            dx[e]  = beta * b24 * u10[e-1] + beta * b26 * u10[e];
            dx[e] += beta * b27 * u00[e-1] + beta * b29 * u00[e];

            dx[e] += gamma * b30 * boundary(sn, tn10, condition);
            dx[e] += beta  * ht_ht*f(sn,tn10);

            //dx[e] += gamma * b30 * (w1*value + w2*boundary(sn, tn10, condition) + w3*boundary(sn, tn00, condition));
            //dx[e] += beta * (ht_ht*w1*f(sn, tn20) + ht_ht*w2*f(sn, tn10) + ht_ht*w3*f(sn, tn00));
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Robin)
        {
            e = N;

            //ax[e] = -beta;
            //bx[e] = hx*alpha+beta;
            //cx[e] = 0.0;
            //dx[e] = hx*gamma*value;

            ax[e]  = beta * b21;
            bx[e]  = beta * b23 + alpha * b22;
            cx[e]  = 0.0;

            dx[e]  = beta * b24 * u10[e-1] + alpha * b25 * u10[e] + beta * b26 * u10[e];
            dx[e] += beta * b27 * u00[e-1] + alpha * b28 * u00[e] + beta * b29 * u00[e];

            dx[e] += gamma * b30 * boundary(sn, tn10, condition);
            dx[e] += beta  * ht_ht*f(sn,tn10);

            //dx[e] += gamma * b30 * (w1*value + w2*boundary(sn, tn10, condition) + w3*boundary(sn, tn00, condition));
            //dx[e] += beta * (ht_ht*w1*f(sn, tn20) + ht_ht*w2*f(sn, tn10) + ht_ht*w3*f(sn, tn00));
        }

        tomasAlgorithmLeft2Right(ax+s, bx+s, cx+s, dx+s, rx+s, e-s+1);
        for (unsigned int n=s; n<=e; n++) u20[n] = rx[n];
        layerInfo(u20, tn20);

        for (unsigned int i=0; i<=N; i++)
        {
            u00[i] = u10[i];
            u10[i] = u20[i];
        }
    }

    u00.clear();
    u10.clear();
    u20.clear();

    free(rx);
    free(dx);
    free(cx);
    free(bx);
    free(ax);
}

void IWaveEquationFBVP::explicit_calculate_D2V1() const
{
    const Dimension _spaceDimensionX = spaceDimensionX();
    const Dimension _spaceDimensionY = spaceDimensionX();
    const Dimension _timeDimension = timeDimension();

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

void IWaveEquationFBVP::implicit_calculate_D2V1() const
{
    const Dimension _spaceDimensionX = spaceDimensionX();
    const Dimension _spaceDimensionY = spaceDimensionX();
    const Dimension _timeDimension = timeDimension();

    const unsigned int N = static_cast<unsigned int>( _spaceDimensionX.size() - 1 );
    const unsigned int M = static_cast<unsigned int>( _spaceDimensionY.size() - 1 );
    const unsigned int L = static_cast<unsigned int>( _timeDimension.size() - 1 );

    const double hx = _spaceDimensionX.step();
    const double hy = _spaceDimensionY.step();
    const double ht = _timeDimension.step();

    const double _lambda = weight();
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

double IWaveEquationFBVP::weight() const { return 0.25; }

//--------------------------------------------------------------------------------------------------------------//
