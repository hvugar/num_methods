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

        BoundaryConditionPDE condition; double alpha, beta, /*gamma,*/ value;

        for (unsigned int n=1; n<=N-1; n++)
        {
            sn.i = static_cast<int>(n); sn.x = n*hx;

            u20[n]  = 2.0*u10[n] - u00[n] + wvdsp_ht_05*u00[n] + ht_ht*f(sn, tn10);
            u20[n] += (u10[n-1] - 2.0*u10[n] + u10[n+1])*p_aa_htht__hxhx;
            u20[n] /= (1.0+wvdsp_ht_05);
        }

        sn.i = static_cast<int>(0); sn.x = 0*hx;
        value = boundary(sn, tn10, condition);
        if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
        {
            u20[0] = /*(gamma/alpha)**/value;
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
        {
            u20[0]  = /*beta **/(2.0*u10[0] - u00[0] + wvdsp_ht_05*u00[0]);
            u20[0] += /*beta **/(2.0*p_aa_htht__hxhx*(u10[1]-u10[0]));
            u20[0] += /*beta **/(ht_ht*f(sn,tn10));
            u20[0] += /*gamma**/(-2.0*p_aa_htht__hx)*value;
            u20[0] /= /*beta **/(1.0+wvdsp_ht_05);
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Robin)
        {
            alpha = condition.alpha();
            beta  = condition.beta();

            u20[0]  = beta *(2.0*u10[0] - u00[0] + wvdsp_ht_05*u00[0]);
            u20[0] += beta *(2.0*p_aa_htht__hxhx*(u10[1]-u10[0]));
            u20[0] += beta *(ht_ht*f(sn,tn10));
            u20[0] += /*gamma**/(-2.0*p_aa_htht__hx)*value;
            u20[0] += alpha*(+2.0*p_aa_htht__hx)*u10[0];
            u20[0] /= beta*(1.0+wvdsp_ht_05);
        }

        sn.i = static_cast<int>(N); sn.x = N*hx;
        value = boundary(sn, tn10, condition);
        if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
        {
            u20[N] = /*(gamma/alpha)**/value;
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
        {
            u20[N]  = /*beta **/(2.0*u10[N] - u00[N] + wvdsp_ht_05*u00[N]);
            u20[N] += /*beta **/(2.0*p_aa_htht__hxhx*(u10[N-1]-u10[N]));
            u20[N] += /*beta **/(ht_ht*f(sn,tn10));
            u20[N] += /*gamma**/(+2.0*p_aa_htht__hx)*value;
            u20[N] /= /*beta **/(1.0+wvdsp_ht_05);
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Robin)
        {
            alpha = condition.alpha();
            beta  = condition.beta();
            /*gamma = condition.gamma();*/

            u20[N]  = beta *(2.0*u10[N] - u00[N] + wvdsp_ht_05*u00[N]);
            u20[N] += beta *(2.0*p_aa_htht__hxhx*(u10[N-1]-u10[N]));
            u20[N] += beta *(ht_ht*f(sn,tn10));
            u20[N] += /*gamma**/(+2.0*p_aa_htht__hx)*value;
            u20[N] += alpha*(-2.0*p_aa_htht__hx)*u10[N];
            u20[N] /= beta *(1.0+wvdsp_ht_05);
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
    const double b__2_0hx = b/(2.0*hx);
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
        double fdt = initial(sn, InitialCondition::InitialFirstDerivative);

        double sdx = 0.0;
        double fdx = 0.0;
        if (n==xmin)
        {
            sdx = aa__hxhx*(+2.0*u00[0]-5.0*u00[1]+4.0*u00[2]-1.0*u00[3]);
            fdx = b__2_0hx*(-3.0*u00[0]+4.0*u00[1]-u00[2]);
        }
        else if (n==xmax)
        {
            sdx = aa__hxhx*(-1.0*u00[N-3]+4.0*u00[N-2] -5.0*u00[N-1] +2.0*u00[N]);
            fdx = b__2_0hx*(u00[N-2]-4.0*u00[N-1]+3.0*u00[N]);
        }
        else
        {
            sdx = aa__hxhx*(u00[i-1]-2.0*u00[i]+u00[i+1]);
            fdx = b__2_0hx*(u00[i+1]-u00[i-1]);
        }

        double sdt = sdx + fdx + c*u00[i] + f(sn, tn00) - d*fdt;

        u10[i] = u00[i] + fdt*ht + sdt*htht_05;
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
        BoundaryConditionPDE condition; double alpha, beta, /*gamma,*/ value;

        sn.i = xmin; sn.x = xmin*hx;
        value = boundary(sn, tn20, condition);

        if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
        {
            s = 1;

            u20[0] = /*(gamma/alpha)**/value;
            dx[1] -= k11 * u20[0];
            ax[1] = ax[0] = bx[0] = cx[0] = dx[0] = rx[0] = 0.0;
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
        {
            s = 0;

            ax[s]  = 0.0;
            bx[s]  = /*beta **/ b11 + alpha * b12;
            cx[s]  = /*beta **/ b13;

            dx[s]  = /*beta **/ b14 * u10[s] + beta * b16 * u10[s+1];
            dx[s] += /*beta **/ b17 * u00[s] + beta * b19 * u00[s+1];

            dx[s] += /*gamma **/ b20 * boundary(sn, tn10, condition);
            dx[s] += /*beta  **/ ht_ht*f(sn,tn10);

            //dx[s] += gamma * b20 * (w1*value + w2*boundary(sn, tn10, condition) + w3*boundary(sn, tn00, condition));
            //dx[s] += beta * (ht_ht*w1*f(sn, tn20) + ht_ht*w2*f(sn, tn10) + ht_ht*w3*f(sn, tn00));
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Robin)
        {
            s = 0;

            alpha = condition.alpha();
            beta  = condition.beta();
            /*gamma = condition.gamma();*/

            //ax[s] = 0.0;
            //bx[s] = hx*alpha-beta;
            //cx[s] = beta;
            //dx[s] = hx*gamma*value;

            ax[s]  = 0.0;
            bx[s]  = beta * b11 + alpha * b12;
            cx[s]  = beta * b13;

            dx[s]  = beta * b14 * u10[s] + alpha * b15 * u10[s] + beta * b16 * u10[s+1];
            dx[s] += beta * b17 * u00[s] + alpha * b18 * u00[s] + beta * b19 * u00[s+1];

            dx[s] += /*gamma * */b20 * boundary(sn, tn10, condition);
            dx[s] += beta  * ht_ht*f(sn,tn10);

            //dx[s] += gamma * b20 * (w1*value + w2*boundary(sn, tn10, condition) + w3*boundary(sn, tn00, condition));
            //dx[s] += beta * (ht_ht*w1*f(sn, tn20) + ht_ht*w2*f(sn, tn10) + ht_ht*w3*f(sn, tn00));
        }

        sn.i = xmax; sn.x = xmax*hx;
        value = boundary(sn, tn20, condition);

        if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
        {
            e = N-1;

            u20[N] = /*(gamma/alpha)**/value;
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

            dx[e] += /*gamma * */b30 * boundary(sn, tn10, condition);
            dx[e] += beta  * ht_ht*f(sn,tn10);

            //dx[e] += gamma * b30 * (w1*value + w2*boundary(sn, tn10, condition) + w3*boundary(sn, tn00, condition));
            //dx[e] += beta * (ht_ht*w1*f(sn, tn20) + ht_ht*w2*f(sn, tn10) + ht_ht*w3*f(sn, tn00));
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Robin)
        {
            e = N;

            alpha = condition.alpha();
            beta  = condition.beta();
            /*gamma = condition.gamma();*/

            //ax[e] = -beta;
            //bx[e] = hx*alpha+beta;
            //cx[e] = 0.0;
            //dx[e] = hx*gamma*value;

            ax[e]  = beta * b21;
            bx[e]  = beta * b23 + alpha * b22;
            cx[e]  = 0.0;

            dx[e]  = beta * b24 * u10[e-1] + alpha * b25 * u10[e] + beta * b26 * u10[e];
            dx[e] += beta * b27 * u00[e-1] + alpha * b28 * u00[e] + beta * b29 * u00[e];

            dx[e] += /*gamma * */b30 * boundary(sn, tn10, condition);
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
    const unsigned int N = static_cast<unsigned int>(spaceDimensionX().size()) - 1;
    const unsigned int M = static_cast<unsigned int>(spaceDimensionY().size()) - 1;
    const unsigned int L = static_cast<unsigned int>(timeDimension().size()) - 1;

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

    double ht_max = 1.0/sqrt((a1*a1)/(hx*hx)+(a2*a2)/(hy*hy));
    if (ht > ht_max) { throw std::runtime_error("Differential scheme is conditionally steady."); }

    // common parameters
    const double htht10 = ht*ht;
    const double htht05 = ht*ht*0.5;

    double k1 = +1.0 + 0.5*d*ht;
    double k2 = -1.0 + 0.5*d*ht;
    double k3 = (a1*a1)*((ht*ht)/(hx*hx)) - b1*((ht*ht)/(2.0*hx));
    double k4 = (a1*a1)*((ht*ht)/(hx*hx)) + b1*((ht*ht)/(2.0*hx));
    double k5 = (a2*a2)*((ht*ht)/(hy*hy)) - b2*((ht*ht)/(2.0*hy));
    double k6 = (a2*a2)*((ht*ht)/(hy*hy)) + b2*((ht*ht)/(2.0*hy));
    double k7 = (-2.0*a1*a1)*((ht*ht)/(hx*hx)) + (-2.0*a2*a2)*((ht*ht)/(hy*hy)) + c*ht*ht + 2.0;

    // initial condition parameters
    const double a1a1__hxhx = (a1*a1)/(hx*hx);
    const double a2a2__hyhy = (a2*a2)/(hy*hy);
    const double b1__20_0hx = (b1)/(2.0*hx);
    const double b2__20_0hy = (b2)/(2.0*hy);

    double **u00 = static_cast<double**>(malloc(sizeof(double*)*(M+1)));
    double **u10 = static_cast<double**>(malloc(sizeof(double*)*(M+1)));
    double **u20 = static_cast<double**>(malloc(sizeof(double*)*(M+1)));
    for (unsigned int i=0; i<=M; i++)
    {
        u00[i] = static_cast<double*>(malloc(sizeof(double)*(N+1)));
        u10[i] = static_cast<double*>(malloc(sizeof(double)*(N+1)));
        u20[i] = static_cast<double*>(malloc(sizeof(double)*(N+1)));
    }

    //DoubleMatrix u00(M+1, N+1);
    //DoubleMatrix u10(M+1, N+1);
    //DoubleMatrix u20(M+1, N+1);

    auto calculate_border = [&xmin, &ymin, &xmax, &ymax, &hx, &hy, &N, &M](SpaceNodePDE sn, TimeNodePDE tn, double **u, unsigned int i, unsigned int j, const IWaveEquationIBVP* w)
    {
        BoundaryConditionPDE condition; double value, alpha, beta, gamma;

        value = w->boundary(sn, tn, condition);

        if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
        {
            u[j][i] = (gamma/alpha)*value;
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
        {
            if (sn.i==xmin && sn.j==ymin)
            {

            }
            else if (sn.i==xmin && sn.j==ymax)
            {

            }
            else if (sn.i==xmax && sn.j==ymax)
            {

            }
            else if (sn.i==xmax && sn.j==ymin)
            {

            }
            else if (sn.i==xmin)
            {
                u[j][0] = (3.5*u[j][1] - 2.0*u[j][2] + 0.5*u[j][3] + hx*(gamma/beta)*value)/(2.0);
            }
            else if (sn.i==xmax)
            {
                u[j][N] = (3.5*u[j][N-1] - 2.0*u[j][N-2] + 0.5*u[j][N-3] + hx*(gamma/beta)*value)/(2.0);
            }
            else if (sn.j==ymin)
            {
                u[0][i] = (3.5*u[1][i] - 2.0*u[2][i] + 0.5*u[3][i] + hy*(gamma/beta)*value)/(2.0);
            }
            else if (sn.j==ymax)
            {
                u[M][i] = (3.5*u[M-1][i] - 2.0*u[M-2][i] + 0.5*u[M-3][i] + hy*(gamma/beta)*value)/(2.0);
            }
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Robin)
        {
            alpha = condition.alpha();
            beta  = condition.beta();
            /*gamma = condition.gamma();*/

            if (sn.i==xmin && sn.j==ymin)
            {

            }
            else if (sn.i==xmin && sn.j==ymax)
            {

            }
            else if (sn.i==xmax && sn.j==ymax)
            {

            }
            else if (sn.i==xmax && sn.j==ymin)
            {

            }
            else if (sn.i==xmin)
            {
                u[j][0] = (3.5*u[j][1] - 2.0*u[j][2] + 0.5*u[j][3] + hx*(gamma/beta)*value)/(2.0 + (alpha/beta)*hx);
            }
            else if (sn.i==xmax)
            {
                u[j][N] = (3.5*u[j][N-1] - 2.0*u[j][N-2] + 0.5*u[j][N-3] + hx*(gamma/beta)*value)/(2.0 + (alpha/beta)*hx);
            }
            else if (sn.j==ymin)
            {
                u[0][i] = (3.5*u[1][i] - 2.0*u[2][i] + 0.5*u[3][i] + hy*(gamma/beta)*value)/(2.0 + (alpha/beta)*hy);
            }
            else if (sn.j==ymax)
            {
                u[M][i] = (3.5*u[M-1][i] - 2.0*u[M-2][i] + 0.5*u[M-3][i] + hy*(gamma/beta)*value)/(2.0 + (alpha/beta)*hy);
            }
        }
    };

    /***********************************************************************************************/
    /***********************************************************************************************/

    TimeNodePDE tn00; tn00.i = 0; tn00.t = 0.5*tn00.i*ht;
    TimeNodePDE tn10; tn10.i = 2; tn10.t = 0.5*tn10.i*ht;
    SpaceNodePDE sn;

    unsigned int i=0, j=0;
    int n = 0, m = 0;
    //BoundaryConditionPDE condition; double value, alpha, beta, gamma;

    for (m=ymin, sn.j = m, sn.y = m*hy, j=0; m<=ymax; ++m, sn.j=m, sn.y=m*hy, ++j)
    {
        for (n=xmin, sn.i = n, sn.x = n*hx, i=0; n<=xmax; ++n, sn.i=n, sn.x=n*hx, ++i)
        {
            u00[j][i] = initial(sn, InitialCondition::InitialValue);
        }
    }

    layerInfo(DoubleMatrix(u00, M+1, N+1), tn00);

    /***********************************************************************************************/

    for (m=ymin, sn.j = m, sn.y = m*hy, j=0; m<=ymax; ++m, sn.j=m, sn.y=m*hy, ++j)
    {
        for (n=xmin, sn.i = n, sn.x = n*hx, i=0; n<=xmax; ++n, sn.i=n, sn.x=n*hx, ++i)
        {
            double fdt = initial(sn, InitialCondition::InitialFirstDerivative);

            double sdxy = 0.0;
            double fdxy = 0.0;
            if (j==0)
            {
                sdxy += a2a2__hyhy*(+2.0*u00[0][i]-5.0*u00[1][i]+4.0*u00[2][i]-1.0*u00[3][i]);
                fdxy += b2__20_0hy*(-3.0*u00[0][i]+4.0*u00[1][i]-1.0*u00[2][i]);
            }
            else if (j==M)
            {
                sdxy += a2a2__hyhy*(-1.0*u00[M-3][i]+4.0*u00[M-2][i]-5.0*u00[M-1][i]+2.0*u00[M][i]);
                fdxy += b2__20_0hy*(u00[M-2][i]-4.0*u00[M-1][i]+3.0*u00[M][i]);
            }
            else
            {
                sdxy += a2a2__hyhy*(u00[j-1][i]-2.0*u00[j][i]+u00[j+1][i]);
                fdxy += b2__20_0hy*(u00[j+1][i]-u00[j-1][i]);
            }

            if (i==0)
            {
                sdxy += a1a1__hxhx*(+2.0*u00[j][0]-5.0*u00[j][1]+4.0*u00[j][2]-1.0*u00[j][3]);
                fdxy += b1__20_0hx*(-3.0*u00[j][0]+4.0*u00[j][1]-1.0*u00[j][2]);

            }
            else if (i==N)
            {
                sdxy += a1a1__hxhx*(-1.0*u00[j][N-3]+4.0*u00[j][N-2]-5.0*u00[j][N-1]+2.0*u00[j][N]);
                fdxy += b1__20_0hx*(u00[j][N-2]-4.0*u00[j][N-1]+3.0*u00[j][N]);
            }
            else
            {
                sdxy += a1a1__hxhx*(u00[j][i-1]-2.0*u00[j][i]+u00[j][i+1]);
                fdxy += b1__20_0hx*(u00[j][i+1]-u00[j][i-1]);
            }

            double sdt = sdxy + fdxy + c*u00[j][i] + f(sn,tn00) - d*fdt;
            u10[j][i] = u00[j][i] + fdt*ht + sdt*htht05;
        }
    }

    layerInfo(DoubleMatrix(u10, M+1, N+1), tn10);

    /***********************************************************************************************/

    for (unsigned int ln=2; ln<=L; ln++)
    {
        //TimeNodePDE tn00; tn00.i = ln-2; tn00.t = tn00.i*ht;
        TimeNodePDE tn10; tn10.i = 2*(ln-1); tn10.t = 0.5*tn10.i*ht;
        TimeNodePDE tn20; tn20.i = 2*(ln-0); tn20.t = 0.5*tn20.i*ht;

        for (m=ymin+1, sn.j=m, sn.y=m*hy, j=1; m<=ymax-1; ++m, sn.j=m, sn.y=m*hy, ++j)
        {
            for (n=xmin+1, sn.i=n, sn.x=n*hx, i=1; n<=xmax-1; ++n, sn.i=n, sn.x=n*hx, ++i)
            {
                u20[j][i] = k2*u00[j][i] + k3*u10[j][i-1] + k4*u10[j][i+1] + k5*u10[j-1][i] + k6*u10[j+1][i] + k7*u10[j][i] + htht10*f(sn, tn10);
                u20[j][i] /= k1;
            }
        }

        /**************************************************** border conditions ***************************************************/

        SpaceNodePDE sn1, sn2;
        sn1.j = ymin; sn1.y = ymin*hy;
        sn2.j = ymax; sn2.y = ymax*hy;
        for (n=xmin+1, sn1.i=sn2.i=n, sn1.x=sn2.x=n*hx, i=1; n<=xmax-1; ++n, sn1.i=sn2.i=n, sn1.x=sn2.x=n*hx, ++i)
        {
            calculate_border(sn1, tn20, u20, i, 0, this);
            calculate_border(sn2, tn20, u20, i, M, this);
        }

        sn1.i = xmin; sn1.x = xmin*hx;
        sn2.i = xmax; sn2.x = xmax*hx;
        for (m=ymin+1, sn1.j=sn2.j=m, sn1.y=sn2.y=m*hy, j=1; m<=ymax-1; ++m, sn1.j=sn2.j=m, sn1.y=sn2.y=m*hy, ++j)
        {
            calculate_border(sn1, tn20, u20, 0, j, this);
            calculate_border(sn2, tn20, u20, N, j, this);
        }

        sn.i = xmin; sn.x = xmin*hx; sn.j = ymin; sn.y = ymin*hy;
        calculate_border(sn, tn20, u20, 0, 0, this);
        sn.i = xmin; sn.x = xmin*hx; sn.j = ymax; sn.y = ymax*hy;
        calculate_border(sn, tn20, u20, 0, M, this);
        sn.i = xmax; sn.x = xmax*hx; sn.j = ymax; sn.y = ymax*hy;
        calculate_border(sn, tn20, u20, N, M, this);
        sn.i = xmax; sn.x = xmax*hx; sn.j = ymin; sn.y = ymin*hy;
        calculate_border(sn, tn20, u20, N, 0, this);

        /**************************************************** border conditions ***************************************************/

        layerInfo(DoubleMatrix(u20, M+1, N+1), tn20);

        double **_tmp = u00; u00 = u10; u10 = u20; u20 = _tmp;

        //        for (unsigned int m=0; m<=M; m++)
        //        {
        //            for (unsigned int n=0; n<=N; n++)
        //            {
        //                u00[m][n] = u10[m][n];
        //                u10[m][n] = u20[m][n];
        //            }
        //        }
    }


    for (unsigned int i=0; i<=M; i++) { free(u00[i]); free(u10[i]); free(u20[i]); }
    free(u00); free(u10); free(u20);

    //u00.clear();
    //u10.clear();
    //u20.clear();
}

void IWaveEquationIBVP::implicit_calculate_D2V1() const
{
    const unsigned int N = static_cast<unsigned int>(spaceDimensionX().size()) - 1;
    const unsigned int M = static_cast<unsigned int>(spaceDimensionY().size()) - 1;
    const unsigned int L = static_cast<unsigned int>(timeDimension().size()) - 1;

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
    const double htht_0500 = ht*ht*0.500;
    const double htht_0250 = ht*ht*0.250;
    const double htht_0125 = ht*ht*0.125;
    const double ht_0500 = ht*0.500;
    //const double ht_0250 = ht*0.250;

    // equation parameters
    const double k101 = -(0.25*a1*a1)*((ht*ht)/(hx*hx))*w1 + (0.125*b1)*((ht*ht)/(hx))*w1;           // i-1
    const double k102 = +1.0 + (0.5*a1*a1)*((ht*ht)/(hx*hx))*w1 - (0.25*c)*(ht*ht)*w1 + 0.25*d*ht;   // i
    const double k103 = -(0.25*a1*a1)*((ht*ht)/(hx*hx))*w1 - (0.125*b1)*((ht*ht)/(hx))*w1;           // i+1
    const double k104 = +(0.25*a1*a1)*((ht*ht)/(hx*hx))*w2 - (0.125*b1)*((ht*ht)/(hx))*w2;           // i-1
    const double k105 = +2.0 - (0.5*a1*a1)*((ht*ht)/(hx*hx))*w2 + (0.25*c)*(ht*ht)*w2;               // i
    const double k106 = +(0.25*a1*a1)*((ht*ht)/(hx*hx))*w2 + (0.125*b1)*((ht*ht)/(hx))*w2;           // i+1
    const double k107 = +(0.25*a1*a1)*((ht*ht)/(hx*hx))*w3 - (0.125*b1)*((ht*ht)/(hx))*w3;           // i-1
    const double k108 = -1.0 - (0.5*a1*a1)*((ht*ht)/(hx*hx))*w3 + (0.25*c)*(ht*ht)*w3 + 0.25*d*ht;   // i
    const double k109 = +(0.25*a1*a1)*((ht*ht)/(hx*hx))*w3 + (0.125*b1)*((ht*ht)/(hx))*w3;           // i+1
    const double k110 = +(0.25*a2*a2)*((ht*ht)/(hy*hy)) - (0.125*b2)*((ht*ht)/(hy));                 // j-1
    const double k111 = -(0.50*a2*a2)*((ht*ht)/(hy*hy));                                             // j
    const double k112 = +(0.25*a2*a2)*((ht*ht)/(hy*hy)) + (0.125*b2)*((ht*ht)/(hy));                 // j+1

    const double k201 = -(0.25*a2*a2)*((ht*ht)/(hy*hy))*w1 + (0.125*b2)*((ht*ht)/(hy))*w1;           // j-1
    const double k202 = +1.0 + (0.5*a2*a2)*((ht*ht)/(hy*hy))*w1 - (0.25*c)*(ht*ht)*w1 + 0.25*d*ht;   // j
    const double k203 = -(0.25*a2*a2)*((ht*ht)/(hy*hy))*w1 - (0.125*b2)*((ht*ht)/(hy))*w1;           // j+1
    const double k204 = +(0.25*a2*a2)*((ht*ht)/(hy*hy))*w2 - (0.125*b2)*((ht*ht)/(hy))*w2;           // j-1
    const double k205 = +2.0 - (0.5*a2*a2)*((ht*ht)/(hy*hy))*w2 + (0.25*c)*(ht*ht)*w2;               // j
    const double k206 = +(0.25*a2*a2)*((ht*ht)/(hy*hy))*w2 + (0.125*b2)*((ht*ht)/(hy))*w2;           // j+1
    const double k207 = +(0.25*a2*a2)*((ht*ht)/(hy*hy))*w3 - (0.125*b2)*((ht*ht)/(hy))*w3;           // j-1
    const double k208 = -1.0 - (0.5*a2*a2)*((ht*ht)/(hy*hy))*w3 + (0.25*c)*(ht*ht)*w3 + 0.25*d*ht;   // j
    const double k209 = +(0.25*a2*a2)*((ht*ht)/(hy*hy))*w3 + (0.125*b2)*((ht*ht)/(hy))*w3;           // j+1
    const double k210 = +(0.25*a1*a1)*((ht*ht)/(hx*hx)) - (0.125*b1)*((ht*ht)/(hx));                 // i-1
    const double k211 = -(0.50*a1*a1)*((ht*ht)/(hx*hx));                                             // i
    const double k212 = +(0.25*a1*a1)*((ht*ht)/(hx*hx)) + (0.125*b1)*((ht*ht)/(hx));                 // i+1

    // left border condition parameters x = xmin
    const double b111 = +1.0 + (0.25*a1*a1*((ht*ht)/(0.5*hx*hx)))*w1 - ((0.25*c)*(ht*ht))*w1 + 0.25*d*ht;
    const double b112 = +(0.25*a1*a1*(ht*ht)/(0.5*hx))*w1 - (0.25*b1*ht*ht)*w1;
    const double b113 = -((0.25*a1*a1*ht*ht)/(0.5*hx*hx))*w1;
    const double b114 = +2.0 - (0.25*a1*a1*((ht*ht)/(0.5*hx*hx)))*w2 + ((0.25*c)*(ht*ht))*w2;
    const double b115 = -(0.25*a1*a1*(ht*ht)/(0.5*hx))*w2 + (0.25*b1*ht*ht)*w2;
    const double b116 = +((0.25*a1*a1*ht*ht)/(0.5*hx*hx))*w2;
    const double b117 = -1.0 - (0.25*a1*a1*((ht*ht)/(0.5*hx*hx)))*w3 + ((0.25*c)*(ht*ht))*w3 + 0.25*d*ht;
    const double b118 = -(0.25*a1*a1*(ht*ht)/(0.5*hx))*w3 + (0.25*b1*ht*ht)*w3;
    const double b119 = +((0.25*a1*a1*ht*ht)/(0.5*hx*hx))*w3;
    const double b120 = +((0.25*a2*a2*ht*ht)/(hy*hy));
    const double b121 = +((0.25*b2*ht*ht)/(2.0*hy));
    const double b122 = +(0.25*a1*a1*ht*ht)/(0.5*hx) - (0.25*b1*ht*ht);

    // right border condition parameters x = xmax
    const double b211 = -((0.25*a1*a1*ht*ht)/(0.5*hx*hx))*w1;
    const double b212 = +(0.25*a1*a1*(ht*ht)/(0.5*hx))*w1 + (0.25*b1*ht*ht)*w1;
    const double b213 = +1.0 + (0.25*a1*a1*((ht*ht)/(0.5*hx*hx)))*w1 - ((0.25*c)*(ht*ht))*w1 + 0.25*d*ht;
    const double b214 = +((0.25*a1*a1*ht*ht)/(0.5*hx*hx))*w2;
    const double b215 = -(0.25*a1*a1*(ht*ht)/(0.5*hx))*w2 - (0.25*b1*ht*ht)*w2;
    const double b216 = +2.0 - (0.25*a1*a1*((ht*ht)/(0.5*hx*hx)))*w2 + ((0.25*c)*(ht*ht))*w2;
    const double b217 = +((0.25*a1*a1*ht*ht)/(0.5*hx*hx))*w3;
    const double b218 = -(0.25*a1*a1*(ht*ht)/(0.5*hx))*w3 - (0.25*b1*ht*ht)*w3;
    const double b219 = -1.0 - (0.25*a1*a1*((ht*ht)/(0.5*hx*hx)))*w3 + ((0.25*c)*(ht*ht))*w3 + 0.25*d*ht;
    const double b220 = +((0.25*a2*a2*ht*ht)/(hy*hy));
    const double b221 = +((0.25*b2*ht*ht)/(2.0*hy));
    const double b222 = +(0.25*a1*a1*ht*ht)/(0.5*hx) + (0.25*b1*ht*ht);

    // bottom border condition parameters y = ymin
    const double b311 = +1.0 + (0.25*a2*a2*((ht*ht)/(0.5*hy*hy)))*w1 - ((0.25*c)*(ht*ht))*w1 + 0.25*d*ht;
    const double b312 = +(0.25*a2*a2*(ht*ht)/(0.5*hy))*w1 - (0.25*b2*ht*ht)*w1;
    const double b313 = -((0.25*a2*a2*ht*ht)/(0.5*hy*hy))*w1;
    const double b314 = +2.0 - (0.25*a2*a2*((ht*ht)/(0.5*hy*hy)))*w2 + ((0.25*c)*(ht*ht))*w2;
    const double b315 = -(0.25*a2*a2*(ht*ht)/(0.5*hy))*w2 + (0.25*b2*ht*ht)*w2;
    const double b316 = +((0.25*a2*a2*ht*ht)/(0.5*hy*hy))*w2;
    const double b317 = -1.0 - (0.25*a2*a2*((ht*ht)/(0.5*hy*hy)))*w3 + ((0.25*c)*(ht*ht))*w3 + 0.25*d*ht;
    const double b318 = -(0.25*a2*a2*(ht*ht)/(0.5*hy))*w3 + (0.25*b2*ht*ht)*w3;
    const double b319 = +((0.25*a2*a2*ht*ht)/(0.5*hy*hy))*w3;
    const double b320 = +((0.25*a1*a1*ht*ht)/(hx*hx));
    const double b321 = +((0.25*b1*ht*ht)/(2.0*hx));
    const double b322 = +(0.25*a2*a2*ht*ht)/(0.5*hy) - (0.25*b2*ht*ht);

    // upper border condition parameters y = ymax
    const double b411 = -((0.25*a2*a2*ht*ht)/(0.5*hy*hy))*w1;
    const double b412 = +(0.25*a2*a2*(ht*ht)/(0.5*hy))*w1 + (0.25*b2*ht*ht)*w1;
    const double b413 = +1.0 + (0.25*a2*a2*((ht*ht)/(0.5*hy*hy)))*w1 - ((0.25*c)*(ht*ht))*w1 + 0.25*d*ht;
    const double b414 = +((0.25*a2*a2*ht*ht)/(0.5*hy*hy))*w2;
    const double b415 = -(0.25*a2*a2*(ht*ht)/(0.5*hy))*w2 - (0.25*b2*ht*ht)*w2;
    const double b416 = +2.0 - (0.25*a2*a2*((ht*ht)/(0.5*hy*hy)))*w2 + ((0.25*c)*(ht*ht))*w2;
    const double b417 = +((0.25*a2*a2*ht*ht)/(0.5*hy*hy))*w3;
    const double b418 = -(0.25*a2*a2*(ht*ht)/(0.5*hy))*w3 - (0.25*b1*ht*ht)*w3;
    const double b419 = -1.0 - (0.25*a2*a2*((ht*ht)/(0.5*hy*hy)))*w3 + ((0.25*c)*(ht*ht))*w3 + 0.25*d*ht;
    const double b420 = +((0.25*a1*a1*ht*ht)/(hx*hx));
    const double b421 = +((0.25*b1*ht*ht)/(2.0*hx));
    const double b422 = +(0.25*a2*a2*ht*ht)/(0.5*hy) + (0.25*b2*ht*ht);

    // initial condition parameters
    const double a1a1__hxhx = (a1*a1)/(hx*hx);
    const double a2a2__hyhy = (a2*a2)/(hy*hy);
    const double b1__20_0hx = ((b1)/(2.0*hx));
    const double b2__20_0hy = ((b2)/(2.0*hy));


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

    //    DoubleMatrix u00(M+1, N+1);
    //    DoubleMatrix u05(M+1, N+1);
    //    DoubleMatrix u10(M+1, N+1);
    //    DoubleMatrix u15(M+1, N+1);
    //    DoubleMatrix u20(M+1, N+1);

    double **u00 = static_cast<double**>(malloc(sizeof(double*)*(M+1)));
    double **u05 = static_cast<double**>(malloc(sizeof(double*)*(M+1)));
    double **u10 = static_cast<double**>(malloc(sizeof(double*)*(M+1)));
    double **u15 = static_cast<double**>(malloc(sizeof(double*)*(M+1)));
    double **u20 = static_cast<double**>(malloc(sizeof(double*)*(M+1)));
    for (unsigned int i=0; i<=M; i++)
    {
        u00[i] = static_cast<double*>(malloc(sizeof(double)*(N+1)));
        u05[i] = static_cast<double*>(malloc(sizeof(double)*(N+1)));
        u10[i] = static_cast<double*>(malloc(sizeof(double)*(N+1)));
        u15[i] = static_cast<double*>(malloc(sizeof(double)*(N+1)));
        u20[i] = static_cast<double*>(malloc(sizeof(double)*(N+1)));
    }

    /***********************************************************************************************/
    /***********************************************************************************************/

    TimeNodePDE tn00; tn00.i = 0; tn00.t = 0.5*tn00.i*ht;
    TimeNodePDE tn05; tn05.i = 1; tn05.t = 0.5*tn05.i*ht;
    TimeNodePDE tn10; tn10.i = 2; tn10.t = 0.5*tn10.i*ht;
    SpaceNodePDE sn;

    unsigned int i=0, j=0, s=0, e=N;
    int n = 0, m = 0;
    BoundaryConditionPDE condition; double value/*, alpha, beta, gamma*/;

    for (m=ymin, sn.j = m, sn.y = m*hy, j=0; m<=ymax; ++m, sn.j=m, sn.y=m*hy, ++j)
    {
        for (n=xmin, sn.i = n, sn.x = n*hx, i=0; n<=xmax; ++n, sn.i=n, sn.x=n*hx, ++i)
        {
            u00[j][i] = initial(sn, InitialCondition::InitialValue);
        }
    }

    layerInfo(DoubleMatrix(u00, M+1, N+1), tn00);

    /***********************************************************************************************/

    for (m=ymin, sn.j = m, sn.y = m*hy, j=0; m<=ymax; ++m, sn.j=m, sn.y=m*hy, ++j)
    {
        for (n=xmin, sn.i = n, sn.x = n*hx, i=0; n<=xmax; ++n, sn.i=n, sn.x=n*hx, ++i)
        {
            double fdt = initial(sn, InitialCondition::InitialFirstDerivative);

            double sdxy = 0.0;
            double fdxy = 0.0;
            if (j==0)
            {
                sdxy += a2a2__hyhy*(+2.0*u00[0][i]-5.0*u00[1][i]+4.0*u00[2][i]-1.0*u00[3][i]);
                fdxy += b2__20_0hy*(-3.0*u00[0][i]+4.0*u00[1][i]-1.0*u00[2][i]);
            }
            else if (j==M)
            {
                sdxy += a2a2__hyhy*(-1.0*u00[M-3][i]+4.0*u00[M-2][i]-5.0*u00[M-1][i]+2.0*u00[M][i]);
                fdxy += b2__20_0hy*(u00[M-2][i]-4.0*u00[M-1][i]+3.0*u00[M][i]);
            }
            else
            {
                sdxy += a2a2__hyhy*(u00[j-1][i]-2.0*u00[j][i]+u00[j+1][i]);
                fdxy += b2__20_0hy*(u00[j+1][i]-u00[j-1][i]);
            }

            if (i==0)
            {
                sdxy += a1a1__hxhx*(+2.0*u00[j][0]-5.0*u00[j][1]+4.0*u00[j][2]-1.0*u00[j][3]);
                fdxy += b1__20_0hx*(-3.0*u00[j][0]+4.0*u00[j][1]-1.0*u00[j][2]);
            }
            else if (i==N)
            {
                sdxy += a1a1__hxhx*(-1.0*u00[j][N-3]+4.0*u00[j][N-2]-5.0*u00[j][N-1]+2.0*u00[j][N]);
                fdxy += b1__20_0hx*(u00[j][N-2]-4.0*u00[j][N-1]+3.0*u00[j][N]);
            }
            else
            {
                sdxy += a1a1__hxhx*(u00[j][i-1]-2.0*u00[j][i]+u00[j][i+1]);
                fdxy += b1__20_0hx*(u00[j][i+1]-u00[j][i-1]);
            }

            double sdt = sdxy + fdxy + c*u00[j][i] + f(sn, tn00) - d*fdt;

            u05[j][i] = u00[j][i] + fdt*ht_0500 + sdt*htht_0125;
            u10[j][i] = u00[j][i] + fdt*ht      + sdt*htht_0500;
        }
    }

    //layerInfo(u05, tn05);
    //layerInfo(u10, tn10);
    layerInfo(DoubleMatrix(u05, M+1, N+1), tn00);
    layerInfo(DoubleMatrix(u10, M+1, N+1), tn00);

    /***********************************************************************************************/

    for (unsigned int ln=2; ln<=L; ln++)
    {
        //TimeNodePDE tn05; tn05.i = 2*ln-3; tn05.t = 0.5*tn05.i*ht;
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
                dx[i] += htht_0250 * f(sn, tn10);
            }

            sn.i = xmin; sn.x = xmin*hx;
            value = boundary(sn, tn15, condition);

            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                s = 1;
                u15[j][0] = /*(gamma/alpha)**/value;
                dx[1] -= k101 * u15[j][0];
                ax[1] = ax[0] = bx[0] = cx[0] = dx[0] = rx[0] = 0.0;
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
                s = 0;
                ax[0]  = 0.0;
                bx[0]  = /*beta  **/ b111;
                cx[0]  = /*beta  **/ b113;

                dx[0]  = /*beta  **/ b114 * u10[j][0] + /*beta  **/ b116 * u10[j][1];
                dx[0] += /*beta  **/ b117 * u05[j][0] + /*beta  **/ b119 * u05[j][1];

                dx[0] += /*beta  **/ b120 * (u10[j+1][0]-2.0*u10[j][0]+u10[j-1][0]);
                dx[0] += /*beta  **/ b121 * (u10[j+1][0]-u10[j-1][0]);

                dx[0] += /*beta  **/ htht_0250 * f(sn, tn10);
                dx[0] += /*gamma **/ b122 * boundary(sn, tn10, condition);
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                s = 0;

                double alpha = condition.alpha();
                double beta  = condition.beta();
                /*double gamma = condition.gamma();*/

                ax[0]  = 0.0;
                bx[0]  = beta  * b111 + alpha * b112;
                cx[0]  = beta  * b113;

                dx[0]  = beta  * b114 * u10[j][0] + alpha  * b115 * u10[j][0] + beta  * b116 * u10[j][1];
                dx[0] += beta  * b117 * u05[j][0] + alpha  * b118 * u05[j][0] + beta  * b119 * u05[j][1];

                dx[0] += beta  * b120 * (u10[j+1][0]-2.0*u10[j][0]+u10[j-1][0]);
                dx[0] += beta  * b121 * (u10[j+1][0]-u10[j-1][0]);

                dx[0] += beta  * htht_0250 * f(sn, tn10);
                dx[0] += /*gamma **/ b122 * boundary(sn, tn10, condition);
            }

            sn.i = xmax; sn.x = xmax*hx;
            value = boundary(sn, tn15, condition);

            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                e = N-1;
                u15[j][N] = /*(gamma/alpha)**/value;
                dx[N-1] -= k103 * u15[j][N];
                cx[N-1] = ax[N] = bx[N] = cx[N] = dx[N] = rx[N] = 0.0;
            }
            if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
                e = N;
                ax[N]  = /*beta  **/ b211;
                bx[N]  = /*beta  **/ b213;
                cx[N]  = 0.0;

                dx[N]  = /*beta  **/ b214 * u10[j][N-1] + /*beta  **/ b216 * u10[j][N];
                dx[N] += /*beta  **/ b217 * u05[j][N-1] + /*beta  **/ b219 * u05[j][N];

                dx[N] += /*beta  **/ b220 * (u10[j+1][N]-2.0*u10[j][N]+u10[j-1][N]);
                dx[N] += /*beta  **/ b221 * (u10[j+1][N]-u10[j-1][N]);

                dx[N] += /*beta  **/ htht_0250 * f(sn, tn10);
                dx[N] += /*gamma **/ b222 * boundary(sn, tn10, condition);
            }
            if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                e = N;

                double alpha = condition.alpha();
                double beta  = condition.beta();
                /*double gamma = condition.gamma();*/

                ax[N]  = beta  * b211;
                bx[N]  = beta  * b213 + alpha * b212;
                cx[N]  = 0.0;

                dx[N]  = beta  * b214 * u10[j][N-1] + alpha  * b215 * u10[j][N] + beta  * b216 * u10[j][N];
                dx[N] += beta  * b217 * u05[j][N-1] + alpha  * b218 * u05[j][N] + beta  * b219 * u05[j][N];

                dx[N] += beta  * b220 * (u10[j+1][N]-2.0*u10[j][N]+u10[j-1][N]);
                dx[N] += beta  * b221 * (u10[j+1][N]-u10[j-1][N]);

                dx[N] += beta  * htht_0250 * f(sn, tn10);
                dx[N] += /*gamma **/ b222 * boundary(sn, tn10, condition);
            }

            tomasAlgorithmLeft2Right(ax+s, bx+s, cx+s, dx+s, rx+s, e-s+1);
            for (unsigned int i=s; i<=e; i++) u15[j][i] = rx[i];
        }

        for (n=xmin, sn.i=n, sn.x=n*hx, i=0; n<=xmax; ++n, sn.i=n, sn.x=n*hx, ++i)
        {
            sn.j = ymin; sn.y = ymin*hy;
            value = boundary(sn, tn15, condition);

            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                u15[0][i] = /*(gamma/alpha)**/value;
            }
            if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
                u15[0][i] = (3.5*u15[1][i] - 2.0*u15[2][i] + 0.5*u15[3][i] + hy*/*(gamma/beta)**/value)/(2.0);
                //u05[0][i] = (u05[1][i] + hy*(gamma/beta)*value));
            }
            if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                double alpha = condition.alpha();
                double beta  = condition.beta();
                /*gamma = condition.gamma()*/;

                u15[0][i] = (3.5*u15[1][i] - 2.0*u15[2][i] + 0.5*u15[3][i] + hy*/*(gamma/beta)**/value)/(2.0 + (alpha/beta)*hy);
                //u05[0][i] = (u05[1][i] + hy*(gamma/beta)*value)/(1.0+hy*(alpha/beta));
            }

            sn.j = ymax; sn.y = ymax*hy;
            value = boundary(sn, tn15, condition);

            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                u15[M][i] = /*(gamma/alpha)**/value;
            }
            if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
                u15[M][i] = (3.5*u15[M-1][i] - 2.0*u15[M-2][i] + 0.5*u15[M-3][i] + hy*/*(gamma/beta)**/value)/(2.0);
                //u15[M][i] = (u15[M-1][i] + hy*(gamma/beta)*value));
            }
            if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                double alpha = condition.alpha();
                double beta  = condition.beta();
                /*gamma = condition.gamma();*/

                u15[M][i] = (3.5*u15[M-1][i] - 2.0*u15[M-2][i] + 0.5*u15[M-3][i] + hy*(/*gamma*/1.0/beta)*value)/(2.0 + (alpha/beta)*hy);
                //u15[M][i] = (u15[M-1][i] + hy*(gamma/beta)*value)/(1.0 + hy*(alpha/beta));
            }
        }

        //layerInfo(u15, tn15);
        layerInfo(DoubleMatrix(u15, M+1, N+1), tn00);

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
                //dy[j] += htht_0250 * f(sn, tn15);
                dy[j] += htht_0250 * (0.5*(f(sn, tn10)+f(sn, tn20)));
            }

            sn.j = ymin; sn.y = ymin*hy;
            value = boundary(sn, tn20, condition);

            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                s = 1;
                u20[0][i] = /*(gamma/alpha)**/value;
                dy[1] -= k201 * u20[0][i];
                ay[1] = ay[0] = by[0] = cy[0] = dy[0] = ry[0] = 0.0;
            }
            if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
                s = 0;
                ay[0]  = 0.0;
                by[0]  = /*beta  **/ b311;
                cy[0]  = /*beta  **/ b313;

                dy[0]  = /*beta  **/ b314 * u15[0][i] + /*beta  **/ b316 * u15[1][i];
                dy[0] += /*beta  **/ b317 * u10[0][i] + /*beta  **/ b319 * u10[1][i];

                dy[0] += /*beta  **/ b320 * (u15[0][i+1]-2.0*u15[0][i]+u15[0][i-1]);
                dy[0] += /*beta  **/ b321 * (u15[0][i+1]-u15[0][i-1]);

                dy[0] += /*beta  **/ htht_0250 * f(sn, tn15);
                dy[0] += /*gamma **/ b322 * boundary(sn, tn15, condition);
            }
            if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                s = 0;

                double alpha = condition.alpha();
                double beta  = condition.beta();
                /*gamma = condition.gamma();*/

                ay[0]  = 0.0;
                by[0]  = /*beta  **/ b311 + alpha * b312;
                cy[0]  = /*beta  **/ b313;

                dy[0]  = /*beta  **/ b314 * u15[0][i] + alpha  * b315 * u15[0][i] + beta  * b316 * u15[1][i];
                dy[0] += /*beta  **/ b317 * u10[0][i] + alpha  * b318 * u10[0][i] + beta  * b319 * u10[1][i];

                dy[0] += /*beta  **/ b320 * (u15[0][i+1]-2.0*u15[0][i]+u15[0][i-1]);
                dy[0] += /*beta  **/ b321 * (u15[0][i+1]-u15[0][i-1]);

                dy[0] += /*beta  **/ htht_0250 * f(sn, tn15);
                dy[0] += /*gamma **/ b322 * boundary(sn, tn15, condition);
            }

            sn.j = ymax; sn.y = ymax*hy;
            value = boundary(sn, tn20, condition);

            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                e = M-1;
                u20[M][i] = /*(gamma/alpha)**/value;
                dy[M-1] -= k203 * u20[M][i];
                cy[M-1] = ay[M] = by[M] = cy[M] = dy[M] = ry[M] = 0.0;
            }
            if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
                e = M;
                ay[M]  = /*beta  **/ b411;
                by[M]  = /*beta  **/ b413;
                cy[M]  = 0.0;

                dy[M]  = /*beta  **/ b414 * u15[M-1][i] + /*beta  **/ b416 * u15[M][i];
                dy[M] += /*beta  **/ b417 * u10[M-1][i] + /*beta  **/ b419 * u10[M][i];

                dy[M] += /*beta  **/ b420 * (u15[M][i+1]-2.0*u15[M][i]+u15[M][i-1]);
                dy[M] += /*beta  **/ b421 * (u15[M][i+1]-u15[M][i-1]);

                dy[M] += /*beta  **/ htht_0250 * f(sn, tn15);
                dy[M] += /*gamma **/ b422 * boundary(sn, tn15, condition);
            }
            if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                e = M;

                double alpha = condition.alpha();
                double beta  = condition.beta();
                /*gamma = condition.gamma();*/

                ay[M]  = beta  * b411;
                by[M]  = beta  * b413 + alpha * b412;
                cy[M]  = 0.0;

                dy[M]  = beta  * b414 * u15[M-1][i] + alpha  * b415 * u15[M][i] + beta  * b416 * u15[M][i];
                dy[M] += beta  * b417 * u10[M-1][i] + alpha  * b418 * u10[M][i] + beta  * b419 * u10[M][i];

                dy[M] += beta  * b420 * (u15[M][i+1]-2.0*u15[M][i]+u15[M][i-1]);
                dy[M] += beta  * b421 * (u15[M][i+1]-u15[M][i-1]);

                dy[M] += beta  * htht_0250 * f(sn, tn15);
                dy[M] += /*gamma **/ b422 * boundary(sn, tn15, condition);
            }

            tomasAlgorithmLeft2Right(ay+s, by+s, cy+s, dy+s, ry+s, e-s+1);
            for (unsigned int j=s; j<=e; j++) u20[j][i] = ry[j];
        }

        for (m=ymin, sn.j=m, sn.y=m*hy, j=0; m<=ymax; ++m, sn.j=m, sn.y=m*hy, ++j)
        {
            sn.i = xmin; sn.x = xmin*hx;
            value = boundary(sn, tn20, condition);

            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                u20[j][0] = /*(gamma/alpha)**/value;
            }
            if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
                u20[j][0] = (3.5*u20[j][1] - 2.0*u20[j][2] + 0.5*u20[j][3] + hx*/*(gamma/beta)**/value)/(2.0);
                //u20[j][0] = (u20[j][1] + hx*(gamma/beta)*value);
            }
            if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                double alpha = condition.alpha();
                double beta  = condition.beta();
                /*gamma = condition.gamma();*/

                u20[j][0] = (3.5*u20[j][1] - 2.0*u20[j][2] + 0.5*u20[j][3] + hx*(/*gamma*/1.0/beta)*value)/(2.0 + (alpha/beta)*hx);
                //u20[j][0] = (u20[j][1] + hx*(gamma/beta)*value)/(1.0+hx*(alpha/beta));
            }

            sn.i = xmax; sn.x = xmax*hx;
            value = boundary(sn, tn20, condition);

            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                u20[j][N] = /*(gamma/alpha)**/value;
            }
            if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
                u20[j][N] = (3.5*u20[j][N-1] - 2.0*u20[j][N-2] + 0.5*u20[j][N-3] + hx*/*(gamma/beta)**/value)/(2.0);
                //u20[j][N] = (u20[j][N-1] + hx*(gamma/beta)*value);
            }
            if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                double alpha = condition.alpha();
                double beta  = condition.beta();
                /*gamma = condition.gamma();*/

                u20[j][N] = (3.5*u20[j][N-1] - 2.0*u20[j][N-2] + 0.5*u20[j][N-3] + hx*(/*gamma*/1.0/beta)*value)/(2.0 + (alpha/beta)*hx);
                //u20[j][N] = (u20[j][N-1] + hx*(gamma/beta)*value)/(1.0 + hx*(alpha/beta));
            }
        }

        //layerInfo(u20, tn20);
        layerInfo(DoubleMatrix(u20, M+1, N+1), tn00);

        /**************************************************** y direction apprx ***************************************************/

        //        for (unsigned int m=0; m<=M; m++)
        //        {
        //            for (unsigned int n=0; n<=N; n++)
        //            {
        //                u00[m][n] = u10[m][n];
        //                u05[m][n] = u15[m][n];
        //                u10[m][n] = u20[m][n];
        //            }
        //        }

        u00 = u10; u05 = u15; u10 = u20;
    }

    for (unsigned int i=0; i<=M; i++) { free(u00[i]); free(u05[i]); free(u10[i]); free(u15[i]); free(u20[i]); }
    free(u00); free(u05); free(u10); free(u15); free(u20);

    //    u00.clear();
    //    u05.clear();
    //    u10.clear();
    //    u15.clear();
    //    u20.clear();

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
    const double b__2_0hx = (b/(2.0*hx));

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

    DoubleVector p00(N+1);
    DoubleVector p10(N+1);
    DoubleVector p20(N+1);

    /***********************************************************************************************/
    /***********************************************************************************************/

    TimeNodePDE tn00; tn00.i = L;   tn00.t = tn00.i*ht;
    TimeNodePDE tn10; tn10.i = L-1; tn10.t = tn10.i*ht;

    SpaceNodePDE sn;
    unsigned int i = 0;
    for (int n=xmin; n<=xmax; n++, i++)
    {
        sn.i = n; sn.x = n*hx;
        p00[i] = final(sn, FinalCondition::FinalValue);
    }
    layerInfo(p00, tn00);

    /***********************************************************************************************/

    i = 0;
    for (int n=xmin; n<=xmax; n++, i++)
    {
        sn.i = n; sn.x = n*hx;
        double fdt = final(sn, FinalCondition::FinalFirstDerivative);

        double sdx = 0.0;
        double fdx = 0.0;
        if (n==xmin)
        {
            sdx = aa__hxhx*(+2.0*p00[0]-5.0*p00[1]+4.0*p00[2]-1.0*p00[3]);
            fdx = b__2_0hx*(-3.0*p00[0]+4.0*p00[1]-p00[2]);
        }
        else if (n==xmax)
        {
            sdx = aa__hxhx*(-1.0*p00[N-3]+4.0*p00[N-2] -5.0*p00[N-1] +2.0*p00[N]);
            fdx = b__2_0hx*(p00[N-2]-4.0*p00[N-1]+3.0*p00[N]);
        }
        else
        {
            sdx = aa__hxhx*(p00[i-1]-2.0*p00[i]+p00[i+1]);
            fdx = b__2_0hx*(p00[i+1]-p00[i-1]);
        }

        double sdt = sdx + fdx + c*p00[i] + f(sn, tn00) - d*fdt;

        p10[i] = p00[i] - fdt*ht + sdt*htht_05;
    }

    layerInfo(p10, tn10);

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
            dx[i] += k21 * p10[i-1] + k22 * p10[i] + k23 * p10[i+1];
            dx[i] += k31 * p00[i-1] + k32 * p00[i] + k33 * p00[i+1];

            dx[i] += ht_ht*f(sn, tn10);
            //dx[i] += ht_ht*w1*f(sn, tn20) + ht_ht*w2*f(sn, tn10) + ht_ht*w3*f(sn, tn00);
        }

        unsigned int s=0, e=N;
        BoundaryConditionPDE condition; double alpha, beta, gamma, value;

        sn.i = xmin; sn.x = xmin*hx;
        value = boundary(sn, tn20, condition);

        if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
        {
            s = 1;

            p20[0] = /*(gamma/alpha)**/value;
            dx[1] -= k11 * p20[0];
            ax[1] = ax[0] = bx[0] = cx[0] = dx[0] = rx[0] = 0.0;
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
        {
            s = 0;

            ax[s]  = 0.0;
            bx[s]  = beta * b11 + alpha * b12;
            cx[s]  = beta * b13;

            dx[s]  = beta * b14 * p10[s] + beta * b16 * p10[s+1];
            dx[s] += beta * b17 * p00[s] + beta * b19 * p00[s+1];

            dx[s] += gamma * b20 * boundary(sn, tn10, condition);
            dx[s] += beta  * ht_ht*f(sn,tn10);

            //dx[s] += gamma * b20 * (w1*value + w2*boundary(sn, tn10, condition) + w3*boundary(sn, tn00, condition));
            //dx[s] += beta * (ht_ht*w1*f(sn, tn20) + ht_ht*w2*f(sn, tn10) + ht_ht*w3*f(sn, tn00));
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Robin)
        {
            s = 0;

            double alpha = condition.alpha();
            double beta  = condition.beta();
            /*gamma = condition.gamma();*/

            //ax[s] = 0.0;
            //bx[s] = hx*alpha-beta;
            //cx[s] = beta;
            //dx[s] = hx*gamma*value;

            ax[s]  = 0.0;
            bx[s]  = beta * b11 + alpha * b12;
            cx[s]  = beta * b13;

            dx[s]  = beta * b14 * p10[s] + alpha * b15 * p10[s] + beta * b16 * p10[s+1];
            dx[s] += beta * b17 * p00[s] + alpha * b18 * p00[s] + beta * b19 * p00[s+1];

            dx[s] += gamma * b20 * boundary(sn, tn10, condition);
            dx[s] += beta  * ht_ht*f(sn,tn10);

            //dx[s] += gamma * b20 * (w1*value + w2*boundary(sn, tn10, condition) + w3*boundary(sn, tn00, condition));
            //dx[s] += beta * (ht_ht*w1*f(sn, tn20) + ht_ht*w2*f(sn, tn10) + ht_ht*w3*f(sn, tn00));
        }

        sn.i = xmax; sn.x = xmax*hx;
        value = boundary(sn, tn20, condition);

        if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
        {
            e = N-1;

            p20[N] = (gamma/alpha)*value;
            dx[N-1] -= k13 * p20[N];
            cx[N-1] = ax[N] = bx[N] = cx[N] = dx[N] = rx[N] = 0.0;
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
        {
            e = N;

            ax[e]  = beta * b21;
            bx[e]  = beta * b23 + alpha * b22;
            cx[e]  = 0.0;

            dx[e]  = beta * b24 * p10[e-1] + beta * b26 * p10[e];
            dx[e] += beta * b27 * p00[e-1] + beta * b29 * p00[e];

            dx[e] += gamma * b30 * boundary(sn, tn10, condition);
            dx[e] += beta  * ht_ht*f(sn,tn10);

            //dx[e] += gamma * b30 * (w1*value + w2*boundary(sn, tn10, condition) + w3*boundary(sn, tn00, condition));
            //dx[e] += beta * (ht_ht*w1*f(sn, tn20) + ht_ht*w2*f(sn, tn10) + ht_ht*w3*f(sn, tn00));
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Robin)
        {
            e = N;

            alpha = condition.alpha();
            beta  = condition.beta();
            /*gamma = condition.gamma();*/

            //ax[e] = -beta;
            //bx[e] = hx*alpha+beta;
            //cx[e] = 0.0;
            //dx[e] = hx*gamma*value;

            ax[e]  = beta * b21;
            bx[e]  = beta * b23 + alpha * b22;
            cx[e]  = 0.0;

            dx[e]  = beta * b24 * p10[e-1] + alpha * b25 * p10[e] + beta * b26 * p10[e];
            dx[e] += beta * b27 * p00[e-1] + alpha * b28 * p00[e] + beta * b29 * p00[e];

            dx[e] += /*gamma **/ b30 * boundary(sn, tn10, condition);
            dx[e] += beta  * ht_ht*f(sn,tn10);

            //dx[e] += gamma * b30 * (w1*value + w2*boundary(sn, tn10, condition) + w3*boundary(sn, tn00, condition));
            //dx[e] += beta * (ht_ht*w1*f(sn, tn20) + ht_ht*w2*f(sn, tn10) + ht_ht*w3*f(sn, tn00));
        }

        tomasAlgorithmLeft2Right(ax+s, bx+s, cx+s, dx+s, rx+s, e-s+1);
        for (unsigned int n=s; n<=e; n++) p20[n] = rx[n];
        layerInfo(p20, tn20);

        for (unsigned int i=0; i<=N; i++)
        {
            p00[i] = p10[i];
            p10[i] = p20[i];
        }
    }

    p00.clear();
    p10.clear();
    p20.clear();

    free(rx);
    free(dx);
    free(cx);
    free(bx);
    free(ax);
}

void IWaveEquationFBVP::explicit_calculate_D2V1() const
{
    const unsigned int N = static_cast<unsigned int>(spaceDimensionX().size()) - 1;
    const unsigned int M = static_cast<unsigned int>(spaceDimensionY().size()) - 1;
    const unsigned int L = static_cast<unsigned int>(timeDimension().size()) - 1;

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

    double ht_max = 1.0/sqrt((a1*a1)/(hx*hx)+(a2*a2)/(hy*hy));
    if (ht > ht_max) { throw std::runtime_error("error"); }

    // common parameters
    const double htht10 = ht*ht;
    const double htht05 = ht*ht*0.5;

    double k1 = +1.0 - 0.5*d*ht;
    double k2 = -1.0 - 0.5*d*ht;
    double k3 = (a1*a1)*((ht*ht)/(hx*hx)) - b1*((ht*ht)/(2.0*hx));
    double k4 = (a1*a1)*((ht*ht)/(hx*hx)) + b1*((ht*ht)/(2.0*hx));
    double k5 = (a2*a2)*((ht*ht)/(hy*hy)) - b2*((ht*ht)/(2.0*hy));
    double k6 = (a2*a2)*((ht*ht)/(hy*hy)) + b2*((ht*ht)/(2.0*hy));
    double k7 = (-2.0*a1*a1)*((ht*ht)/(hx*hx)) + (-2.0*a2*a2)*((ht*ht)/(hy*hy)) + c*ht*ht + 2.0;

    // initial condition parameters
    const double a1a1__hxhx = (a1*a1)/(hx*hx);
    const double a2a2__hyhy = (a2*a2)/(hy*hy);
    const double b1__20_0hx = ((b1)/(2.0*hx));
    const double b2__20_0hy = ((b2)/(2.0*hy));

    //DoubleMatrix p00(M+1, N+1);
    //DoubleMatrix p10(M+1, N+1);
    //DoubleMatrix p20(M+1, N+1);

    double **p00 = static_cast<double**>(malloc(sizeof(double*)*(M+1)));
    double **p10 = static_cast<double**>(malloc(sizeof(double*)*(M+1)));
    double **p20 = static_cast<double**>(malloc(sizeof(double*)*(M+1)));
    for (unsigned int i=0; i<=M; i++)
    {
        p00[i] = static_cast<double*>(malloc(sizeof(double)*(N+1)));
        p10[i] = static_cast<double*>(malloc(sizeof(double)*(N+1)));
        p20[i] = static_cast<double*>(malloc(sizeof(double)*(N+1)));
    }

    auto calculate_border = [](SpaceNodePDE sn, TimeNodePDE tn, double **p, unsigned int i, unsigned int j, const IWaveEquationFBVP* w)
    {
        BoundaryConditionPDE condition; double value, alpha, beta, gamma;

        value = w->boundary(sn, tn, condition);

        if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
        {
            p[j][i] = /*(gamma/alpha)**/value;
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
        {
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Robin)
        {
            alpha = condition.alpha();
            beta  = condition.beta();
            /*gamma = condition.gamma();*/
        }
    };

    /***********************************************************************************************/
    /***********************************************************************************************/

    TimeNodePDE tn00; tn00.i = 2*(L-0); tn00.t = 0.5*tn00.i*ht;
    TimeNodePDE tn10; tn10.i = 2*(L-1); tn10.t = 0.5*tn10.i*ht;
    SpaceNodePDE sn;

    unsigned int i=0, j=0;
    int n = 0, m = 0;
    //BoundaryConditionPDE condition; double value, alpha, beta, gamma;

    for (m=ymin, sn.j = m, sn.y = m*hy, j=0; m<=ymax; ++m, sn.j=m, sn.y=m*hy, ++j)
    {
        for (n=xmin, sn.i = n, sn.x = n*hx, i=0; n<=xmax; ++n, sn.i=n, sn.x=n*hx, ++i)
        {
            p00[j][i] = final(sn, FinalCondition::FinalValue);
        }
    }

    layerInfo(DoubleMatrix(p00, M+1, N+1), tn00);

    /***********************************************************************************************/

    for (m=ymin, sn.j = m, sn.y = m*hy, j=0; m<=ymax; ++m, sn.j=m, sn.y=m*hy, ++j)
    {
        for (n=xmin, sn.i = n, sn.x = n*hx, i=0; n<=xmax; ++n, sn.i=n, sn.x=n*hx, ++i)
        {
            double fdt = final(sn, FinalCondition::FinalFirstDerivative);

            double sdxy = 0.0;
            double fdxy = 0.0;
            if (j==0)
            {
                sdxy += a2a2__hyhy*(+2.0*p00[0][i]-5.0*p00[1][i]+4.0*p00[2][i]-1.0*p00[3][i]);
                fdxy += b2__20_0hy*(-3.0*p00[0][i]+4.0*p00[1][i]-1.0*p00[2][i]);
            }
            else if (j==M)
            {
                sdxy += a2a2__hyhy*(-1.0*p00[M-3][i]+4.0*p00[M-2][i]-5.0*p00[M-1][i]+2.0*p00[M][i]);
                fdxy += b2__20_0hy*(p00[M-2][i]-4.0*p00[M-1][i]+3.0*p00[M][i]);
            }
            else
            {
                sdxy += a2a2__hyhy*(p00[j-1][i]-2.0*p00[j][i]+p00[j+1][i]);
                fdxy += b2__20_0hy*(p00[j+1][i]-p00[j-1][i]);
            }

            if (i==0)
            {
                sdxy += a1a1__hxhx*(+2.0*p00[j][0]-5.0*p00[j][1]+4.0*p00[j][2]-1.0*p00[j][3]);
                fdxy += b1__20_0hx*(-3.0*p00[j][0]+4.0*p00[j][1]-1.0*p00[j][2]);
            }
            else if (i==N)
            {
                sdxy += a1a1__hxhx*(-1.0*p00[j][N-3]+4.0*p00[j][N-2]-5.0*p00[j][N-1]+2.0*p00[j][N]);
                fdxy += b1__20_0hx*(p00[j][N-2]-4.0*p00[j][N-1]+3.0*p00[j][N]);
            }
            else
            {
                sdxy += a1a1__hxhx*(p00[j][i-1]-2.0*p00[j][i]+p00[j][i+1]);
                fdxy += b1__20_0hx*(p00[j][i+1]-p00[j][i-1]);
            }

            double sdt = sdxy + fdxy + c*p00[j][i] + f(sn,tn00) - d*fdt;
            p10[j][i] = p00[j][i] - fdt*ht + sdt*htht05;
        }
    }

    //layerInfo(p10, tn10);
    layerInfo(DoubleMatrix(p10, M+1, N+1), tn10);

    /***********************************************************************************************/

    for (unsigned int ln=L-2, size_ln = static_cast<unsigned int>(0)-1; ln != size_ln; ln--)
    {
        //TimeNodePDE tn00; tn00.i = ln+2; tn00.t = tn00.i*ht;
        TimeNodePDE tn10; tn10.i = 2*(ln+1); tn10.t = 0.5*tn10.i*ht;
        TimeNodePDE tn20; tn20.i = 2*(ln+0); tn20.t = 0.5*tn20.i*ht;

        for (m=ymin+1, sn.j=m, sn.y=m*hy, j=1; m<=ymax-1; ++m, sn.j=m, sn.y=m*hy, ++j)
        {
            for (n=xmin+1, sn.i=n, sn.x=n*hx, i=1; n<=xmax-1; ++n, sn.i=n, sn.x=n*hx, ++i)
            {
                p20[j][i] = k2*p00[j][i] + k3*p10[j][i-1] + k4*p10[j][i+1] + k5*p10[j-1][i] + k6*p10[j+1][i] + k7*p10[j][i] + htht10*f(sn, tn10);
                p20[j][i] /= k1;
            }
        }

        /**************************************************** border conditions ***************************************************/

        SpaceNodePDE sn1, sn2;
        sn1.j = ymin; sn1.y = ymin*hy;
        sn2.j = ymax; sn2.y = ymax*hy;
        for (n=xmin+1, sn1.i=sn2.i=n, sn1.x=sn2.x=n*hx, i=1; n<=xmax-1; ++n, sn1.i=sn2.i=n, sn1.x=sn2.x=n*hx, ++i)
        {
            calculate_border(sn1, tn20, p20, i, 0, this);
            calculate_border(sn2, tn20, p20, i, M, this);
        }

        sn1.i = xmin; sn1.x = xmin*hx;
        sn2.i = xmax; sn2.x = xmax*hx;
        for (m=ymin+1, sn1.j=sn2.j=m, sn1.y=sn2.y=m*hy, j=1; m<=ymax-1; ++m, sn1.j=sn2.j=m, sn1.y=sn2.y=m*hy, ++j)
        {
            calculate_border(sn1, tn20, p20, 0, j, this);
            calculate_border(sn2, tn20, p20, N, j, this);
        }

        sn.i = xmin; sn.x = xmin*hx; sn.j = ymin; sn.y = ymin*hy;
        calculate_border(sn, tn20, p20, 0, 0, this);
        sn.i = xmin; sn.x = xmin*hx; sn.j = ymax; sn.y = ymax*hy;
        calculate_border(sn, tn20, p20, 0, M, this);
        sn.i = xmax; sn.x = xmax*hx; sn.j = ymax; sn.y = ymax*hy;
        calculate_border(sn, tn20, p20, N, M, this);
        sn.i = xmax; sn.x = xmax*hx; sn.j = ymin; sn.y = ymin*hy;
        calculate_border(sn, tn20, p20, N, 0, this);

        /**************************************************** border conditions ***************************************************/

        layerInfo(DoubleMatrix(p20, M+1, N+1), tn20);

        //        for (unsigned int m=0; m<=M; m++)
        //        {
        //            for (unsigned int n=0; n<=N; n++)
        //            {
        //                p00[m][n] = p10[m][n];
        //                p10[m][n] = p20[m][n];
        //            }
        //        }

        double **_tmp = p00; p00 = p10; p10 = p20; p20 = _tmp;
    }

    for (unsigned int i=0; i<=M; i++) { free(p00[i]); free(p10[i]); free(p20[i]); }
    free(p00); free(p10); free(p20);

    //p00.clear();
    //p10.clear();
    //p20.clear();
}

void IWaveEquationFBVP::implicit_calculate_D2V1() const
{
    const unsigned int N = static_cast<unsigned int>(spaceDimensionX().size()) - 1;
    const unsigned int M = static_cast<unsigned int>(spaceDimensionY().size()) - 1;
    const unsigned int L = static_cast<unsigned int>(timeDimension().size()) - 1;

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
    const double htht_0500 = ht*ht*0.500;
    const double htht_0250 = ht*ht*0.250;
    const double htht_0125 = ht*ht*0.125;
    const double ht_0500 = ht*0.500;
    //const double ht_0250 = ht*0.250;

    // equation parameters
    const double k101 = -(0.25*a1*a1)*((ht*ht)/(hx*hx))*w1 + (0.125*b1)*((ht*ht)/(hx))*w1;           // i-1
    const double k102 = +1.0 + (0.5*a1*a1)*((ht*ht)/(hx*hx))*w1 - (0.25*c)*(ht*ht)*w1 - 0.25*d*ht;   // i
    const double k103 = -(0.25*a1*a1)*((ht*ht)/(hx*hx))*w1 - (0.125*b1)*((ht*ht)/(hx))*w1;           // i+1
    const double k104 = +(0.25*a1*a1)*((ht*ht)/(hx*hx))*w2 - (0.125*b1)*((ht*ht)/(hx))*w2;           // i-1
    const double k105 = +2.0 - (0.5*a1*a1)*((ht*ht)/(hx*hx))*w2 + (0.25*c)*(ht*ht)*w2;               // i
    const double k106 = +(0.25*a1*a1)*((ht*ht)/(hx*hx))*w2 + (0.125*b1)*((ht*ht)/(hx))*w2;           // i+1
    const double k107 = +(0.25*a1*a1)*((ht*ht)/(hx*hx))*w3 - (0.125*b1)*((ht*ht)/(hx))*w3;           // i-1
    const double k108 = -1.0 - (0.5*a1*a1)*((ht*ht)/(hx*hx))*w3 + (0.25*c)*(ht*ht)*w3 - 0.25*d*ht;   // i
    const double k109 = +(0.25*a1*a1)*((ht*ht)/(hx*hx))*w3 + (0.125*b1)*((ht*ht)/(hx))*w3;           // i+1
    const double k110 = +(0.25*a2*a2)*((ht*ht)/(hy*hy)) - (0.125*b2)*((ht*ht)/(hy));                 // j-1
    const double k111 = -(0.50*a2*a2)*((ht*ht)/(hy*hy));                                             // j
    const double k112 = +(0.25*a2*a2)*((ht*ht)/(hy*hy)) + (0.125*b2)*((ht*ht)/(hy));                 // j+1

    const double k201 = -(0.25*a2*a2)*((ht*ht)/(hy*hy))*w1 + (0.125*b2)*((ht*ht)/(hy))*w1;           // j-1
    const double k202 = +1.0 + (0.5*a2*a2)*((ht*ht)/(hy*hy))*w1 - (0.25*c)*(ht*ht)*w1 - 0.25*d*ht;   // j
    const double k203 = -(0.25*a2*a2)*((ht*ht)/(hy*hy))*w1 - (0.125*b2)*((ht*ht)/(hy))*w1;           // j+1
    const double k204 = +(0.25*a2*a2)*((ht*ht)/(hy*hy))*w2 - (0.125*b2)*((ht*ht)/(hy))*w2;           // j-1
    const double k205 = +2.0 - (0.5*a2*a2)*((ht*ht)/(hy*hy))*w2 + (0.25*c)*(ht*ht)*w2;               // j
    const double k206 = +(0.25*a2*a2)*((ht*ht)/(hy*hy))*w2 + (0.125*b2)*((ht*ht)/(hy))*w2;           // j+1
    const double k207 = +(0.25*a2*a2)*((ht*ht)/(hy*hy))*w3 - (0.125*b2)*((ht*ht)/(hy))*w3;           // j-1
    const double k208 = -1.0 - (0.5*a2*a2)*((ht*ht)/(hy*hy))*w3 + (0.25*c)*(ht*ht)*w3 - 0.25*d*ht;   // j
    const double k209 = +(0.25*a2*a2)*((ht*ht)/(hy*hy))*w3 + (0.125*b2)*((ht*ht)/(hy))*w3;           // j+1
    const double k210 = +(0.25*a1*a1)*((ht*ht)/(hx*hx)) - (0.125*b1)*((ht*ht)/(hx));                 // i-1
    const double k211 = -(0.50*a1*a1)*((ht*ht)/(hx*hx));                                             // i
    const double k212 = +(0.25*a1*a1)*((ht*ht)/(hx*hx)) + (0.125*b1)*((ht*ht)/(hx));                 // i+1

    // left border condition parameters x = xmin
    const double b111 = +1.0 + (0.25*a1*a1*((ht*ht)/(0.5*hx*hx)))*w1 - ((0.25*c)*(ht*ht))*w1 - 0.25*d*ht;
    const double b112 = +(0.25*a1*a1*(ht*ht)/(0.5*hx))*w1 - (0.25*b1*ht*ht)*w1;
    const double b113 = -((0.25*a1*a1*ht*ht)/(0.5*hx*hx))*w1;
    const double b114 = +2.0 - (0.25*a1*a1*((ht*ht)/(0.5*hx*hx)))*w2 + ((0.25*c)*(ht*ht))*w2;
    const double b115 = -(0.25*a1*a1*(ht*ht)/(0.5*hx))*w2 + (0.25*b1*ht*ht)*w2;
    const double b116 = +((0.25*a1*a1*ht*ht)/(0.5*hx*hx))*w2;
    const double b117 = -1.0 - (0.25*a1*a1*((ht*ht)/(0.5*hx*hx)))*w3 + ((0.25*c)*(ht*ht))*w3 - 0.25*d*ht;
    const double b118 = -(0.25*a1*a1*(ht*ht)/(0.5*hx))*w3 + (0.25*b1*ht*ht)*w3;
    const double b119 = +((0.25*a1*a1*ht*ht)/(0.5*hx*hx))*w3;
    const double b120 = +((0.25*a2*a2*ht*ht)/(hy*hy));
    const double b121 = +((0.25*b2*ht*ht)/(2.0*hy));
    const double b122 = +(0.25*a1*a1*ht*ht)/(0.5*hx) - (0.25*b1*ht*ht);

    // right border condition parameters x = xmax
    const double b211 = -((0.25*a1*a1*ht*ht)/(0.5*hx*hx))*w1;
    const double b212 = +(0.25*a1*a1*(ht*ht)/(0.5*hx))*w1 + (0.25*b1*ht*ht)*w1;
    const double b213 = +1.0 + (0.25*a1*a1*((ht*ht)/(0.5*hx*hx)))*w1 - ((0.25*c)*(ht*ht))*w1 - 0.25*d*ht;
    const double b214 = +((0.25*a1*a1*ht*ht)/(0.5*hx*hx))*w2;
    const double b215 = -(0.25*a1*a1*(ht*ht)/(0.5*hx))*w2 - (0.25*b1*ht*ht)*w2;
    const double b216 = +2.0 - (0.25*a1*a1*((ht*ht)/(0.5*hx*hx)))*w2 + ((0.25*c)*(ht*ht))*w2;
    const double b217 = +((0.25*a1*a1*ht*ht)/(0.5*hx*hx))*w3;
    const double b218 = -(0.25*a1*a1*(ht*ht)/(0.5*hx))*w3 - (0.25*b1*ht*ht)*w3;
    const double b219 = -1.0 - (0.25*a1*a1*((ht*ht)/(0.5*hx*hx)))*w3 + ((0.25*c)*(ht*ht))*w3 - 0.25*d*ht;
    const double b220 = +((0.25*a2*a2*ht*ht)/(hy*hy));
    const double b221 = +((0.25*b2*ht*ht)/(2.0*hy));
    const double b222 = +(0.25*a1*a1*ht*ht)/(0.5*hx) + (0.25*b1*ht*ht);

    // bottom border condition parameters y = ymin
    const double b311 = +1.0 + (0.25*a2*a2*((ht*ht)/(0.5*hy*hy)))*w1 - ((0.25*c)*(ht*ht))*w1 - 0.25*d*ht;
    const double b312 = +(0.25*a2*a2*(ht*ht)/(0.5*hy))*w1 - (0.25*b2*ht*ht)*w1;
    const double b313 = -((0.25*a2*a2*ht*ht)/(0.5*hy*hy))*w1;
    const double b314 = +2.0 - (0.25*a2*a2*((ht*ht)/(0.5*hy*hy)))*w2 + ((0.25*c)*(ht*ht))*w2;
    const double b315 = -(0.25*a2*a2*(ht*ht)/(0.5*hy))*w2 + (0.25*b2*ht*ht)*w2;
    const double b316 = +((0.25*a2*a2*ht*ht)/(0.5*hy*hy))*w2;
    const double b317 = -1.0 - (0.25*a2*a2*((ht*ht)/(0.5*hy*hy)))*w3 + ((0.25*c)*(ht*ht))*w3 - 0.25*d*ht;
    const double b318 = -(0.25*a2*a2*(ht*ht)/(0.5*hy))*w3 + (0.25*b2*ht*ht)*w3;
    const double b319 = +((0.25*a2*a2*ht*ht)/(0.5*hy*hy))*w3;
    const double b320 = +((0.25*a1*a1*ht*ht)/(hx*hx));
    const double b321 = +((0.25*b1*ht*ht)/(2.0*hx));
    const double b322 = +(0.25*a2*a2*ht*ht)/(0.5*hy) - (0.25*b2*ht*ht);

    // upper border condition parameters y = ymax
    const double b411 = -((0.25*a2*a2*ht*ht)/(0.5*hy*hy))*w1;
    const double b412 = +(0.25*a2*a2*(ht*ht)/(0.5*hy))*w1 + (0.25*b2*ht*ht)*w1;
    const double b413 = +1.0 + (0.25*a2*a2*((ht*ht)/(0.5*hy*hy)))*w1 - ((0.25*c)*(ht*ht))*w1 - 0.25*d*ht;
    const double b414 = +((0.25*a2*a2*ht*ht)/(0.5*hy*hy))*w2;
    const double b415 = -(0.25*a2*a2*(ht*ht)/(0.5*hy))*w2 - (0.25*b2*ht*ht)*w2;
    const double b416 = +2.0 - (0.25*a2*a2*((ht*ht)/(0.5*hy*hy)))*w2 + ((0.25*c)*(ht*ht))*w2;
    const double b417 = +((0.25*a2*a2*ht*ht)/(0.5*hy*hy))*w3;
    const double b418 = -(0.25*a2*a2*(ht*ht)/(0.5*hy))*w3 - (0.25*b1*ht*ht)*w3;
    const double b419 = -1.0 - (0.25*a2*a2*((ht*ht)/(0.5*hy*hy)))*w3 + ((0.25*c)*(ht*ht))*w3 - 0.25*d*ht;
    const double b420 = +((0.25*a1*a1*ht*ht)/(hx*hx));
    const double b421 = +((0.25*b1*ht*ht)/(2.0*hx));
    const double b422 = +(0.25*a2*a2*ht*ht)/(0.5*hy) + (0.25*b2*ht*ht);

    // initial condition parameters
    const double a1a1__hxhx = (a1*a1)/(hx*hx);
    const double a2a2__hyhy = (a2*a2)/(hy*hy);
    const double b1__20_0hx = ((b1)/(2.0*hx));
    const double b2__20_0hy = ((b2)/(2.0*hy));


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

    //DoubleMatrix p00(M+1, N+1);
    //DoubleMatrix p05(M+1, N+1);
    //DoubleMatrix p10(M+1, N+1);
    //DoubleMatrix p15(M+1, N+1);
    //DoubleMatrix p20(M+1, N+1);

    double **p00 = static_cast<double**>(malloc(sizeof(double*)*(M+1)));
    double **p05 = static_cast<double**>(malloc(sizeof(double*)*(M+1)));
    double **p10 = static_cast<double**>(malloc(sizeof(double*)*(M+1)));
    double **p15 = static_cast<double**>(malloc(sizeof(double*)*(M+1)));
    double **p20 = static_cast<double**>(malloc(sizeof(double*)*(M+1)));
    for (unsigned int i=0; i<=M; i++)
    {
        p00[i] = static_cast<double*>(malloc(sizeof(double)*(N+1)));
        p05[i] = static_cast<double*>(malloc(sizeof(double)*(N+1)));
        p10[i] = static_cast<double*>(malloc(sizeof(double)*(N+1)));
        p15[i] = static_cast<double*>(malloc(sizeof(double)*(N+1)));
        p20[i] = static_cast<double*>(malloc(sizeof(double)*(N+1)));
    }

    /***********************************************************************************************/
    /***********************************************************************************************/

    TimeNodePDE tn00; tn00.i = 2*L-0; tn00.t = 0.5*tn00.i*ht;
    TimeNodePDE tn05; tn05.i = 2*L-1; tn05.t = 0.5*tn05.i*ht;
    TimeNodePDE tn10; tn10.i = 2*L-2; tn10.t = 0.5*tn10.i*ht;
    SpaceNodePDE sn;

    unsigned int i=0, j=0, s=0, e=N;
    int n = 0, m = 0;
    BoundaryConditionPDE condition; double value/*, alpha, beta, gamma*/;

    for (m=ymin, sn.j = m, sn.y = m*hy, j=0; m<=ymax; ++m, sn.j=m, sn.y=m*hy, ++j)
    {
        for (n=xmin, sn.i = n, sn.x = n*hx, i=0; n<=xmax; ++n, sn.i=n, sn.x=n*hx, ++i)
        {
            p00[j][i] = final(sn, FinalCondition::FinalValue);
        }
    }
    layerInfo(DoubleMatrix(p00, M+1, N+1), tn00);

    /***********************************************************************************************/

    for (m=ymin, sn.j = m, sn.y = m*hy, j=0; m<=ymax; ++m, sn.j=m, sn.y=m*hy, ++j)
    {
        for (n=xmin, sn.i = n, sn.x = n*hx, i=0; n<=xmax; ++n, sn.i=n, sn.x=n*hx, ++i)
        {
            double fdt = final(sn, FinalCondition::FinalFirstDerivative);

            double sdxy = 0.0;
            double fdxy = 0.0;
            if (j==0)
            {
                sdxy += a2a2__hyhy*(+2.0*p00[0][i]-5.0*p00[1][i]+4.0*p00[2][i]-1.0*p00[3][i]);
                fdxy += b2__20_0hy*(-3.0*p00[0][i]+4.0*p00[1][i]-1.0*p00[2][i]);
            }
            else if (j==M)
            {
                sdxy += a2a2__hyhy*(-1.0*p00[M-3][i]+4.0*p00[M-2][i]-5.0*p00[M-1][i]+2.0*p00[M][i]);
                fdxy += b2__20_0hy*(p00[M-2][i]-4.0*p00[M-1][i]+3.0*p00[M][i]);
            }
            else
            {
                sdxy += a2a2__hyhy*(p00[j-1][i]-2.0*p00[j][i]+p00[j+1][i]);
                fdxy += b2__20_0hy*(p00[j+1][i]-p00[j-1][i]);
            }

            if (i==0)
            {
                sdxy += a1a1__hxhx*(+2.0*p00[j][0]-5.0*p00[j][1]+4.0*p00[j][2]-1.0*p00[j][3]);
                fdxy += b1__20_0hx*(-3.0*p00[j][0]+4.0*p00[j][1]-1.0*p00[j][2]);
            }
            else if (i==N)
            {
                sdxy += a1a1__hxhx*(-1.0*p00[j][N-3]+4.0*p00[j][N-2]-5.0*p00[j][N-1]+2.0*p00[j][N]);
                fdxy += b1__20_0hx*(p00[j][N-2]-4.0*p00[j][N-1]+3.0*p00[j][N]);
            }
            else
            {
                sdxy += a1a1__hxhx*(p00[j][i-1]-2.0*p00[j][i]+p00[j][i+1]);
                fdxy += b1__20_0hx*(p00[j][i+1]-p00[j][i-1]);
            }

            double sdt = sdxy + fdxy + c*p00[j][i] + f(sn,tn00) - d*fdt;

            p05[j][i] = p00[j][i] - fdt*ht_0500 + sdt*htht_0125;
            p10[j][i] = p00[j][i] - fdt*ht      + sdt*htht_0500;
        }
    }

    //layerInfo(p05, tn05);
    //layerInfo(p10, tn10);
    layerInfo(DoubleMatrix(p05, M+1, N+1), tn00);
    layerInfo(DoubleMatrix(p10, M+1, N+1), tn00);

    /***********************************************************************************************/

    for (unsigned int ln=L-2, size_ln = static_cast<unsigned int>(0)-1; ln !=size_ln; ln--)
    {
        //TimeNodePDE tn05; tn05.i = 2*ln+3; tn05.t = 0.5*tn05.i*ht;
        TimeNodePDE tn10; tn10.i = 2*ln+2; tn10.t = 0.5*tn10.i*ht;
        TimeNodePDE tn15; tn15.i = 2*ln+1; tn15.t = 0.5*tn15.i*ht;
        TimeNodePDE tn20; tn20.i = 2*ln+0; tn20.t = 0.5*tn20.i*ht;

        /**************************************************** x direction apprx ***************************************************/

        for (m=ymin+1, sn.j=m, sn.y=m*hy, j=1; m<=ymax-1; ++m, sn.j=m, sn.y=m*hy, ++j)
        {
            for (n=xmin+1, sn.i=n, sn.x=n*hx, i=1; n<=xmax-1; ++n, sn.i=n, sn.x=n*hx, ++i)
            {
                dx[i] = 0.0;
                dx[i] += k104*p10[j][i-1] + k105*p10[j][i] + k106*p10[j][i+1];
                dx[i] += k107*p05[j][i-1] + k108*p05[j][i] + k109*p05[j][i+1];
                dx[i] += k110*p10[j-1][i] + k111*p10[j][i] + k112*p10[j+1][i];
                dx[i] += htht_0250 * f(sn, tn10);
            }

            sn.i = xmin; sn.x = xmin*hx;
            value = boundary(sn, tn15, condition);

            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                s = 1;
                p15[j][0] = /*(gamma/alpha)**/value;
                dx[1] -= k101 * p15[j][0];
                ax[1] = ax[0] = bx[0] = cx[0] = dx[0] = rx[0] = 0.0;
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
                s = 0;
                ax[0]  = 0.0;
                bx[0]  = /*beta  **/ b111;
                cx[0]  = /*beta  **/ b113;

                dx[0]  = /*beta  **/ b114 * p10[j][0] + /*beta  **/ b116 * p10[j][1];
                dx[0] += /*beta  **/ b117 * p05[j][0] + /*beta  **/ b119 * p05[j][1];

                dx[0] += /*beta  **/ b120 * (p10[j+1][0]-2.0*p10[j][0]+p10[j-1][0]);
                dx[0] += /*beta  **/ b121 * (p10[j+1][0]-p10[j-1][0]);

                dx[0] += /*beta  **/ htht_0250 * f(sn, tn10);
                dx[0] += /*gamma **/ b122 * boundary(sn, tn10, condition);
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                s = 0;

                const double alpha = condition.alpha();
                const double beta  = condition.beta();
                /*gamma = condition.gamma();*/
                const double gamma = 1.0;

                ax[0]  = 0.0;
                bx[0]  = beta  * b111 + alpha * b112;
                cx[0]  = beta  * b113;

                dx[0]  = beta  * b114 * p10[j][0] + alpha  * b115 * p10[j][0] + beta  * b116 * p10[j][1];
                dx[0] += beta  * b117 * p05[j][0] + alpha  * b118 * p05[j][0] + beta  * b119 * p05[j][1];

                dx[0] += beta  * b120 * (p10[j+1][0]-2.0*p10[j][0]+p10[j-1][0]);
                dx[0] += beta  * b121 * (p10[j+1][0]-p10[j-1][0]);

                dx[0] += beta  * htht_0250 * f(sn, tn10);
                dx[0] += gamma * b122 * boundary(sn, tn10, condition);
            }

            sn.i = xmax; sn.x = xmax*hx;
            value = boundary(sn, tn15, condition);

            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                e = N-1;
                p15[j][N] = /*(gamma/alpha)**/value;
                dx[N-1] -= k103 * p15[j][N];
                cx[N-1] = ax[N] = bx[N] = cx[N] = dx[N] = rx[N] = 0.0;
            }
            if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
                e = N;
                ax[N]  = /*beta  **/ b211;
                bx[N]  = /*beta  **/ b213;
                cx[N]  = 0.0;

                dx[N]  = /*beta  **/ b214 * p10[j][N-1] + /*beta  **/ b216 * p10[j][N];
                dx[N] += /*beta  **/ b217 * p05[j][N-1] + /*beta  **/ b219 * p05[j][N];

                dx[N] += /*beta  **/ b220 * (p10[j+1][N]-2.0*p10[j][N]+p10[j-1][N]);
                dx[N] += /*beta  **/ b221 * (p10[j+1][N]-p10[j-1][N]);

                dx[N] += /*beta  **/ htht_0250 * f(sn, tn10);
                dx[N] += /*gamma **/ b222 * boundary(sn, tn10, condition);
            }
            if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                e = N;

                const double alpha = condition.alpha();
                const double beta  = condition.beta();
                /*gamma = condition.gamma();*/
                const double gamma = 1.0;

                ax[N]  = beta  * b211;
                bx[N]  = beta  * b213 + alpha * b212;
                cx[N]  = 0.0;

                dx[N]  = beta  * b214 * p10[j][N-1] + alpha  * b215 * p10[j][N] + beta  * b216 * p10[j][N];
                dx[N] += beta  * b217 * p05[j][N-1] + alpha  * b218 * p05[j][N] + beta  * b219 * p05[j][N];

                dx[N] += beta  * b220 * (p10[j+1][N]-2.0*p10[j][N]+p10[j-1][N]);
                dx[N] += beta  * b221 * (p10[j+1][N]-p10[j-1][N]);

                dx[N] += beta  * htht_0250 * f(sn, tn10);
                dx[N] += gamma * b222 * boundary(sn, tn10, condition);
            }

            tomasAlgorithmLeft2Right(ax+s, bx+s, cx+s, dx+s, rx+s, e-s+1);
            for (unsigned int i=s; i<=e; i++) p15[j][i] = rx[i];
        }

        for (n=xmin, sn.i=n, sn.x=n*hx, i=0; n<=xmax; ++n, sn.i=n, sn.x=n*hx, ++i)
        {
            sn.j = ymin; sn.y = ymin*hy;
            value = boundary(sn, tn15, condition);

            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                p15[0][i] = /*(gamma/alpha)**/value;
            }
            if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
                p15[0][i] = (3.5*p15[1][i] - 2.0*p15[2][i] + 0.5*p15[3][i] + hy*/*(gamma/beta)**/value)/(2.0);
                //p15[0][i] = (p15[1][i] + hy*(gamma/beta)*value));
            }
            if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                const double alpha = condition.alpha();
                const double beta  = condition.beta();
                /*gamma = condition.gamma();*/
                const double gamma = 1.0;

                p15[0][i] = (3.5*p15[1][i] - 2.0*p15[2][i] + 0.5*p15[3][i] + hy*(gamma/beta)*value)/(2.0 + (alpha/beta)*hy);
                //p15[0][i] = (p15[1][i] + hy*(gamma/beta)*value)/(1.0+hy*(alpha/beta));
            }

            sn.j = ymax; sn.y = ymax*hy;
            value = boundary(sn, tn15, condition);

            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                p15[M][i] = /*(gamma/alpha)**/value;
            }
            if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
                p15[M][i] = (3.5*p15[M-1][i] - 2.0*p15[M-2][i] + 0.5*p15[M-3][i] + hy*/*(gamma/beta)**/value)/(2.0);
                //p15[M][i] = (p15[M-1][i] + hy*(gamma/beta)*value));
            }
            if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                const double alpha = condition.alpha();
                const double beta  = condition.beta();
                /*gamma = condition.gamma();*/
                const double gamma = 1.0;

                p15[M][i] = (3.5*p15[M-1][i] - 2.0*p15[M-2][i] + 0.5*p15[M-3][i] + hy*(gamma/beta)*value)/(2.0 + (alpha/beta)*hy);
                //p15[M][i] = (p15[M-1][i] + hy*(gamma/beta)*value)/(1.0 + hy*(alpha/beta));
            }
        }

        layerInfo(DoubleMatrix(p15, M+1, N+1), tn15);

        /**************************************************** x direction apprx ***************************************************/

        /**************************************************** y direction apprx ***************************************************/

        for (n=xmin+1, sn.i=n, sn.x=n*hx, i=1; n<=xmax-1; ++n, sn.i=n, sn.x=n*hx, ++i)
        {
            for (m=ymin+1, sn.j=m, sn.y=m*hy, j=1; m<=ymax-1; ++m, sn.j=m, sn.y=m*hy, ++j)
            {
                dy[j] = 0.0;
                dy[j] += k204*p15[j-1][i] + k205*p15[j][i] + k206*p15[j+1][i];
                dy[j] += k207*p10[j-1][i] + k208*p10[j][i] + k209*p10[j+1][i];
                dy[j] += k210*p15[j][i-1] + k211*p15[j][i] + k212*p15[j][i+1];
                //dy[j] += htht_0250 * f(sn, tn15);
                dy[j] += htht_0250 * (0.5*(f(sn, tn10)+f(sn, tn20)));
            }

            sn.j = ymin; sn.y = ymin*hy;
            value = boundary(sn, tn20, condition);

            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                s = 1;
                p20[0][i] = /*(gamma/alpha)**/value;
                dy[1] -= k201 * p20[0][i];
                ay[1] = ay[0] = by[0] = cy[0] = dy[0] = ry[0] = 0.0;
            }
            if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
                s = 0;
                ay[0]  = 0.0;
                by[0]  = /*beta  **/ b311;
                cy[0]  = /*beta  **/ b313;

                dy[0]  = /*beta  **/ b314 * p15[0][i] + /*beta  **/ b316 * p15[1][i];
                dy[0] += /*beta  **/ b317 * p10[0][i] + /*beta  **/ b319 * p10[1][i];

                dy[0] += /*beta  **/ b320 * (p15[0][i+1]-2.0*p15[0][i]+p15[0][i-1]);
                dy[0] += /*beta  **/ b321 * (p15[0][i+1]-p15[0][i-1]);

                dy[0] += /*beta  **/ htht_0250 * f(sn, tn15);
                dy[0] += /*gamma **/ b322 * boundary(sn, tn15, condition);
            }
            if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                s = 0;

                const double alpha = condition.alpha();
                const double beta  = condition.beta();
                /*gamma = condition.gamma();*/
                const double gamma = 1.0;

                ay[0]  = 0.0;
                by[0]  = beta  * b311 + alpha * b312;
                cy[0]  = beta  * b313;

                dy[0]  = beta  * b314 * p15[0][i] + alpha  * b315 * p15[0][i] + beta  * b316 * p15[1][i];
                dy[0] += beta  * b317 * p10[0][i] + alpha  * b318 * p10[0][i] + beta  * b319 * p10[1][i];

                dy[0] += beta  * b320 * (p15[0][i+1]-2.0*p15[0][i]+p15[0][i-1]);
                dy[0] += beta  * b321 * (p15[0][i+1]-p15[0][i-1]);

                dy[0] += beta  * htht_0250 * f(sn, tn15);
                dy[0] += gamma * b322 * boundary(sn, tn15, condition);
            }

            sn.j = ymax; sn.y = ymax*hy;
            value = boundary(sn, tn20, condition);

            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                e = M-1;
                p20[M][i] = /*(gamma/alpha)**/value;
                dy[M-1] -= k203 * p20[M][i];
                cy[M-1] = ay[M] = by[M] = cy[M] = dy[M] = ry[M] = 0.0;
            }
            if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
                e = M;
                ay[M]  = /*beta  **/ b411;
                by[M]  = /*beta  **/ b413;
                cy[M]  = 0.0;

                dy[M]  = /*beta  **/ b414 * p15[M-1][i] + /*beta  **/ b416 * p15[M][i];
                dy[M] += /*beta  **/ b417 * p10[M-1][i] + /*beta  **/ b419 * p10[M][i];

                dy[M] += /*beta  **/ b420 * (p15[M][i+1]-2.0*p15[M][i]+p15[M][i-1]);
                dy[M] += /*beta  **/ b421 * (p15[M][i+1]-p15[M][i-1]);

                dy[M] += /*beta  **/ htht_0250 * f(sn, tn15);
                dy[M] += /*gamma **/ b422 * boundary(sn, tn15, condition);
            }
            if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                e = M;

                const double alpha = condition.alpha();
                const double beta  = condition.beta();
                /*gamma = condition.gamma();*/
                const double gamma = 1.0;

                ay[M]  = beta  * b411;
                by[M]  = beta  * b413 + alpha * b412;
                cy[M]  = 0.0;

                dy[M]  = beta  * b414 * p15[M-1][i] + alpha  * b415 * p15[M][i] + beta  * b416 * p15[M][i];
                dy[M] += beta  * b417 * p10[M-1][i] + alpha  * b418 * p10[M][i] + beta  * b419 * p10[M][i];

                dy[M] += beta  * b420 * (p15[M][i+1]-2.0*p15[M][i]+p15[M][i-1]);
                dy[M] += beta  * b421 * (p15[M][i+1]-p15[M][i-1]);

                dy[M] += beta  * htht_0250 * f(sn, tn15);
                dy[M] += gamma * b422 * boundary(sn, tn15, condition);
            }

            tomasAlgorithmLeft2Right(ay+s, by+s, cy+s, dy+s, ry+s, e-s+1);
            for (unsigned int j=s; j<=e; j++) p20[j][i] = ry[j];
        }

        for (m=ymin, sn.j=m, sn.y=m*hy, j=0; m<=ymax; ++m, sn.j=m, sn.y=m*hy, ++j)
        {
            sn.i = xmin; sn.x = xmin*hx;
            value = boundary(sn, tn20, condition);

            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                p20[j][0] = /*(gamma/alpha)**/value;
            }
            if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
                p20[j][0] = (3.5*p20[j][1] - 2.0*p20[j][2] + 0.5*p20[j][3] + hx*/*(gamma/beta)**/value)/(2.0);
                //p20[j][0] = (p20[j][1] + hx*(gamma/beta)*value);
            }
            if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                const double alpha = condition.alpha();
                const double beta  = condition.beta();
                /*gamma = condition.gamma();*/
                const double gamma = 1.0;

                p20[j][0] = (3.5*p20[j][1] - 2.0*p20[j][2] + 0.5*p20[j][3] + hx*(gamma/beta)*value)/(2.0 + (alpha/beta)*hx);
                //p20[j][0] = (p20[j][1] + hx*(gamma/beta)*value)/(1.0+hx*(alpha/beta));
            }

            sn.i = xmax; sn.x = xmax*hx;
            value = boundary(sn, tn20, condition);

            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                p20[j][N] = /*(gamma/alpha)**/value;
            }
            if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
                p20[j][N] = (3.5*p20[j][N-1] - 2.0*p20[j][N-2] + 0.5*p20[j][N-3] + hx*/*(gamma/beta)**/value)/(2.0);
                //p20[j][N] = (p20[j][N-1] + hx*(gamma/beta)*value);
            }
            if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                const double alpha = condition.alpha();
                const double beta  = condition.beta();
                /*gamma = condition.gamma();*/
                const double gamma = 1.0;

                p20[j][N] = (3.5*p20[j][N-1] - 2.0*p20[j][N-2] + 0.5*p20[j][N-3] + hx*(gamma/beta)*value)/(2.0 + (alpha/beta)*hx);
                //p20[j][N] = (p20[j][N-1] + hx*(gamma/beta)*value)/(1.0 + hx*(alpha/beta));
            }
        }

        //layerInfo(p20, tn20);
        layerInfo(DoubleMatrix(p20, M+1, N+1), tn15);

        /**************************************************** y direction apprx ***************************************************/

        //        for (unsigned int m=0; m<=M; m++)
        //        {
        //            for (unsigned int n=0; n<=N; n++)
        //            {
        //                p00[m][n] = p10[m][n];
        //                p05[m][n] = p15[m][n];
        //                p10[m][n] = p20[m][n];
        //            }
        //        }

        p00 = p10; p05 = p15; p10 = p20;
    }

    for (unsigned int i=0; i<=M; i++) { free(p00[i]); free(p05[i]); free(p10[i]); free(p15[i]); free(p20[i]); }
    free(p00); free(p05); free(p10); free(p15); free(p20);

    //    p00.clear();
    //    p05.clear();
    //    p10.clear();
    //    p15.clear();
    //    p20.clear();

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
