#include "pibvp.h"

void IParabolicIBVP::layerInfo(const DoubleVector &, const TimeNodePDE &) const {}

void IParabolicIBVP::layerInfo(const DoubleMatrix &, const TimeNodePDE &) const {}

//--------------------------------------------------------------------------------------------------------------//xd

void IParabolicFBVP::layerInfo(const DoubleVector &, const TimeNodePDE &) const {}

void IParabolicFBVP::layerInfo(const DoubleMatrix &, const TimeNodePDE &) const {}

//--------------------------------------------------------------------------------------------------------------//

IHeatEquationIBVP::IHeatEquationIBVP(double thermalDiffusivity, double thermalConductivity, double thermalConvection) : IParabolicIBVP(),
    _thermalDiffusivity(thermalDiffusivity), _thermalConductivity(thermalConductivity), _thermalConvection(thermalConvection) {}

IHeatEquationIBVP::IHeatEquationIBVP(const IHeatEquationIBVP &other) : IParabolicIBVP(other),
    _thermalDiffusivity(other._thermalDiffusivity), _thermalConductivity(other._thermalConductivity), _thermalConvection(other._thermalConvection) {}

IHeatEquationIBVP & IHeatEquationIBVP::operator=(const IHeatEquationIBVP &other)
{
    if (this == &other) { return *this; }

    this->_thermalDiffusivity = other._thermalDiffusivity;
    this->_thermalConductivity = other._thermalConductivity;
    this->_thermalConvection = other._thermalConvection;
    return *this;
}

IHeatEquationIBVP::~IHeatEquationIBVP() {}

double IHeatEquationIBVP::thermalDiffusivity() const { return _thermalDiffusivity; }

void IHeatEquationIBVP::setThermalDiffusivity(double thermalDiffusivity) { this->_thermalDiffusivity = thermalDiffusivity; }

double IHeatEquationIBVP::thermalConvection() const { return _thermalConvection; }

void IHeatEquationIBVP::setThermalConvection(double thermalConvection) { this->_thermalConvection = thermalConvection; }

double IHeatEquationIBVP::thermalConductivity() const { return _thermalConductivity; }

void IHeatEquationIBVP::setThermalConductivity(double thermalConductivity) { this->_thermalConductivity = thermalConductivity; }

void IHeatEquationIBVP::explicit_calculate_D1V1() const {}

void IHeatEquationIBVP::implicit_calculate_D1V1() const
{
    const unsigned int N = static_cast<unsigned int>(spaceDimensionX().size()) - 1;
    const unsigned int L = static_cast<unsigned int>(timeDimension().size()) - 1;

    const int xmin = spaceDimensionX().min();
    const int xmax = spaceDimensionX().max();
    //const int tmin = timeDimension().min();
    //const int tmax = timeDimension().max();

    const double hx = spaceDimensionX().step();
    const double ht = timeDimension().step();

    const double a = thermalDiffusivity();
    const double b = thermalConductivity();
    const double c = thermalConvection();
    const double w1 = weight();
    const double w2 = 1.0 - w1;

    // equation parameters
    const double k11 = -((a*ht)/(hx*hx))*w1 + ((b*ht)/(2.0*hx))*w1;
    const double k12 = +1.0 + ((2.0*a*ht)/(hx*hx))*w1 - c*ht*w1;
    const double k13 = -((a*ht)/(hx*hx))*w1 - ((b*ht)/(2.0*hx))*w1;
    const double k21 = +((a*ht)/(hx*hx))*w2 - ((b*ht)/(2.0*hx))*w2;
    const double k22 = +1.0 - ((2.0*a*ht)/(hx*hx))*w2 + c*ht*w2;
    const double k23 = +((a*ht)/(hx*hx))*w2 + ((b*ht)/(2.0*hx))*w2;
    
    // left border condition parameters
    const double b11 = +1.0 + ((2.0*a*ht)/(hx*hx))*w1 - c*ht*w1;
    const double b12 = -((2.0*a*ht)/hx)*w1 + b*ht*w1;
    const double b13 = -((2.0*a*ht)/(hx*hx))*w1;
    const double b14 = +1.0 - ((2.0*a*ht)/(hx*hx))*w2 + c*ht*w2;
    const double b15 = +((2.0*a*ht)/hx)*w2 - b*ht*w2;
    const double b16 = +((2.0*a*ht)/(hx*hx))*w2;
    const double b17 = -((2.0*a*ht)/hx)*w1 + b*ht*w1;
    const double b18 = -((2.0*a*ht)/hx)*w2 + b*ht*w2;

    // right border condition parameters
    const double b21 = -((2.0*a*ht)/(hx*hx))*w1;
    const double b22 = +((2.0*a*ht)/hx)*w1 + b*ht*w1;
    const double b23 = +1.0 + ((2.0*a*ht)/(hx*hx))*w1 - c*ht*w1;
    const double b24 = +((2.0*a*ht)/(hx*hx))*w2;
    const double b25 = -((2.0*a*ht)/hx)*w2 - b*ht*w2;
    const double b26 = +1.0 - ((2.0*a*ht)/(hx*hx))*w2 + c*ht*w2;
    const double b27 = +((2.0*a*ht)/hx)*w1 + b*ht*w1;
    const double b28 = +((2.0*a*ht)/hx)*w2 + b*ht*w2;

    // common parameters
    const double ht_w1 = +ht*w1;
    const double ht_w2 = +ht*w2;

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

    /***********************************************************************************************/
    /***********************************************************************************************/

    TimeNodePDE tn0; tn0.i = 0; tn0.t = tn0.i*ht;
    SpaceNodePDE sn;
    unsigned int i = 0;
    for (int n=xmin; n<=xmax; n++, i++)
    {
        sn.i = n; sn.x = n*hx;
        u00[i] = initial(sn, InitialCondition::InitialValue);
    }
    layerInfo(u00, tn0);

    /***********************************************************************************************/
    /***********************************************************************************************/

    for (unsigned int ln=1; ln<=L; ln++)
    {
        TimeNodePDE tn00; tn00.i = ln-1; tn00.t = tn00.i*ht;
        TimeNodePDE tn10; tn10.i = ln;   tn10.t = tn10.i*ht;

        unsigned int i = 1;
        for (int n=xmin+1; n<=xmax-1; n++, i++)
        {
            sn.i = n; sn.x = n*hx;
            dx[i] = k21*u00[i-1] + k22*u00[i] + k23*u00[i+1] + (ht_w1*f(sn, tn10)+ht_w2*f(sn, tn00));
        }

        unsigned int s=0, e=N;
        BoundaryConditionPDE condition; double value, alpha, beta, gamma;
        BoundaryConditionPDE condition0; double value0;

        sn.i = xmin; sn.x = xmin*hx;
        value = boundary(sn, tn10, condition);
        value0 = boundary(sn, tn00, condition);
        alpha = condition.alpha();
        beta  = condition.beta();
        gamma = condition.gamma();
        if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
        {
            s = 1;
            u10[0] = (gamma/alpha)*value;
            dx[1] -= k11 * u10[0];
            ax[1] = ax[0] = bx[0] = cx[0] = dx[0] = rx[0] = 0.0;
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
        {
            s = 0;
            ax[s]  = 0.0;
            bx[s]  = beta*b11;
            cx[s]  = beta*b13;
            dx[s]  = beta*b14*u00[s] + beta*b16*u00[s+1]
                    + beta*(ht_w1*f(sn, tn10)+ht_w2*f(sn, tn00)) + gamma*(b17*value+b18*value0);
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Robin)
        {
            s = 0;
            ax[s]  = 0.0;
            bx[s]  = beta*b11 + alpha*b12;
            cx[s]  = beta*b13;
            dx[s]  = beta*b14*u00[s] + alpha*b15*u00[s] + beta*b16*u00[s+1]
                    + beta*(ht_w1*f(sn, tn10)+ht_w2*f(sn, tn00)) + gamma*(b17*value+b18*value0);
        }

        sn.i = xmax; sn.x = xmax*hx;
        value = boundary(sn, tn10, condition);
        value0 = boundary(sn, tn00, condition);
        alpha = condition.alpha();
        beta  = condition.beta();
        gamma = condition.gamma();
        if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
        {
            e = N-1;
            u10[N] = (gamma/alpha)*value;
            dx[N-1] -= k13 * u10[N];
            cx[N-1] = ax[N] = bx[N] = cx[N] = dx[N] = rx[N] = 0.0;
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
        {
            e = N;
            ax[e]  = beta*b21;
            bx[e]  = beta*b23;
            cx[e]  = 0.0;
            dx[e]  = beta*b24*u00[e-1] + beta*b26*u00[e]
                    + beta*(ht_w1*f(sn, tn10)+ht_w2*f(sn, tn00)) + gamma*(b27*value+b28*value0);
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Robin)
        {
            e = N;
            ax[e]  = beta*b21;
            bx[e]  = alpha*b22 + beta*b23;
            cx[e]  = 0.0;
            dx[e]  = beta*b24*u00[e-1] + alpha*b25*u00[e] + beta*b26*u00[e]
                    + beta*(ht_w1*f(sn, tn10)+ht_w2*f(sn, tn00)) + gamma*(b27*value+b28*value0);
        }

        tomasAlgorithm(ax+s, bx+s, cx+s, dx+s, rx+s, e-s+1);
        for (unsigned int n=s; n<=e; n++) u10[n] = rx[n];
        layerInfo(u10, tn10);

        u00 = u10;
    }

    u00.clear();
    u10.clear();

    free(rx);
    free(dx);
    free(cx);
    free(bx);
    free(ax);
}

void IHeatEquationIBVP::explicit_calculate_D2V1() const
{
    const unsigned int N = static_cast<unsigned int>( spaceDimensionX().size() );
    const unsigned int M = static_cast<unsigned int>( spaceDimensionY().size() );
    const unsigned int L = static_cast<unsigned int>( timeDimension().size() );

    const double hx = spaceDimensionX().step();
    const double hy = spaceDimensionY().step();
    const double ht = timeDimension().step();

    const double td = thermalDiffusivity();
    const double tc = thermalConductivity();

    if (ht > 0.5/(1.0/(hx*hx)+1.0/(hy*hy))) throw std::runtime_error("Differential scheme not steady!");

    const double p_td_ht__hxhx = +((td*ht)/(hx*hx));
    const double p_td_ht__hyhy = +((td*ht)/(hy*hy));
    const double ht_tc = ht*tc;

    DoubleMatrix u00(M+1, N+1);
    DoubleMatrix u10(M+1, N+1);

    SpaceNodePDE sn;
    TimeNodePDE tn;

    /***********************************************************************************************/
    /***********************************************************************************************/

    for (unsigned int m=0; m<=M; m++)
    {
        sn.j = static_cast<int>(m); sn.y = sn.j*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = static_cast<int>(n); sn.x = sn.i*hx;
            u00[m][n] = initial(sn, InitialCondition::InitialValue);
        }
    }
    layerInfo(u00, tn);

    /***********************************************************************************************/

    for (unsigned int ln=1; ln<=L; ln++)
    {
        tn.i = ln; tn.t = tn.i*ht;

        /**************************************************** border conditions ***************************************************/

        SpaceNodePDE sn0, sn1;
        BoundaryConditionPDE condition;

        sn0.i = static_cast<int>(0); sn0.x = 0*hx;
        sn1.i = static_cast<int>(N); sn1.x = N*hx;
        for (unsigned int m=0; m<=M; m++)
        {
            sn0.j = sn1.j = static_cast<int>(m); sn0.y = sn1.y = m*hy;
            u10[m][0] = boundary(sn0, tn, condition);
            u10[m][N] = boundary(sn1, tn, condition);
        }

        sn0.j = static_cast<int>(0); sn0.y = 0*hy;
        sn1.j = static_cast<int>(M); sn1.y = M*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn0.i = sn1.i = static_cast<int>(n); sn0.x = sn1.x = n*hx;
            u10[0][n] = boundary(sn0, tn, condition);
            u10[M][n] = boundary(sn1, tn, condition);
        }

        /**************************************************** border conditions ***************************************************/

        /**************************************************************************************************************************/
        for (unsigned int m=1; m<=M-1; m++)
        {
            sn.j = static_cast<int>(m); sn.y = m*hy;
            for (unsigned int n=1; n<=N-1; n++)
            {
                sn.i = static_cast<int>(n); sn.x = n*hx;
                u10[m][n] = u00[m][n]
                        + p_td_ht__hxhx*(u00[m][n-1] - 2.0*u00[m][n] + u00[m][n+1])
                        + p_td_ht__hyhy*(u00[m-1][n] - 2.0*u00[m][n] + u00[m+1][n])
                        - ht_tc*u00[m][n] + ht*f(sn, tn);
            }
        }
        layerInfo(u10, tn);
        /**************************************************************************************************************************/

        u00 = u10;
    }

    u00.clear();
    u10.clear();
}

void IHeatEquationIBVP::implicit_calculate_D2V1() const
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

    const double a1 = thermalDiffusivity();
    const double a2 = thermalDiffusivity();
    const double b1 = thermalConductivity();
    const double b2 = thermalConductivity();
    const double c  = thermalConvection();

    //if (w1 >= 0.5-(0.25/ht)/(1.0/(hx*hx)+1.0/(hy*hy))) throw std::runtime_error("Differential scheme is conditionally steady.");

    // common parameters
    const double ht05 = 0.5*ht;

    // equation parameters
    const double k11 = -(a1*ht)/(2.0*hx*hx) + (b1*ht)/(4.0*hx);
    const double k12 = +1.0 + (a1*ht)/(hx*hx);
    const double k13 = -(a1*ht)/(2.0*hx*hx) - (b1*ht)/(4.0*hx);
    const double k14 = +(a2*ht)/(2.0*hy*hy) - (b2*ht)/(4.0*hy);
    const double k15 = +1.0 - (a2*ht)/(hy*hy) + 0.5*c*ht;
    const double k16 = +(a2*ht)/(2.0*hy*hy) + (b2*ht)/(4.0*hy);

    const double k21 = -(a2*ht)/(2.0*hy*hy) + (b2*ht)/(4.0*hy);
    const double k22 = +1.0 + (a2*ht)/(hy*hy) - 0.5*c*ht;
    const double k23 = -(a2*ht)/(2.0*hy*hy) - (b2*ht)/(4.0*hy);
    const double k24 = +(a1*ht)/(2.0*hx*hx) - (b1*ht)/(4.0*hx);
    const double k25 = +1.0 - (a1*ht)/(hx*hx);
    const double k26 = +(a1*ht)/(2.0*hx*hx) + (b1*ht)/(4.0*hx);

    // left border condition parameters
    const double b111 = +1.0 + ((a1*ht)/(hx*hx));
    const double b112 = +(a1*ht)/hx - 0.5*b1*ht;
    const double b113 = -(a1*ht)/(hx*hx);
    const double b114 = +1.0 + 0.5*c*ht;
    const double b115 = +0.5*a2*ht;
    const double b116 = +0.5*b2*ht;
    const double b117 = +(a1*ht)/hx - 0.5*b1*ht;

    // right border condition parameters
    const double b121 = -(a1*ht)/(hx*hx);
    const double b122 = +(a1*ht)/hx + 0.5*b1*ht;
    const double b123 = +1.0 + ((a1*ht)/(hx*hx));
    const double b124 = +1.0 + 0.5*c*ht;
    const double b125 = +0.5*a2*ht;
    const double b126 = +0.5*b2*ht;
    const double b127 = +(a1*ht)/hx + 0.5*b1*ht;

    // bottom border condition parameters
    const double b211 = +1.0 + ((a2*ht)/(hy*hy)) - 0.5*c*ht;
    const double b212 = +(a2*ht)/hy - 0.5*b2*ht;
    const double b213 = -(a2*ht)/(hy*hy);
    const double b214 = +1.0;
    const double b215 = +0.5*a1*ht;
    const double b216 = +0.5*b1*ht;
    const double b217 = +(a2*ht)/hy - 0.5*b2*ht;

    // top border condition parameters
    const double b221 = -(a2*ht)/(hy*hy);
    const double b222 = +(a2*ht)/hy + 0.5*b2*ht;
    const double b223 = +1.0 + ((a2*ht)/(hy*hy)) - 0.5*c*ht;
    const double b224 = +1.0;
    const double b225 = +0.5*a1*ht;
    const double b226 = +0.5*b1*ht;
    const double b227 = +(a2*ht)/hy + 0.5*b2*ht;

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

    double *ay = static_cast<double*>(malloc(sizeof(double)*(M+1)));
    double *by = static_cast<double*>(malloc(sizeof(double)*(M+1)));
    double *cy = static_cast<double*>(malloc(sizeof(double)*(M+1)));
    double *dy = static_cast<double*>(malloc(sizeof(double)*(M+1)));
    double *ry = static_cast<double*>(malloc(sizeof(double)*(M+1)));

    for (unsigned int m=0; m<=M; m++)
    {
        ay[m] = k21;
        by[m] = k22;
        cy[m] = k23;
    }
    ay[0] = 0.0; cy[M] = 0.0;

    DoubleMatrix u00(M+1, N+1);
    DoubleMatrix u05(M+1, N+1);
    DoubleMatrix u10(M+1, N+1);

    /***********************************************************************************************/
    /***********************************************************************************************/

    TimeNodePDE tn00, tn05, tn10;
    SpaceNodePDE sn;
    unsigned int i=0, j=0, s=0, e=N;
    int n = 0, m = 0;
    BoundaryConditionPDE condition; double value, alpha, beta, gamma;

    for (m=ymin, sn.j=m, sn.y=m*hy, j=0; m<=ymax; ++m, sn.j=m, sn.y=m*hy, ++j)
    {
        for (n=xmin, sn.i=n, sn.x=n*hx, i=0; n<=xmax; ++n, sn.i=n, sn.x=n*hx, ++i)
        {
            u00[j][i] = initial(sn, InitialCondition::InitialValue);
        }
    }

    tn00.i = 0; tn00.t = 0.5*tn00.i*ht;
    layerInfo(u00, tn00);

    /***********************************************************************************************/

    for (unsigned int ln=1; ln<=L; ln++)
    {
        tn00.i = 2*ln-2; tn00.t = 0.5*tn00.i*ht;
        tn05.i = 2*ln-1; tn05.t = 0.5*tn05.i*ht;
        tn10.i = 2*ln-0; tn10.t = 0.5*tn10.i*ht;

        /**************************************************** x direction apprx ***************************************************/

        for (m=ymin+1, sn.j=m, sn.y=m*hy, j=1; m<=ymax-1; ++m, sn.j=m, sn.y=m*hy, ++j)
        {
            for (n=xmin+1, sn.i=n, sn.x=n*hx, i=1, dx[i]=0.0; n<=xmax-1; ++n, sn.i=n, sn.x=n*hx, ++i, dx[i]=0.0)
            {
                dx[i] += k14*u00[j-1][i] + k15*u00[j][i] + k16*u00[j+1][i] + ht05*f(sn, tn00);
            }

            sn.i = xmin; sn.x = xmin*hx;
            value = boundary(sn, tn05, condition);
            alpha = condition.alpha();
            beta  = condition.beta();
            gamma = condition.gamma();

            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                s = 1;
                u05[j][0] = (gamma/alpha)*value;
                dx[1] -= k11 * u05[j][0];
                ax[1] = ax[0] = bx[0] = cx[0] = dx[0] = rx[0] = 0.0;
            }
            if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
                s = 0;
                ax[s]  = 0.0;
                bx[s]  = beta  * b111;
                cx[s]  = beta  * b113;
                dx[s]  = beta  * b114 * u00[j][s];
                dx[s] += beta  * b115 * ((u00[j+1][s]-2.0*u00[j][s]+u00[j-1][s])/(hy*hy));
                dx[s] += beta  * b116 * ((u00[j+1][s]-u00[j-1][s])/(2.0*hy));
                dx[s] += beta  * ht05 * f(sn, tn00);
                dx[s] += gamma * b117 * value;
            }
            if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                s = 0;
                ax[s]  = 0.0;
                bx[s]  = beta  * b111 + alpha * b112;
                cx[s]  = beta  * b113;
                dx[s]  = beta  * b114 * u00[j][s];
                dx[s] += beta  * b115 * ((u00[j+1][s]-2.0*u00[j][s]+u00[j-1][s])/(hy*hy));
                dx[s] += beta  * b116 * ((u00[j+1][s]-u00[j-1][s])/(2.0*hy));
                dx[s] += beta  * ht05 * f(sn, tn00);
                dx[s] += gamma * b117 * value;
            }

            sn.i = xmax; sn.x = xmax*hx;
            value = boundary(sn, tn05, condition);
            alpha = condition.alpha();
            beta  = condition.beta();
            gamma = condition.gamma();

            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                e = N-1;
                u05[j][N] = (gamma/alpha)*value;
                dx[N-1] -= k13 * u05[j][N];
                cx[N-1] = ax[N] = bx[N] = cx[N] = dx[N] = rx[N] = 0.0;
            }
            if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
                e = N;
                ax[e]  = beta  * b121;
                bx[e]  = beta  * b123 + alpha * b122;
                cx[e]  = 0.0;
                dx[e]  = beta  * b124 * u00[j][e];
                dx[e] += beta  * b125 * ((u00[j+1][e]-2.0*u00[j][e]+u00[j-1][e])/(hy*hy));
                dx[e] += beta  * b126 * ((u00[j+1][e]-u00[j-1][e])/(2.0*hy));
                dx[e] += beta  * ht05 * f(sn, tn00);
                dx[e] += gamma * b127 * value;
            }
            if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                e = N;
                ax[e]  = beta  * b121;
                bx[e]  = beta  * b123 + alpha * b122;
                cx[e]  = 0.0;
                dx[e]  = beta  * b124 * u00[j][e];
                dx[e] += beta  * b125 * ((u00[j+1][e]-2.0*u00[j][e]+u00[j-1][e])/(hy*hy));
                dx[e] += beta  * b126 * ((u00[j+1][e]-u00[j-1][e])/(2.0*hy));
                dx[e] += beta  * ht05 * f(sn, tn00);
                dx[e] += gamma * b127 * value;
            }

            tomasAlgorithm(ax+s, bx+s, cx+s, dx+s, rx+s, e-s+1);
            for (unsigned int i=s; i<=e; i++) u05[j][i] = rx[i];
        }

        for (n=xmin, sn.i=n, sn.x=n*hx, i=0, dx[i]=0.0; n<=xmax; ++n, sn.i=n, sn.x=n*hx, ++i, dx[i]=0.0)
        {
            sn.j = ymin; sn.y = ymin*hy;
            value = boundary(sn, tn05, condition);
            alpha = condition.alpha();
            beta  = condition.beta();
            gamma = condition.gamma();

            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                u05[0][i] = (gamma/alpha)*value;
            }
            if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
                u05[0][i] = (3.5*u05[1][i] - 2.0*u05[2][i] + 0.5*u05[3][i] + hy*(gamma/beta)*value)/(2.0);
                //u05[0][i] = (u05[1][i] + hy*(gamma/beta)*value));
            }
            if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                u05[0][i] = (3.5*u05[1][i] - 2.0*u05[2][i] + 0.5*u05[3][i] + hy*(gamma/beta)*value)/(2.0 + (alpha/beta)*hy);
                //u05[0][i] = (u05[1][i] + hy*(gamma/beta)*value)/(1.0+hy*(alpha/beta));
            }

            sn.j = ymax; sn.y = ymax*hy;
            value = boundary(sn, tn05, condition);
            alpha = condition.alpha();
            beta  = condition.beta();
            gamma = condition.gamma();

            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                u05[M][i] = (gamma/alpha)*value;
            }
            if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
                u05[M][i] = (3.5*u05[M-1][i] - 2.0*u05[M-2][i] + 0.5*u05[M-3][i] + hy*(gamma/beta)*value)/(2.0);
                //u05[M][i] = (u05[M-1][i] + hy*(gamma/beta)*value));
            }
            if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                u05[M][i] = (3.5*u05[M-1][i] - 2.0*u05[M-2][i] + 0.5*u05[M-3][i] + hy*(gamma/beta)*value)/(2.0 + (alpha/beta)*hy);
                //u05[M][i] = (u05[M-1][i] + hy*(gamma/beta)*value)/(1.0 + hy*(alpha/beta));
            }
        }

        layerInfo(u05, tn05);

        /**************************************************** x direction apprx ***************************************************/

        /**************************************************** y direction apprx ***************************************************/

        for (n=xmin+1, sn.i=n, sn.x=n*hx, i=1; n<=xmax-1; ++n, sn.i=n, sn.x=n*hx, ++i)
        {
            for (m=ymin+1, sn.j=m, sn.y=m*hy, j=1, dy[j]=0.0; m<=ymax-1; ++m, sn.j=m, sn.y=m*hy, ++j, dy[j]=0.0)
            {
                dy[j] += k24*u05[j][i-1] + k25*u05[j][i] + k26*u05[j][i+1] + ht05*f(sn, tn10);
            }

            sn.j = ymin; sn.y = ymin*hy;
            value = boundary(sn, tn10, condition);
            alpha = condition.alpha();
            beta  = condition.beta();
            gamma = condition.gamma();

            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                s = 1;
                u10[0][i] = (gamma/alpha)*value;
                dy[1] -= k21 * u10[0][i];
                ay[1] = ay[0] = by[0] = cy[0] = dy[0] = ry[0] = 0.0;
            }
            if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
                s = 0;
                ay[s]  = 0.0;
                by[s]  = beta  * b211;
                cy[s]  = beta  * b213;
                dy[s]  = beta  * b214 * u05[s][i];
                dy[s] += beta  * b215 * ((u05[s][i+1]-2.0*u05[s][i]+u05[s][i-1])/(hx*hx));
                dy[s] += beta  * b216 * ((u05[s][i+1]-u05[s][i-1])/(2.0*hx));
                dy[s] += beta  * ht05 * f(sn, tn10);
                dy[s] += gamma * b217 * value;
            }
            if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                s = 0;
                ay[s]  = 0.0;
                by[s]  = beta  * b211 + alpha * b212;
                cy[s]  = beta  * b213;
                dy[s]  = beta  * b214 * u05[s][i];
                dy[s] += beta  * b215 * ((u05[s][i+1]-2.0*u05[s][i]+u05[s][i-1])/(hx*hx));
                dy[s] += beta  * b216 * ((u05[s][i+1]-u05[s][i-1])/(2.0*hx));
                dy[s] += beta  * ht05 * f(sn, tn10);
                dy[s] += gamma * b217 * value;
            }

            sn.j = ymax; sn.y = ymax*hy;
            value = boundary(sn, tn10, condition);
            alpha = condition.alpha();
            beta  = condition.beta();
            gamma = condition.gamma();

            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                e = M-1;
                u10[M][i] = (gamma/alpha)*value;
                dy[M-1] -= k23 * u10[M][i];
                cy[M-1] = ay[M] = by[M] = cy[M] = dy[M] = ry[M] = 0.0;
            }
            if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
                e = M;
                ay[e]  = beta  * b221;
                by[e]  = beta  * b223;
                cy[e]  = 0.0;
                dy[e]  = beta  * b224 * u05[e][i];
                dy[e] += beta  * b225 * ((u05[e][i+1]-2.0*u05[e][i]+u05[e][i-1])/(hx*hx));
                dy[e] += beta  * b226 * ((u05[e][i+1]-u05[e][i-1])/(2.0*hx));
                dy[e] += beta  * ht05 * f(sn, tn10);
                dy[e] += gamma * b227 * value;
            }
            if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                e = M;
                ay[e]  = beta  * b221;
                by[e]  = beta  * b223 + alpha * b222;
                cy[e]  = 0.0;
                dy[e]  = beta  * b224 * u05[e][i];
                dy[e] += beta  * b225 * ((u05[e][i+1]-2.0*u05[e][i]+u05[e][i-1])/(hx*hx));
                dy[e] += beta  * b226 * ((u05[e][i+1]-u05[e][i-1])/(2.0*hx));
                dy[e] += beta  * ht05 * f(sn, tn10);
                dy[e] += gamma * b227 * value;
            }

            tomasAlgorithm(ay+s, by+s, cy+s, dy+s, ry+s, e-s+1);
            for (unsigned int j=s; j<=e; j++) u10[j][i] = ry[j];
        }

        for (m=ymin, sn.j=m, sn.y=m*hy, j=0, dy[j]=0.0; m<=ymax; ++m, sn.j=m, sn.y=m*hy, ++j, dy[j]=0.0)
        {
            sn.i = xmin; sn.x = xmin*hx;
            value = boundary(sn, tn10, condition);
            alpha = condition.alpha();
            beta  = condition.beta();
            gamma = condition.gamma();

            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                u10[j][0] = (gamma/alpha)*value;
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
            value = boundary(sn, tn10, condition);
            alpha = condition.alpha();
            beta  = condition.beta();
            gamma = condition.gamma();

            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                u10[j][N] = (gamma/alpha)*value;
            }
            if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
                u10[j][N] = (3.5*u10[j][N-1] - 2.0*u10[j][N-2] + 0.5*u10[j][N-3] + hx*(gamma/beta)*value)/(2.0);
                //u10[j][N] = (u10[j][N-1] + hx*(gamma/beta)*value));
            }
            if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                u10[j][N] = (3.5*u10[j][N-1] - 2.0*u10[j][N-2] + 0.5*u10[j][N-3] + hx*(gamma/beta)*value)/(2.0 + (alpha/beta)*hx);
                //u10[j][N] = (u10[j][N-1] + hx*(gamma/beta)*value)/(1.0 + hx*(alpha/beta));
            }
        }

        layerInfo(u10, tn10);

        /**************************************************** y direction apprx ***************************************************/

        for (unsigned int j=0; j<=M; j++)
        {
            for (unsigned int i=0; i<=N; i++)
            {
                u00[j][i] = u10[j][i];
            }
        }
    }

    u00.clear();
    u05.clear();
    u10.clear();

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

double IHeatEquationIBVP::weight() const { return 0.5; }

//--------------------------------------------------------------------------------------------------------------//

IHeatEquationFBVP::IHeatEquationFBVP(double thermalDiffusivity, double thermalConductivity, double thermalConvection) : IParabolicFBVP(),
    _thermalDiffusivity(thermalDiffusivity), _thermalConductivity(thermalConductivity), _thermalConvection(thermalConvection) {}

IHeatEquationFBVP::IHeatEquationFBVP(const IHeatEquationFBVP &other) : IParabolicFBVP(other),
    _thermalDiffusivity(other._thermalDiffusivity), _thermalConductivity(other._thermalConductivity), _thermalConvection(other._thermalConvection) {}

IHeatEquationFBVP & IHeatEquationFBVP::operator=(const IHeatEquationFBVP &other)
{
    if (this == &other) { return *this; }

    this->_thermalDiffusivity = other._thermalDiffusivity;
    this->_thermalConductivity = other._thermalConductivity;
    this->_thermalConvection = other._thermalConvection;
    return *this;
}

IHeatEquationFBVP::~IHeatEquationFBVP() {}

double IHeatEquationFBVP::thermalDiffusivity() const { return _thermalDiffusivity; }

void IHeatEquationFBVP::setThermalDiffusivity(double thermalDiffusivity) { this->_thermalDiffusivity = thermalDiffusivity; }

double IHeatEquationFBVP::thermalConductivity() const { return _thermalConductivity; }

void IHeatEquationFBVP::setThermalConductivity(double convection) { this->_thermalConductivity = convection; }

double IHeatEquationFBVP::thermalConvection() const { return _thermalConvection; }

void IHeatEquationFBVP::setThermalConvection(double thermalConvection) { this->_thermalConvection = thermalConvection; }

void IHeatEquationFBVP::explicit_calculate_D1V1() const {}

void IHeatEquationFBVP::implicit_calculate_D1V1() const
{
    const unsigned int N = static_cast<unsigned int>(spaceDimensionX().size()) - 1;
    const unsigned int M = static_cast<unsigned int>(timeDimension().size()) - 1;

    const int xmin = spaceDimensionX().min();
    const int xmax = spaceDimensionX().max();
    //const int tmin = timeDimension().min();
    //const int tmax = timeDimension().max();

    const double hx = spaceDimensionX().step();
    const double ht = timeDimension().step();

    const double a = thermalDiffusivity();
    const double b = thermalConductivity();
    const double c = thermalConvection();
    const double w1 = weight();
    const double w2 = 1.0 - w1;

    // equation parameters
    const double k11 = -((a*ht)/(hx*hx))*w1 + ((b*ht)/(2.0*hx))*w1;
    const double k12 = -1.0 + ((2.0*a*ht)/(hx*hx))*w1 - c*ht*w1;
    const double k13 = -((a*ht)/(hx*hx))*w1 - ((b*ht)/(2.0*hx))*w1;
    const double k21 = +((a*ht)/(hx*hx))*w2 - ((b*ht)/(2.0*hx))*w2;
    const double k22 = -1.0 - ((2.0*a*ht)/(hx*hx))*w2 + c*ht*w2;
    const double k23 = +((a*ht)/(hx*hx))*w2 + ((b*ht)/(2.0*hx))*w2;

    // left border condition parameters
    const double b11 = -1.0 + ((2.0*a*ht)/(hx*hx))*w1 - c*ht*w1;
    const double b12 = -((2.0*a*ht)/hx)*w1 + b*ht*w1;
    const double b13 = -((2.0*a*ht)/(hx*hx))*w1;
    const double b14 = -1.0 - ((2.0*a*ht)/(hx*hx))*w2 + c*ht*w2;
    const double b15 = +((2.0*a*ht)/hx)*w2 - b*ht*w2;
    const double b16 = +((2.0*a*ht)/(hx*hx))*w2;
    const double b17 = -((2.0*a*ht)/hx)*w1 + b*ht*w1;
    const double b18 = -((2.0*a*ht)/hx)*w2 + b*ht*w2;

    // right border condition parameters
    const double b21 = -((2.0*a*ht)/(hx*hx))*w1;
    const double b22 = +((2.0*a*ht)/hx)*w1 + b*ht*w1;
    const double b23 = -1.0 + ((2.0*a*ht)/(hx*hx))*w1 - c*ht*w1;
    const double b24 = +((2.0*a*ht)/(hx*hx))*w2;
    const double b25 = -((2.0*a*ht)/hx)*w2 - b*ht*w2;
    const double b26 = -1.0 - ((2.0*a*ht)/(hx*hx))*w2 + c*ht*w2;
    const double b27 = +((2.0*a*ht)/hx)*w1 + b*ht*w1;
    const double b28 = +((2.0*a*ht)/hx)*w2 + b*ht*w2;

    // common parameters
    const double ht_w1 = +ht*w1;
    const double ht_w2 = +ht*w2;

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

    /***********************************************************************************************/
    /***********************************************************************************************/

    TimeNodePDE tn0; tn0.i = M; tn0.t = tn0.i*ht;
    SpaceNodePDE sn;
    unsigned int i = 0;
    for (int n=xmin; n<=xmax; n++, i++)
    {
        sn.i = n; sn.x = sn.i*hx;
        u00[i] = final(sn, FinalCondition::FinalValue);
    }
    layerInfo(u00, tn0);

    /***********************************************************************************************/
    /***********************************************************************************************/

    for (unsigned int ln=M-1, size_ln = static_cast<unsigned int>(0)-1; ln != size_ln; ln--)
    {
        TimeNodePDE tn00; tn00.i = ln+1; tn00.t = tn00.i*ht;
        TimeNodePDE tn10; tn10.i = ln;   tn10.t = tn10.i*ht;

        unsigned int i = 1;
        for (int n=xmin+1; n<=xmax-1; n++, i++)
        {
            sn.i = n; sn.x = n*hx;
            dx[i] = k21*u00[i-1] + k22*u00[i] + k23*u00[i+1] + (ht_w1*f(sn, tn10)+ht_w2*f(sn, tn00));
        }

        unsigned int s=0, e=N;
        BoundaryConditionPDE condition; double alpha, beta, gamma, value;
        BoundaryConditionPDE condition0; double value0;

        sn.i = xmin; sn.x = xmin*hx;
        value = boundary(sn, tn10, condition);
        value0 = boundary(sn, tn00, condition);
        alpha = condition.alpha();
        beta  = condition.beta();
        gamma = condition.gamma();
        if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
        {
            s = 1;
            u10[0] = (gamma/alpha)*value;
            dx[1] -= k11 * u10[0];
            ax[1] = ax[0] = bx[0] = cx[0] = dx[0] = rx[0] = 0.0;
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
        {
            s = 0;
            ax[s]  = 0.0;
            bx[s]  = beta*b11;
            cx[s]  = beta*b13;
            dx[s]  = beta*b14*u00[s] + beta*b16*u00[s+1]
                    + beta*(ht_w1*f(sn, tn10)+ht_w2*f(sn, tn00)) + gamma*(b17*value+b18*value0);
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Robin)
        {
            s = 0;
            ax[s]  = 0.0;
            bx[s]  = beta*b11 + alpha*b12;
            cx[s]  = beta*b13;
            dx[s]  = beta*b14*u00[s] + alpha*b15*u00[s] + beta*b16*u00[s+1]
                    + beta*(ht_w1*f(sn, tn10)+ht_w2*f(sn, tn00)) + gamma*(b17*value+b18*value0);
        }

        sn.i = xmax; sn.x = xmax*hx;
        value = boundary(sn, tn10, condition);
        value0 = boundary(sn, tn00, condition);
        alpha = condition.alpha();
        beta  = condition.beta();
        gamma = condition.gamma();
        if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
        {
            e = N-1;
            u10[N] = (gamma/alpha)*value;
            dx[N-1] -= k13 * u10[N];
            cx[N-1] = ax[N] = bx[N] = cx[N] = dx[N] = rx[N] = 0.0;
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
        {
            e = N;
            ax[e]  = beta*b21;
            bx[e]  = beta*b23;
            cx[e]  = 0.0;
            dx[e]  = beta*b24*u00[e-1] + beta*b26*u00[e]
                    + beta*(ht_w1*f(sn, tn10)+ht_w2*f(sn, tn00)) + gamma*(b27*value+b28*value0);
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Robin)
        {
            e = N;
            ax[e]  = beta*b21;
            bx[e]  = alpha*b22 + beta*b23;
            cx[e]  = 0.0;
            dx[e]  = beta*b24*u00[e-1] + alpha*b25*u00[e] + beta*b26*u00[e]
                    + beta*(ht_w1*f(sn, tn10)+ht_w2*f(sn, tn00)) + gamma*(b27*value+b28*value0);
        }

        tomasAlgorithm(ax+s, bx+s, cx+s, dx+s, rx+s, e-s+1);
        for (unsigned int n=s; n<=e; n++) u10[n] = rx[n];
        layerInfo(u10, tn10);

        u00 = u10;
    }

    u00.clear();
    u10.clear();

    free(rx);
    free(dx);
    free(cx);
    free(bx);
    free(ax);
}

void IHeatEquationFBVP::implicit_calculate_D2V1() const
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

    const double a1 = thermalDiffusivity();
    const double a2 = thermalDiffusivity();
    const double b1 = thermalConductivity();
    const double b2 = thermalConductivity();
    const double c  = thermalConvection();

    //if (w1 >= 0.5-(0.25/ht)/(1.0/(hx*hx)+1.0/(hy*hy))) throw std::runtime_error("Differential scheme is conditionally steady.");

    // common parameters
    const double ht05 = 0.5*ht;

    // equation parameters
    const double k11 = -(a1*ht)/(2.0*hx*hx) + (b1*ht)/(4.0*hx);
    const double k12 = -1.0 + (a1*ht)/(hx*hx);
    const double k13 = -(a1*ht)/(2.0*hx*hx) - (b1*ht)/(4.0*hx);
    const double k14 = +(a2*ht)/(2.0*hy*hy) - (b2*ht)/(4.0*hy);
    const double k15 = -1.0 - (a2*ht)/(hy*hy) + 0.5*c*ht;
    const double k16 = +(a2*ht)/(2.0*hy*hy) + (b2*ht)/(4.0*hy);

    const double k21 = -(a2*ht)/(2.0*hy*hy) + (b2*ht)/(4.0*hy);
    const double k22 = -1.0 + (a2*ht)/(hy*hy) - 0.5*c*ht;
    const double k23 = -(a2*ht)/(2.0*hy*hy) - (b2*ht)/(4.0*hy);
    const double k24 = +(a1*ht)/(2.0*hx*hx) - (b1*ht)/(4.0*hx);
    const double k25 = -1.0 - (a1*ht)/(hx*hx);
    const double k26 = +(a1*ht)/(2.0*hx*hx) + (b1*ht)/(4.0*hx);

    // left border condition parameters
    const double b111 = -1.0 + ((a1*ht)/(hx*hx));
    const double b112 = +(a1*ht)/hx - 0.5*b1*ht;
    const double b113 = -(a1*ht)/(hx*hx);
    const double b114 = -1.0 + 0.5*c*ht;
    const double b115 = +0.5*a2*ht;
    const double b116 = +0.5*b2*ht;
    const double b117 = +(a1*ht)/hx - 0.5*b1*ht;

    // right border condition parameters
    const double b121 = -(a1*ht)/(hx*hx);
    const double b122 = +(a1*ht)/hx + 0.5*b1*ht;
    const double b123 = -1.0 + ((a1*ht)/(hx*hx));
    const double b124 = -1.0 + 0.5*c*ht;
    const double b125 = +0.5*a2*ht;
    const double b126 = +0.5*b2*ht;
    const double b127 = +(a1*ht)/hx + 0.5*b1*ht;

    // bottom border condition parameters
    const double b211 = -1.0 + ((a2*ht)/(hy*hy)) - 0.5*c*ht;
    const double b212 = +(a2*ht)/hy - 0.5*b2*ht;
    const double b213 = -(a2*ht)/(hy*hy);
    const double b214 = -1.0;
    const double b215 = +0.5*a1*ht;
    const double b216 = +0.5*b1*ht;
    const double b217 = +(a2*ht)/hy - 0.5*b2*ht;

    // top border condition parameters
    const double b221 = -(a2*ht)/(hy*hy);
    const double b222 = +(a2*ht)/hy + 0.5*b2*ht;
    const double b223 = -1.0 + ((a2*ht)/(hy*hy)) - 0.5*c*ht;
    const double b224 = -1.0;
    const double b225 = +0.5*a1*ht;
    const double b226 = +0.5*b1*ht;
    const double b227 = +(a2*ht)/hy + 0.5*b2*ht;

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

    double *ay = static_cast<double*>(malloc(sizeof(double)*(M+1)));
    double *by = static_cast<double*>(malloc(sizeof(double)*(M+1)));
    double *cy = static_cast<double*>(malloc(sizeof(double)*(M+1)));
    double *dy = static_cast<double*>(malloc(sizeof(double)*(M+1)));
    double *ry = static_cast<double*>(malloc(sizeof(double)*(M+1)));

    for (unsigned int m=0; m<=M; m++)
    {
        ay[m] = k21;
        by[m] = k22;
        cy[m] = k23;
    }
    ay[0] = 0.0; cy[M] = 0.0;

    DoubleMatrix u00(M+1, N+1);
    DoubleMatrix u05(M+1, N+1);
    DoubleMatrix u10(M+1, N+1);

    /***********************************************************************************************/
    /***********************************************************************************************/

    TimeNodePDE tn00, tn05, tn10;
    SpaceNodePDE sn;
    unsigned int i=0, j=0, s=0, e=N;
    int n = 0, m = 0;
    BoundaryConditionPDE condition; double value, alpha, beta, gamma;

    for (m=ymin, sn.j=m, sn.y=m*hy, j=0; m<=ymax; ++m, sn.j=m, sn.y=m*hy, ++j)
    {
        for (n=xmin, sn.i=n, sn.x=n*hx, i=0; n<=xmax; ++n, sn.i=n, sn.x=n*hx, ++i)
        {
            u00[j][i] = final(sn, FinalCondition::FinalValue);
        }
    }

    tn00.i = 2*L; tn00.t = 0.5*tn00.i*ht;
    layerInfo(u00, tn00);

    /***********************************************************************************************/

    for (unsigned int ln=L-1, size_ln = static_cast<unsigned int>(0)-1; ln != size_ln; ln--)
    {
        tn00.i = 2*ln+2; tn00.t = 0.5*tn00.i*ht;
        tn05.i = 2*ln+1; tn05.t = 0.5*tn05.i*ht;
        tn10.i = 2*ln-0; tn10.t = 0.5*tn10.i*ht;

        /**************************************************** x direction apprx ***************************************************/

        for (m=ymin+1, sn.j=m, sn.y=m*hy, j=1; m<=ymax-1; ++m, sn.j=m, sn.y=m*hy, ++j)
        {
            for (n=xmin+1, sn.i=n, sn.x=n*hx, i=1, dx[i]=0.0; n<=xmax-1; ++n, sn.i=n, sn.x=n*hx, ++i, dx[i]=0.0)
            {
                dx[i] += k14*u00[j-1][i] + k15*u00[j][i] + k16*u00[j+1][i] + ht05*f(sn, tn00);
            }

            sn.i = xmin; sn.x = xmin*hx;
            value = boundary(sn, tn05, condition);
            alpha = condition.alpha();
            beta  = condition.beta();
            gamma = condition.gamma();

            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                s = 1;
                u05[j][0] = (gamma/alpha)*value;
                dx[1] -= k11 * u05[j][0];
                ax[1] = ax[0] = bx[0] = cx[0] = dx[0] = rx[0] = 0.0;
            }
            if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
                s = 0;
                ax[s]  = 0.0;
                bx[s]  = beta  * b111;
                cx[s]  = beta  * b113;
                dx[s]  = beta  * b114 * u00[j][s];
                dx[s] += beta  * b115 * ((u00[j+1][s]-2.0*u00[j][s]+u00[j-1][s])/(hy*hy));
                dx[s] += beta  * b116 * ((u00[j+1][s]-u00[j-1][s])/(2.0*hy));
                dx[s] += beta  * ht05 * f(sn, tn00);
                dx[s] += gamma * b117 * value;
            }
            if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                s = 0;
                ax[s]  = 0.0;
                bx[s]  = beta  * b111 + alpha * b112;
                cx[s]  = beta  * b113;
                dx[s]  = beta  * b114 * u00[j][s];
                dx[s] += beta  * b115 * ((u00[j+1][s]-2.0*u00[j][s]+u00[j-1][s])/(hy*hy));
                dx[s] += beta  * b116 * ((u00[j+1][s]-u00[j-1][s])/(2.0*hy));
                dx[s] += beta  * ht05 * f(sn, tn00);
                dx[s] += gamma * b117 * value;
            }

            sn.i = xmax; sn.x = xmax*hx;
            value = boundary(sn, tn05, condition);
            alpha = condition.alpha();
            beta  = condition.beta();
            gamma = condition.gamma();

            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                e = N-1;
                u05[j][N] = (gamma/alpha)*value;
                dx[N-1] -= k13 * u05[j][N];
                cx[N-1] = ax[N] = bx[N] = cx[N] = dx[N] = rx[N] = 0.0;
            }
            if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
                e = N;
                ax[e]  = beta  * b121;
                bx[e]  = beta  * b123;
                cx[e]  = 0.0;
                dx[e]  = beta  * b124 * u00[j][e];
                dx[e] += beta  * b125 * ((u00[j+1][e]-2.0*u00[j][e]+u00[j-1][e])/(hy*hy));
                dx[e] += beta  * b126 * ((u00[j+1][e]-u00[j-1][e])/(2.0*hy));
                dx[e] += beta  * ht05 * f(sn, tn00);
                dx[e] += gamma * b127 * value;
            }
            if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                e = N;
                ax[e]  = beta  * b121;
                bx[e]  = beta  * b123 + alpha * b122;
                cx[e]  = 0.0;
                dx[e]  = beta  * b124 * u00[j][e];
                dx[e] += beta  * b125 * ((u00[j+1][e]-2.0*u00[j][e]+u00[j-1][e])/(hy*hy));
                dx[e] += beta  * b126 * ((u00[j+1][e]-u00[j-1][e])/(2.0*hy));
                dx[e] += beta  * ht05 * f(sn, tn00);
                dx[e] += gamma * b127 * value;
            }

            tomasAlgorithm(ax+s, bx+s, cx+s, dx+s, rx+s, e-s+1);
            for (unsigned int i=s; i<=e; i++) u05[j][i] = rx[i];
        }

        for (n=xmin, sn.i=n, sn.x=n*hx, i=0, dx[i]=0.0; n<=xmax; ++n, sn.i=n, sn.x=n*hx, ++i, dx[i]=0.0)
        {
            sn.j = ymin; sn.y = ymin*hy;
            value = boundary(sn, tn05, condition);
            alpha = condition.alpha();
            beta  = condition.beta();
            gamma = condition.gamma();

            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                u05[0][i] = (gamma/alpha)*value;
            }
            if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
                u05[0][i] = (3.5*u05[1][i] - 2.0*u05[2][i] + 0.5*u05[3][i] + hy*(gamma/beta)*value)/(2.0);
                //u05[0][i] = (u05[1][i] + hy*(gamma/beta)*value));
            }
            if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                u05[0][i] = (3.5*u05[1][i] - 2.0*u05[2][i] + 0.5*u05[3][i] + hy*(gamma/beta)*value)/(2.0 + (alpha/beta)*hy);
                //u05[0][i] = (u05[1][i] + hy*(gamma/beta)*value)/(1.0+hy*(alpha/beta));
            }

            sn.j = ymax; sn.y = ymax*hy;
            value = boundary(sn, tn05, condition);
            alpha = condition.alpha();
            beta  = condition.beta();
            gamma = condition.gamma();

            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                u05[M][i] = (gamma/alpha)*value;
            }

            if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                u05[M][i] = (3.5*u05[M-1][i] - 2.0*u05[M-2][i] + 0.5*u05[M-3][i] + hy*(gamma/beta)*value)/(2.0 + (alpha/beta)*hy);
                //u05[M][i] = (u05[M-1][i] + hy*(gamma/beta)*value)/(1.0 + hy*(alpha/beta));
            }
        }

        layerInfo(u05, tn05);

        /**************************************************** x direction apprx ***************************************************/

        /**************************************************** y direction apprx ***************************************************/

        for (n=xmin+1, sn.i=n, sn.x=n*hx, i=1; n<=xmax-1; ++n, sn.i=n, sn.x=n*hx, ++i)
        {
            for (m=ymin+1, sn.j=m, sn.y=m*hy, j=1, dy[j]=0.0; m<=ymax-1; ++m, sn.j=m, sn.y=m*hy, ++j, dy[j]=0.0)
            {
                dy[j] += k24*u05[j][i-1] + k25*u05[j][i] + k26*u05[j][i+1] + ht05*f(sn, tn10);
            }

            sn.j = ymin; sn.y = ymin*hy;
            value = boundary(sn, tn10, condition);
            alpha = condition.alpha();
            beta  = condition.beta();
            gamma = condition.gamma();

            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                s = 1;
                u10[0][i] = (gamma/alpha)*value;
                dy[1] -= k21 * u10[0][i];
                ay[1] = ay[0] = by[0] = cy[0] = dy[0] = ry[0] = 0.0;
            }
            if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
                s = 0;
                ay[s]  = 0.0;
                by[s]  = beta  * b211;
                cy[s]  = beta  * b213;
                dy[s]  = beta  * b214 * u05[s][i];
                dy[s] += beta  * b215 * ((u05[s][i+1]-2.0*u05[s][i]+u05[s][i-1])/(hx*hx));
                dy[s] += beta  * b216 * ((u05[s][i+1]-u05[s][i-1])/(2.0*hx));
                dy[s] += beta  * ht05 * f(sn, tn10);
                dy[s] += gamma * b217 * value;
            }
            if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                s = 0;
                ay[s]  = 0.0;
                by[s]  = beta  * b211 + alpha * b212;
                cy[s]  = beta  * b213;
                dy[s]  = beta  * b214 * u05[s][i];
                dy[s] += beta  * b215 * ((u05[s][i+1]-2.0*u05[s][i]+u05[s][i-1])/(hx*hx));
                dy[s] += beta  * b216 * ((u05[s][i+1]-u05[s][i-1])/(2.0*hx));
                dy[s] += beta  * ht05 * f(sn, tn10);
                dy[s] += gamma * b217 * value;
            }

            sn.j = ymax; sn.y = ymax*hy;
            value = boundary(sn, tn10, condition);
            alpha = condition.alpha();
            beta  = condition.beta();
            gamma = condition.gamma();

            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                e = M-1;
                u10[M][i] = (gamma/alpha)*value;
                dy[M-1] -= k23 * u10[M][i];
                cy[M-1] = ay[M] = by[M] = cy[M] = dy[M] = ry[M] = 0.0;
            }
            if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
                e = M;
                ay[e]  = beta  * b221;
                by[e]  = beta  * b223;
                cy[e]  = 0.0;
                dy[e]  = beta  * b224 * u05[e][i];
                dy[e] += beta  * b225 * ((u05[e][i+1]-2.0*u05[e][i]+u05[e][i-1])/(hx*hx));
                dy[e] += beta  * b226 * ((u05[e][i+1]-u05[e][i-1])/(2.0*hx));
                dy[e] += beta  * ht05 * f(sn, tn10);
                dy[e] += gamma * b227 * value;
            }
            if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                e = M;
                ay[e]  = beta  * b221;
                by[e]  = beta  * b223 + alpha * b222;
                cy[e]  = 0.0;
                dy[e]  = beta  * b224 * u05[e][i];
                dy[e] += beta  * b225 * ((u05[e][i+1]-2.0*u05[e][i]+u05[e][i-1])/(hx*hx));
                dy[e] += beta  * b226 * ((u05[e][i+1]-u05[e][i-1])/(2.0*hx));
                dy[e] += beta  * ht05 * f(sn, tn10);
                dy[e] += gamma * b227 * value;
            }

            tomasAlgorithm(ay+s, by+s, cy+s, dy+s, ry+s, e-s+1);
            for (unsigned int j=s; j<=e; j++) u10[j][i] = ry[j];
        }

        for (m=ymin, sn.j=m, sn.y=m*hy, j=0, dy[j]=0.0; m<=ymax; ++m, sn.j=m, sn.y=m*hy, ++j, dy[j]=0.0)
        {
            sn.i = xmin; sn.x = xmin*hx;
            value = boundary(sn, tn10, condition);
            alpha = condition.alpha();
            beta  = condition.beta();
            gamma = condition.gamma();

            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                u10[j][0] = (gamma/alpha)*value;
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
            value = boundary(sn, tn10, condition);
            alpha = condition.alpha();
            beta  = condition.beta();
            gamma = condition.gamma();

            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                u10[j][N] = (gamma/alpha)*value;
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

        layerInfo(u10, tn10);

        /**************************************************** y direction apprx ***************************************************/

        for (unsigned int j=0; j<=M; j++)
        {
            for (unsigned int i=0; i<=N; i++)
            {
                u00[j][i] = u10[j][i];
            }
        }
    }

    u00.clear();
    u05.clear();
    u10.clear();

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

void IHeatEquationFBVP::explicit_calculate_D2V1() const {}

double IHeatEquationFBVP::weight() const { return 0.5; }
