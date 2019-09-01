#include "pibvp.h"
#include <limits>
#include <exception>
#include <stdexcept>
#include "../cmethods.h"
#include "../linearequation.h"
#include "../printer.h"

IHeatEquationIBVP::IHeatEquationIBVP(double thermalDiffusivity) : _thermalDiffusivity(thermalDiffusivity) {}

IHeatEquationIBVP::IHeatEquationIBVP(const IHeatEquationIBVP &ibvp)
{
    this->_timeDimension = ibvp._timeDimension;
    this->_spaceDimensionX = ibvp._spaceDimensionX;
    this->_spaceDimensionY = ibvp._spaceDimensionY;
    this->_spaceDimensionZ = ibvp._spaceDimensionZ;
    this->_thermalDiffusivity = ibvp._thermalDiffusivity;
}

IHeatEquationIBVP & IHeatEquationIBVP::operator =(const IHeatEquationIBVP &other)
{
    if (this == &other) { return *this; }

    this->_timeDimension = other._timeDimension;
    this->_spaceDimensionX = other._spaceDimensionX;
    this->_spaceDimensionY = other._spaceDimensionY;
    this->_spaceDimensionZ = other._spaceDimensionZ;
    this->_thermalDiffusivity = other._thermalDiffusivity;
    return *this;
}

IHeatEquationIBVP::~IHeatEquationIBVP() {}

double IHeatEquationIBVP::thermalDiffusivity() const { return _thermalDiffusivity; }

void IHeatEquationIBVP::setThermalDiffusivity(double thermalDiffusivity)
{
    this->_thermalDiffusivity = thermalDiffusivity;
}

double IHeatEquationIBVP::thermalConductivity() const { return _thermalConductivity; }

void IHeatEquationIBVP::setThermalConductivity(double thermalConductivity)
{
    this->_thermalConductivity = thermalConductivity;
}

double IHeatEquationIBVP::thermalConvection() const
{
    return _thermalConvection;
}

void IHeatEquationIBVP::setThermalConvection(double thermalConvection)
{
    this->_thermalConvection = thermalConvection;
}

void IHeatEquationIBVP::explicit_calculate_D1V1() const
{}

void IHeatEquationIBVP::implicit_calculate_D1V1() const
{
    const unsigned int N = static_cast<unsigned int>( spaceDimensionX().size() );
    const unsigned int M = static_cast<unsigned int>( timeDimension().size() );

    const double hx = spaceDimensionX().step();
    const double ht = timeDimension().step();

    const double td = thermalDiffusivity();
    const double tc = thermalConductivity();


    const double m_td_ht__hxhx = -((td*ht)/(hx*hx));
    const double b_td_ht__hxhx__tcht = +(1.0 + (2.0*td*ht)/(hx*hx) + ht*tc);


    double *ax = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *bx = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *cx = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *dx = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *rx = static_cast<double*>(malloc(sizeof(double)*(N-1)));

    for (unsigned int n=0; n<=N-2; n++)
    {
        ax[n] = m_td_ht__hxhx;
        bx[n] = b_td_ht__hxhx__tcht;
        cx[n] = m_td_ht__hxhx;
    }
    ax[0] = 0.0; cx[N-2] = 0.0;

    DoubleVector u00(N+1);
    DoubleVector u10(N+1);

    /***********************************************************************************************/
    /***********************************************************************************************/

    TimeNodePDE tn;
    SpaceNodePDE sn;
    BoundaryConditionPDE condition;

    tn.i = 0; tn.t = tn.i*ht;
    for (unsigned int n=0; n<=N; n++)
    {
        sn.i = static_cast<int>(n); sn.x = sn.i*hx;
        u00[n] = initial(sn, InitialCondition::InitialValue);
    }
    layerInfo(u00, tn);

    /***********************************************************************************************/

    for (unsigned int ln=1; ln<=M; ln++)
    {
        tn.i = ln; tn.t = tn.i*ht;

        /**************************************************** border conditions ***************************************************/

        sn.i = static_cast<int>(0); sn.x = 0*hx; u10[0] = boundary(sn, tn, condition);
        sn.i = static_cast<int>(N); sn.x = N*hx; u10[N] = boundary(sn, tn, condition);

        /**************************************************** border conditions ***************************************************/

        /**************************************************** x direction apprx ***************************************************/
        for (unsigned int n=1; n<=N-1; n++)
        {
            sn.i = static_cast<int>(n); sn.x = n*hx;
            dx[n-1]  = u00[n];
            dx[n-1] += ht*f(sn, tn);
        }
        dx[0]   -= u10[0]*m_td_ht__hxhx;
        dx[N-2] -= u10[N]*m_td_ht__hxhx;
        tomasAlgorithm(ax, bx, cx, dx, rx, N-1);
        for (unsigned int n=1; n<=N-1; n++) u10[n] = rx[n-1];
        layerInfo(u10, tn);
        /**************************************************** x direction apprx ***************************************************/

        for (unsigned int n=0; n<=N; n++)
        {
            u00[n] = u10[n];
        }
    }

    u00.clear();
    u10.clear();

    free(rx);
    free(dx);
    free(cx);
    free(bx);
    free(ax);
}

void IHeatEquationIBVP::implicit_calculate_D1V1CN() const
{
    const unsigned int N = static_cast<unsigned int>(spaceDimensionX().size()) - 1;
    const unsigned int M = static_cast<unsigned int>(timeDimension().size());

    const double hx = spaceDimensionX().step();
    const double ht = timeDimension().step();

    const double td = thermalDiffusivity();
    const double tc = thermalConductivity();
    const double cv = thermalConvection();
    const double w1 = weight();
    const double w2 = 1.0 - weight();

    const double m_td_ht__hxhx_w1 = -((td*ht)/(hx*hx))*w1 + ((cv*ht)/(2.0*hx))*w1;
    const double b_td_ht__hxhx_w1_tc_ht_w1 = +(1.0 + ((2.0*td*ht)/(hx*hx))*w1 + tc*ht*w1);

    const double ht_tc_w2 = +ht*tc*w2;
    const double p_td_ht__hxhx_w2 = +((td*ht)/(hx*hx))*w2;
    const double m_cv_ht__20hx_w2 = -((cv*ht)/(2.0*hx))*w2;
    const double p_td_ht__hx_w1 = +((td*ht)/hx)*w1;

    double *ax = static_cast<double*>(malloc(sizeof(double)*(N+1)));
    double *bx = static_cast<double*>(malloc(sizeof(double)*(N+1)));
    double *cx = static_cast<double*>(malloc(sizeof(double)*(N+1)));
    double *dx = static_cast<double*>(malloc(sizeof(double)*(N+1)));
    double *rx = static_cast<double*>(malloc(sizeof(double)*(N+1)));

    for (unsigned int n=0; n<=N; n++)
    {
        ax[n] = m_td_ht__hxhx_w1;
        bx[n] = b_td_ht__hxhx_w1_tc_ht_w1;
        cx[n] = m_td_ht__hxhx_w1;
    }
    ax[0] = 0.0; cx[N] = 0.0;

    DoubleVector u00(N+1);
    DoubleVector u10(N+1);

    /***********************************************************************************************/
    /***********************************************************************************************/

    TimeNodePDE tn0; tn0.i = 0; tn0.t = tn0.i*ht;
    SpaceNodePDE sn;
    for (unsigned int n=0; n<=N; n++)
    {
        sn.i = static_cast<int>(n); sn.x = sn.i*hx;
        u00[n] = initial(sn, InitialCondition::InitialValue);
    }
    layerInfo(u00, tn0);

    /***********************************************************************************************/
    /***********************************************************************************************/

    for (unsigned int ln=1; ln<=M; ln++)
    {
        TimeNodePDE tn00; tn00.i = ln-1; tn00.t = tn00.i*ht;
        TimeNodePDE tn10; tn10.i = ln;   tn10.t = tn10.i*ht;

        for (unsigned int n=1; n<=N-1; n++)
        {
            sn.i = static_cast<int>(n); sn.x = n*hx;
            dx[n] = 0.0;
            dx[n] += u00[n] - u00[n] * ht_tc_w2;
            dx[n] += (u00[n-1] - 2.0*u00[n] + u00[n+1])*p_td_ht__hxhx_w2;
            dx[n] += (u00[n+1] - u00[n-1])*m_cv_ht__20hx_w2;
            dx[n] += ht*(w1*f(sn, tn10)+w2*f(sn, tn00));
        }

        unsigned int s=0, e=N;
        BoundaryConditionPDE condition; double alpha, beta, gamma, value;

        sn.i = static_cast<int>(0); sn.x = 0*hx;
        value = boundary(sn, tn10, condition);
        alpha = condition.alpha();
        beta  = condition.beta();
        gamma = condition.gamma();
        if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
        {
            s = 1;

            u10[0] = (gamma/alpha)*value;
            dx[1] -= u10[0]*m_td_ht__hxhx_w1;
            ax[1] = ax[0] = bx[0] = cx[0] = dx[0] = rx[0] = 0.0;
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
        {
            s = 0;

            ax[s]  = 0.0;
            bx[s]  = beta *(b_td_ht__hxhx_w1_tc_ht_w1);
            cx[s]  = beta *(2.0*m_td_ht__hxhx_w1);

            dx[s]  = u00[s]*(1.0 - ht_tc_w2);
            dx[s] += (2.0*u00[s]-5.0*u00[s+1]+4.0*u00[s+2]-u00[s+3])*p_td_ht__hxhx_w2;
            dx[s] += ht*(w1*f(sn, tn10)+w2*f(sn, tn00));
            dx[s] *= beta;

            dx[s] += gamma*(-2.0*p_td_ht__hx_w1)*value;
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Robin)
        {
            s = 0;

            ax[s]  = 0.0;
            bx[s]  = beta *(b_td_ht__hxhx_w1_tc_ht_w1);
            cx[s]  = beta *(2.0*m_td_ht__hxhx_w1);

            dx[s]  = u00[s]*(1.0 - ht_tc_w2);
            dx[s] += (2.0*u00[s]-5.0*u00[s+1]+4.0*u00[s+2]-u00[s+3])*p_td_ht__hxhx_w2;
            dx[s] += ht*(w1*f(sn, tn10)+w2*f(sn, tn00));
            dx[s] *= beta;

            bx[s] += alpha*(-2.0*p_td_ht__hx_w1);
            dx[s] += gamma*(-2.0*p_td_ht__hx_w1)*value;
        }

        sn.i = static_cast<int>(N); sn.x = N*hx;
        value = boundary(sn, tn10, condition);
        alpha = condition.alpha();
        beta  = condition.beta();
        gamma = condition.gamma();
        if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
        {
            e = N-1;
            u10[N] = (gamma/alpha)*value;
            dx[N-1] -= u10[N]*m_td_ht__hxhx_w1;
            cx[N-1] = ax[N] = bx[N] = cx[N] = dx[N] = rx[N] = 0.0;
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
        {
            e = N;

            ax[e]  = beta *(2.0*m_td_ht__hxhx_w1);
            bx[e]  = beta *(b_td_ht__hxhx_w1_tc_ht_w1);
            cx[e]  = 0.0;

            dx[e]  = u00[e]*(1.0 - ht_tc_w2);
            dx[e] += (-u00[e]+4.0*u00[e-1]-5.0*u00[e-2]+2.0*u00[e-3])*p_td_ht__hxhx_w2;
            dx[e] += ht*(w1*f(sn, tn10)+w2*f(sn, tn00));
            dx[e] *= beta;

            dx[e] += gamma*(+2.0*p_td_ht__hx_w1)*value;
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Robin)
        {
            e = N;

            ax[e]  = beta *(2.0*m_td_ht__hxhx_w1);
            bx[e]  = beta *(b_td_ht__hxhx_w1_tc_ht_w1);
            cx[e]  = 0.0;

            dx[e]  = u00[e]*(1.0 - ht_tc_w2);
            dx[e] += (-u00[e]+4.0*u00[e-1]-5.0*u00[e-2]+2.0*u00[e-3])*p_td_ht__hxhx_w2;
            dx[e] += ht*(w1*f(sn, tn10)+w2*f(sn, tn00));
            dx[e] *= beta;

            bx[e] += alpha*(+2.0*p_td_ht__hx_w1);
            dx[e] += gamma*(+2.0*p_td_ht__hx_w1)*value;
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

    if (ht > 0.5/(1.0/(hx*hx)+1.0/(hy*hy))) throw std::runtime_error("Differential scheme not steady");

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
    const unsigned int N = static_cast<unsigned int>( spaceDimensionX().size() );
    const unsigned int M = static_cast<unsigned int>( spaceDimensionY().size() );
    const unsigned int L = static_cast<unsigned int>( timeDimension().size() );

    const double hx = spaceDimensionX().step();
    const double hy = spaceDimensionY().step();
    const double ht = timeDimension().step();

    const double td = thermalDiffusivity();
    const double tc = thermalConductivity();


    const double ht_050 = ht*0.5;

    const double m_td_ht__hxhx_05 = -((0.5*td*ht)/(hx*hx));
    const double b_td_ht__hxhx_05tcht = +(1.0 + (td*ht)/(hx*hx)) + 0.5*ht*tc;
    const double p_td_ht__hyhy_05 = +((0.5*td*ht)/(hy*hy));


    const double m_td_ht__hyhy_05 = -((0.5*td*ht)/(hy*hy));
    const double b_td_ht__hyhy_05tcht = +(1.0 + (td*ht)/(hy*hy)) + 0.5*ht*tc;
    const double p_td_ht__hxhx_05 = +((0.5*td*ht)/(hx*hx));


    double *ax = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *bx = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *cx = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *dx = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *rx = static_cast<double*>(malloc(sizeof(double)*(N-1)));

    for (unsigned int n=0; n<=N-2; n++)
    {
        ax[n] = m_td_ht__hxhx_05;
        bx[n] = b_td_ht__hxhx_05tcht;
        cx[n] = m_td_ht__hxhx_05;
    }
    ax[0] = 0.0; cx[N-2] = 0.0;

    double *ay = static_cast<double*>(malloc(sizeof(double)*(M-1)));
    double *by = static_cast<double*>(malloc(sizeof(double)*(M-1)));
    double *cy = static_cast<double*>(malloc(sizeof(double)*(M-1)));
    double *dy = static_cast<double*>(malloc(sizeof(double)*(M-1)));
    double *ry = static_cast<double*>(malloc(sizeof(double)*(M-1)));

    for (unsigned int m=0; m<=M-2; m++)
    {
        ay[m] = m_td_ht__hyhy_05;
        by[m] = b_td_ht__hyhy_05tcht;
        cy[m] = m_td_ht__hyhy_05;
    }
    ay[0] = 0.0; cy[M-2] = 0.0;

    DoubleMatrix u00(M+1, N+1);
    DoubleMatrix u05(M+1, N+1);
    DoubleMatrix u10(M+1, N+1);

    /***********************************************************************************************/
    /***********************************************************************************************/

    TimeNodePDE tn00; tn00.i = 0; tn00.t = 0.5*tn00.i*ht;
    TimeNodePDE tn05; tn05.i = 1; tn05.t = 0.5*tn05.i*ht;

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

    for (unsigned int ln=1; ln<=L; ln++)
    {
        TimeNodePDE tn05; tn05.i = 2*ln-1; tn05.t = 0.5*tn05.i*ht;
        TimeNodePDE tn10; tn10.i = 2*ln-0; tn10.t = 0.5*tn10.i*ht;

        /**************************************************** border conditions ***************************************************/

        SpaceNodePDE sn0, sn1;
        BoundaryConditionPDE condition;

        sn0.i = static_cast<int>(0); sn0.x = 0*hx;
        sn1.i = static_cast<int>(N); sn1.x = N*hx;
        for (unsigned int m=0; m<=M; m++)
        {
            sn0.j = sn1.j = static_cast<int>(m); sn0.y = sn1.y = m*hy;
            u05[m][0] = boundary(sn0, tn05, condition); u10[m][0] = boundary(sn0, tn10, condition);
            u05[m][N] = boundary(sn1, tn05, condition); u10[m][N] = boundary(sn1, tn10, condition);
        }

        sn0.j = static_cast<int>(0); sn0.y = 0*hy;
        sn1.j = static_cast<int>(M); sn1.y = M*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn0.i = sn1.i = static_cast<int>(n); sn0.x = sn1.x = n*hx;
            u05[0][n] = boundary(sn0, tn05, condition); u10[0][n] = boundary(sn0, tn10, condition);
            u05[M][n] = boundary(sn1, tn05, condition); u10[M][n] = boundary(sn1, tn10, condition);
        }

        /**************************************************** border conditions ***************************************************/

        /**************************************************** x direction apprx ***************************************************/
        for (unsigned int m=1; m<=M-1; m++)
        {
            sn.j = static_cast<int>(m); sn.y = m*hy;
            for (unsigned int n=1; n<=N-1; n++)
            {
                sn.i = static_cast<int>(n); sn.x = n*hx;
                dx[n-1]  = u00[m][n];
                dx[n-1] += (u00[m-1][n] - 2.0*u00[m][n] + u00[m+1][n])*p_td_ht__hyhy_05;
                dx[n-1] += ht_050*f(sn, tn05);
            }
            dx[0]   -= u05[m][0]*m_td_ht__hxhx_05;
            dx[N-2] -= u05[m][N]*m_td_ht__hxhx_05;
            tomasAlgorithm(ax, bx, cx, dx, rx, N-1);
            for (unsigned int n=1; n<=N-1; n++) u05[m][n] = rx[n-1];
        }
        layerInfo(u05, tn05);
        /**************************************************** x direction apprx ***************************************************/

        /**************************************************** y direction apprx ***************************************************/
        for (unsigned int n=1; n<=N-1; n++)
        {
            sn.i = static_cast<int>(n); sn.x = n*hx;
            for (unsigned int m=1; m<=M-1; m++)
            {
                sn.j = static_cast<int>(m); sn.y = m*hy;
                dy[m-1]  = u05[m][n];
                dy[m-1] += (u05[m][n-1] - 2.0*u05[m][n] + u05[m][n+1])*p_td_ht__hxhx_05;
                dy[m-1] += ht_050*f(sn, tn10);
            }
            dy[0]   -= u10[0][n]*m_td_ht__hyhy_05;
            dy[M-2] -= u10[M][n]*m_td_ht__hyhy_05;
            tomasAlgorithm(ay, by, cy, dy, ry, M-1);
            for (unsigned int m=1; m<=M-1; m++) u10[m][n] = ry[m-1];
        }
        layerInfo(u10, tn10);
        /**************************************************** y direction apprx ***************************************************/

        for (unsigned int m=0; m<=M; m++)
        {
            for (unsigned int n=0; n<=N; n++)
            {
                u00[m][n] = u10[m][n];
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

void IHeatEquationIBVP::implicit_calculate_D2V1CN() const
{
    const unsigned int N = static_cast<unsigned int>( spaceDimensionX().size() );
    const unsigned int M = static_cast<unsigned int>( spaceDimensionY().size() );
    const unsigned int L = static_cast<unsigned int>( timeDimension().size() );

    const double hx = spaceDimensionX().step();
    const double hy = spaceDimensionY().step();
    const double ht = timeDimension().step();

    const double td = thermalDiffusivity();
    const double tc = thermalConductivity();
    const double _lambda = weight();

    if (_lambda >= 0.5-(0.25/ht)/(1.0/(hx*hx)+1.0/(hy*hy))) throw std::runtime_error("Differential scheme is conditionally steady.");

    const double ht_050 = ht*0.5;

    const double m_td_ht__hxhx_05_lambda = -(0.5*td*ht*_lambda)/(hx*hx);
    const double b_td_ht__hxhx_lambda__05tcht = +1.0 + (td*ht*_lambda)/(hx*hx) + 0.5*ht*tc;
    const double p_td_ht__hxhx_05_1mlambda = +(0.5*td*ht*(1.0-_lambda))/(hx*hx);
    const double p_td_ht__hyhy_05 = +((0.5*td*ht)/(hy*hy));

    const double m_td_ht__hyhy_05_lambda = -(0.5*td*ht*_lambda)/(hy*hy);
    const double b_td_ht__hyhy_lambda__conv = +(1.0 + (td*ht*_lambda)/(hy*hy)) + 0.5*ht*tc;
    const double p_td_ht__hyhy_05_1mlambda = +(0.5*td*ht*(1.0-_lambda))/(hy*hy);
    const double p_td_ht__hxhx_05 = +(0.5*td*ht)/(hx*hx);

    double *ax = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *bx = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *cx = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *dx = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *rx = static_cast<double*>(malloc(sizeof(double)*(N-1)));

    for (unsigned int n=0; n<=N-2; n++)
    {
        ax[n] = m_td_ht__hxhx_05_lambda;
        bx[n] = b_td_ht__hxhx_lambda__05tcht;
        cx[n] = m_td_ht__hxhx_05_lambda;
    }
    ax[0] = 0.0; cx[N-2] = 0.0;

    double *ay = static_cast<double*>(malloc(sizeof(double)*(M-1)));
    double *by = static_cast<double*>(malloc(sizeof(double)*(M-1)));
    double *cy = static_cast<double*>(malloc(sizeof(double)*(M-1)));
    double *dy = static_cast<double*>(malloc(sizeof(double)*(M-1)));
    double *ry = static_cast<double*>(malloc(sizeof(double)*(M-1)));

    for (unsigned int m=0; m<=M-2; m++)
    {
        ay[m] = m_td_ht__hyhy_05_lambda;
        by[m] = b_td_ht__hyhy_lambda__conv;
        cy[m] = m_td_ht__hyhy_05_lambda;
    }
    ay[0] = 0.0; cy[M-2] = 0.0;

    DoubleMatrix u00(M+1, N+1);
    DoubleMatrix u05(M+1, N+1);
    DoubleMatrix u10(M+1, N+1);

    /***********************************************************************************************/
    /***********************************************************************************************/

    TimeNodePDE tn00; tn00.i = 0; tn00.t = 0.5*tn00.i*ht;
    TimeNodePDE tn05; tn05.i = 1; tn05.t = 0.5*tn05.i*ht;

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

    for (unsigned int ln=1; ln<=L; ln++)
    {
        TimeNodePDE tn05; tn05.i = 2*ln-1; tn05.t = 0.5*tn05.i*ht;
        TimeNodePDE tn10; tn10.i = 2*ln-0; tn10.t = 0.5*tn10.i*ht;

        /**************************************************** border conditions ***************************************************/

        SpaceNodePDE sn0, sn1;
        BoundaryConditionPDE condition;

        sn0.i = static_cast<int>(0); sn0.x = 0*hx;
        sn1.i = static_cast<int>(N); sn1.x = N*hx;
        for (unsigned int m=0; m<=M; m++)
        {
            sn0.j = sn1.j = static_cast<int>(m); sn0.y = sn1.y = m*hy;
            u05[m][0] = boundary(sn0, tn05, condition); u10[m][0] = boundary(sn0, tn10, condition);
            u05[m][N] = boundary(sn1, tn05, condition); u10[m][N] = boundary(sn1, tn10, condition);
        }

        sn0.j = static_cast<int>(0); sn0.y = 0*hy;
        sn1.j = static_cast<int>(M); sn1.y = M*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn0.i = sn1.i = static_cast<int>(n); sn0.x = sn1.x = n*hx;
            u05[0][n] = boundary(sn0, tn05, condition); u10[0][n] = boundary(sn0, tn10, condition);
            u05[M][n] = boundary(sn1, tn05, condition); u10[M][n] = boundary(sn1, tn10, condition);
        }

        /**************************************************** border conditions ***************************************************/

        /**************************************************** x direction apprx ***************************************************/
        for (unsigned int m=1; m<=M-1; m++)
        {
            sn.j = static_cast<int>(m); sn.y = m*hy;
            for (unsigned int n=1; n<=N-1; n++)
            {
                sn.i = static_cast<int>(n); sn.x = n*hx;
                dx[n-1]  = u00[m][n];
                dx[n-1] += (u00[m][n-1] - 2.0*u00[m][n] + u00[m][n+1])*p_td_ht__hxhx_05_1mlambda;
                dx[n-1] += (u00[m-1][n] - 2.0*u00[m][n] + u00[m+1][n])*p_td_ht__hyhy_05;
                dx[n-1] += ht_050*f(sn, tn05);
                //dx[n-1] += ht_050*(_lambda*f(sn, tn05)+(1.0-_lambda)*f(sn, tn00));
            }
            dx[0]   -= u05[m][0]*m_td_ht__hxhx_05_lambda;
            dx[N-2] -= u05[m][N]*m_td_ht__hxhx_05_lambda;
            tomasAlgorithm(ax, bx, cx, dx, rx, N-1);
            for (unsigned int n=1; n<=N-1; n++) u05[m][n] = rx[n-1];
        }
        layerInfo(u05, tn05);
        /**************************************************** x direction apprx ***************************************************/

        /**************************************************** y direction apprx ***************************************************/
        for (unsigned int n=1; n<=N-1; n++)
        {
            sn.i = static_cast<int>(n); sn.x = n*hx;
            for (unsigned int m=1; m<=M-1; m++)
            {
                sn.j = static_cast<int>(m); sn.y = m*hy;
                dy[m-1]  = u05[m][n];
                dy[m-1] += (u05[m-1][n] - 2.0*u05[m][n] + u05[m+1][n])*p_td_ht__hyhy_05_1mlambda;
                dy[m-1] += (u05[m][n-1] - 2.0*u05[m][n] + u05[m][n+1])*p_td_ht__hxhx_05;
                dy[m-1] += ht_050*f(sn, tn10);
                //dy[m-1] += ht_050*(_lambda*f(sn, tn10)+(1.0-_lambda)*f(sn, tn05));
            }
            dy[0]   -= u10[0][n]*m_td_ht__hyhy_05_lambda;
            dy[M-2] -= u10[M][n]*m_td_ht__hyhy_05_lambda;
            tomasAlgorithm(ay, by, cy, dy, ry, M-1);
            for (unsigned int m=1; m<=M-1; m++) u10[m][n] = ry[m-1];
        }
        layerInfo(u10, tn10);
        /**************************************************** y direction apprx ***************************************************/

        for (unsigned int m=0; m<=M; m++)
        {
            for (unsigned int n=0; n<=N; n++)
            {
                u00[m][n] = u10[m][n];
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

void IHeatEquationIBVP::calculateU(DoubleMatrix &u, double a, double alpha, double lambda)
{
    const unsigned int N = static_cast<unsigned int>(_spaceDimensionX.size());
    const unsigned int M = static_cast<unsigned int>(_spaceDimensionY.size());
    const unsigned int L = static_cast<unsigned int>(_timeDimension.size());
    const double hx = _spaceDimensionX.step();
    const double hy = _spaceDimensionY.step();
    const double ht = _timeDimension.step();

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
            u[m][n] = initial(sn, InitialCondition::InitialValue);
        }
    }
    //layerInfo(u, 0);
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
                    a1X[0] = +0.0;
                    b1X[0] = +2.0 + 2.0*((a*a*ht)/(hx*hx)) + alpha*ht - 2.0*(a*a*lambda*ht/hx);
                    c1X[0] = -2.0*((a*a*ht)/(hx*hx));
                    d1X[0] -= 2.0*(a*a*lambda*ht/hx)*env1(sn, tn);
                }
                else if (n == N)
                {
                    a1X[N] = -2.0*((a*a*ht)/(hx*hx));
                    b1X[N] = +2.0 + 2.0*((a*a*ht)/(hx*hx)) + alpha*ht - 2.0*(a*a*lambda*ht/hx);
                    c1X[N] = +0.0;
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
        //for (unsigned int n=0; n<=N; n++) uh[0][n] = n*hx*n*hx + 0.0*0.0 + tn.t;
        //for (unsigned int n=0; n<=N; n++) uh[M][n] = n*hx*n*hx + 1.0*1.0 + tn.t;
        if (l==1)
        {
            IPrinter::printMatrix(uh);
            IPrinter::printSeperatorLine();
        }return;
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
                    a1Y[0] = +0.0;
                    b1Y[0] = +2.0 + 2.0*((a*a*ht)/(hy*hy)) + alpha*ht - 2.0*(a*a*lambda*ht/hy);
                    c1Y[0] = -2.0*((a*a*ht)/(hy*hy));
                    d1Y[0] -= 2.0*(a*a*lambda*ht/hy)*env1(sn, tn);
                }
                else if (m == M)
                {
                    a1Y[M] = -2.0*((a*a*ht)/(hy*hy));
                    b1Y[M] = +2.0 + 2.0*((a*a*ht)/(hy*hy)) + alpha*ht - 2.0*(a*a*lambda*ht/hy);
                    c1Y[M] = +0.0;
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
        //for (unsigned int m=0; m<=M; m++) u[m][0] = hy*m*hy*m + 0.0*0.0 + tn.t;
        //for (unsigned int m=0; m<=M; m++) u[m][N] = hy*m*hy*m + 1.0*1.0 + tn.t;
        //------------------------------------- approximatin to y direction conditions -------------------------------------//

        //layerInfo(u, l);
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

void IHeatEquationIBVP::gridMethod(DoubleVector &u, double a) const
{
    const Dimension &dimX = spaceDimensionX();//spaceDimension(Dimension::DimensionX);
    const Dimension &time = timeDimension();

    const double ht = time.step();
    const double hx = dimX.step();

    const int minM = time.min();
    const int maxM = time.max();
    const unsigned int M = static_cast<unsigned int>(maxM-minM);

    const int minN = dimX.min();
    const int maxN = dimX.max();
    const unsigned int N = static_cast<unsigned int>(maxN-minN);

    const double alpha = -(a*a*ht)/(hx*hx);
    const double betta = 1.0 - 2.0*alpha;

    u.clear();
    u.resize(N+1);

    double *ka = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *kb = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *kc = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *kd = static_cast<double*>(malloc(sizeof(double)*(N-1)));
    double *rx = static_cast<double*>(malloc(sizeof(double)*(N-1)));

    /* initial condition */
    SpaceNodePDE isn;
    unsigned int i = 0;
    for (int n=minN; n<=maxN; n++, i++)
    {
        isn.i = n;
        isn.x = n*hx;
        u[i] = initial(isn, InitialCondition::InitialValue);
    }
    //layerInfo(u, 0);

    SpaceNodePDE lsn; lsn.i = minN; lsn.x = minN*hx;
    SpaceNodePDE rsn; rsn.i = maxN; rsn.x = maxN*hx;
    BoundaryConditionPDE condition;

    TimeNodePDE tn;
    for (unsigned int m=1; m<=M; m++)
    {
        tn.i = m+static_cast<unsigned int>(minM);
        tn.t = tn.i*ht;

        unsigned int i = 1;
        for (int n=minN+1; n<=maxN-1; n++, i++)
        {
            isn.i = n;
            isn.x = n*hx;

            ka[i-1] = alpha;
            kb[i-1] = betta;
            kc[i-1] = alpha;
            kd[i-1] = u[i] + ht * f(isn, tn);
        }

        ka[0]   = 0.0;
        kc[N-2] = 0.0;

        /* border conditions */
        u[0] = boundary(lsn, tn, condition);
        u[N] = boundary(rsn, tn, condition);

        kd[0]   -= alpha * u[0];
        kd[N-2] -= alpha * u[N];

        tomasAlgorithm(ka, kb, kc, kd, rx, N-1);

        for (unsigned int n=1; n<=N-1; n++) u[n] = rx[n-1];

        //layerInfo(u, m);
    }

    free(ka);
    free(kb);
    free(kc);
    free(kd);
    free(rx);
}

void IHeatEquationIBVP::calculateMVD(DoubleMatrix &u) const
{
    const Dimension &dimX = spaceDimensionX();//mspaceDimension.at(Dimension::DimensionX);
    const Dimension &dimY = spaceDimensionY();//mspaceDimension.at(Dimension::DimensionY);
    const Dimension &time = timeDimension();

    const double ht = time.step();
    const double h1 = dimX.step();
    const double h2 = dimY.step();

    const unsigned int M = time.max();
    const unsigned int N1 = dimX.max();
    const unsigned int N2 = dimY.max();

    double a1 = 1.0;
    double a2 = 1.0;

    //cleaning matrix
    u.clear();
    u.resize(N2+1, N1+1);

    DoubleMatrix v(N2+1, N1+1);

    double* da1 = static_cast<double*>(malloc(sizeof(double)*(N1-1)));
    double* db1 = static_cast<double*>(malloc(sizeof(double)*(N1-1)));
    double* dc1 = static_cast<double*>(malloc(sizeof(double)*(N1-1)));
    double* dd1 = static_cast<double*>(malloc(sizeof(double)*(N1-1)));
    double* rx1 = static_cast<double*>(malloc(sizeof(double)*(N1-1)));

    double* da2 = static_cast<double*>(malloc(sizeof(double)*(N2-1)));
    double* db2 = static_cast<double*>(malloc(sizeof(double)*(N2-1)));
    double* dc2 = static_cast<double*>(malloc(sizeof(double)*(N2-1)));
    double* dd2 = static_cast<double*>(malloc(sizeof(double)*(N2-1)));
    double* rx2 = static_cast<double*>(malloc(sizeof(double)*(N2-1)));

    double x1_a = -(a1*ht)/(2.0*h1*h1);
    double x1_b  = 1.0 + (a1*ht)/(h1*h1);
    double x1_c = (a2*ht)/(2.0*h2*h2);

    double x2_a = -(a2*ht)/(2.0*h2*h2);
    double x2_b  = 1.0 + (a2*ht)/(h2*h2);
    double x2_c = (a1*ht)/(2.0*h1*h1);

    for (unsigned int j=0; j<=N2; j++)
    {
        SpaceNodePDE sn;
        sn.j = j;
        sn.y = j*h2;
        for (unsigned int i=0; i<=N1; i++)
        {
            sn.i = i;
            sn.x = i*h1;
            u[j][i] = initial(sn, InitialCondition::InitialValue);
        }
    }

    TimeNodePDE tn;
    SpaceNodePDE sn;
    BoundaryConditionPDE condition;
    for (unsigned int k=1; k<=M; k++)
    {
        tn.i = 2*k-1;
        tn.t = 0.5*(2*k-1)*ht;

        // Approximation to x1 direction
        for (unsigned int j=1; j<N2; j++)
        {
            sn.j = j;
            sn.y = j*h2;
            for (unsigned int i=1; i<N1; i++)
            {
                sn.i = i;
                sn.x = i*h1;

                da1[i-1] = x1_a;
                db1[i-1] = x1_b;
                dc1[i-1] = x1_a;
                dd1[i-1] = x1_c*(u[j-1][i] - 2.0*u[j][i] + u[j+1][i]) + u[j][i] + 0.5*ht * f(sn, tn);
            }

            da1[0]     = 0.0;
            dc1[N1-2]  = 0.0;

            SpaceNodePDE sn0; sn0.i = 0;  sn0.x = 0.0;   sn0.j = j; sn0.y = j*h2;
            SpaceNodePDE snN; snN.i = N1; snN.x = N1*h1; snN.j = j; snN.y = j*h2;
            v[j][0]  = boundary(sn0, tn, condition);
            v[j][N1] = boundary(snN, tn, condition);

            dd1[0]    -= x1_a * v[j][0];
            dd1[N1-2] -= x1_a * v[j][N1];

            tomasAlgorithm(da1, db1, dc1, dd1, rx1, N1-1);

            for (unsigned int i=1; i<N1; i++) v[j][i] = rx1[i-1];
        }

        for (unsigned int i=0; i<=N1; i++)
        {
            SpaceNodePDE sn0; sn0.i = i; sn0.x = i*h1; sn0.j = 0;  sn0.y = 0.0;
            SpaceNodePDE snN; snN.i = i; snN.x = i*h1; snN.j = N2; snN.y = N2*h2;
            v[0][i]  = boundary(sn0, tn, condition);
            v[N2][i] = boundary(snN, tn, condition);
        }

        tn.i = 2*k;
        tn.t = k*ht;
        // Approximation to x2 direction
        for (unsigned int i=1; i<N1; i++)
        {
            sn.i = i;
            sn.x = i*h1;
            for (unsigned int j=1; j<N2; j++)
            {
                sn.j = j;
                sn.y = j*h2;

                da2[j-1] = x2_a;
                db2[j-1] = x2_b;
                dc2[j-1] = x2_a;
                dd2[j-1] = x2_c*(v[j][i-1] - 2.0*v[j][i] + v[j][i+1]) + v[j][i] + 0.5*ht * f(sn, tn);
            }
            da2[0]     = 0.0;
            dc2[N2-2]  = 0.0;

            SpaceNodePDE sn0; sn0.i = i; sn0.x = i*h1; sn0.j = 0;  sn0.y = 0.0;
            SpaceNodePDE snN; snN.i = i; snN.x = i*h1; snN.j = N2; snN.y = N2*h2;
            u[0][i]  = boundary(sn0, tn, condition);
            u[N2][i] = boundary(snN, tn, condition);

            dd2[0]    -= x2_a * u[0][i];
            dd2[N2-2] -= x2_a * u[N2][i];

            tomasAlgorithm(da2, db2, dc2, dd2, rx2, N2-1);

            for (unsigned int j=1; j<N2; j++) u[j][i] = rx2[j-1];
        }

        for (unsigned int j=0; j<=N2; j++)
        {
            SpaceNodePDE sn0; sn0.i = 0;  sn0.x = 0*h1;  sn0.j = j; sn0.y = j*h2;
            SpaceNodePDE snN; snN.i = N1; snN.x = N1*h1; snN.j = j; snN.y = j*h2;
            u[j][0]  = boundary(sn0, tn, condition);
            u[j][N1] = boundary(snN, tn, condition);
        }
    }

    free(rx2);
    free(dd2);
    free(dc2);
    free(db2);
    free(da2);

    free(rx1);
    free(dd1);
    free(dc1);
    free(db1);
    free(da1);
}

void funcL(const double* a, const double *b, const double *c, const double *d, double *x, unsigned int N);

void funcL(const double* a, const double *b, const double *c, const double *d, double *x, unsigned int N)
{
    double *e = (double*) malloc(sizeof(double)*N);
    for (unsigned int i=0; i<N; i++) e[i] = 0.0;

    double *e0 = (double*)malloc(sizeof(double)*N);
    double *e1 = (double*)malloc(sizeof(double)*N);
    double *e2 = (double*)malloc(sizeof(double)*N);
    for (unsigned int i=0; i<N; i++) e0[i] = e1[i] = e2[i] = 0.0;

    e0[0] = b[0];
    e1[0] = c[0];
    e2[0] = d[0];
    for (unsigned int n=1; n<=N-2; n++)
    {
        e0[n] = -e0[n-1]*(b[n]/a[n]) + e1[n-1];
        e1[n] = -e0[n-1]*(c[n]/a[n]);
        e2[n] = -e0[n-1]*(d[n]/a[n]) + e2[n-1];
    }

    DoubleMatrix M(2,2);
    DoubleVector A(2);
    DoubleVector y(2);

    M[0][0] = e0[N-2];   M[0][1] = e1[N-2];   A[0] = e2[N-2];
    M[1][0] = a[N-1];    M[1][1] = b[N-1];    A[1] = d[N-1];

    LinearEquation::GaussianElimination(M,A,y);

    x[N-1] = y[1];
    x[N-2] = y[0];
    for (unsigned int n=N-3; n!=UINT_MAX; n--)
    {
        x[n] = (e2[n] - e1[n]*x[n+1])/e0[n];
    }

    free(e2);
    free(e1);
    free(e0);
}

void IHeatEquationIBVP::calculateN2L2RD(DoubleMatrix &u) const
{
    double ht = _timeDimension.step();
    unsigned int minM = _timeDimension.min();
    unsigned int maxM = _timeDimension.max();
    unsigned int M = maxM - minM;

    double hx = _spaceDimensionX.step();
    unsigned int minN = _spaceDimensionX.min();
    unsigned int maxN = _spaceDimensionX.max();
    unsigned int N = maxN - minN;

    const unsigned int k = 2;
    double h = ht/(hx*hx);

    u.clear();
    u.resize(M+1, N+1);

    DoubleMatrix A(k, k, 0.0);
    DoubleVector b(k, 0.0);
    DoubleVector x(k, 0.0);
    DoubleMatrix ems(N-k, k);

    /* initial condition */
    SpaceNodePDE isn;
    BoundaryConditionPDE condition;
    for (unsigned int n=0; n<=N; n++)
    {
        isn.i = n+minN;
        isn.x = isn.i*hx;
        u[0][n] = initial(isn, InitialCondition::InitialValue);
    }

    TimeNodePDE tn;
    SpaceNodePDE lsn;
    SpaceNodePDE rsn;
    lsn.i = minN; lsn.x = minN*hx;
    rsn.i = maxN; rsn.x = maxN*hx;
    for (unsigned int m=1; m<=M; m++)
    {
        tn.i = m+minM;
        tn.t = tn.i*ht;

        /* border conditions */
        u[m][0] = boundary(lsn, tn, condition);
        u[m][N] = boundary(rsn, tn, condition);

        /* n=1 */
        isn.i = minN+1;
        isn.x = isn.i*hx;

        double alpha = thermalDiffusivity()*h;
        A[0][0] = -2.0*alpha - 1.0;
        A[0][1] = alpha;
        b[0]    = -u[m-1][1] - alpha*u[m][0] - ht*f(isn,tn);

        A[0][1] /= A[0][0];
        b[0]    /= A[0][0];
        A[0][0] = 1.0;

        ems[0][0] = A[0][1];
        ems[0][1] = b[0];

        for (unsigned int n=2; n<=N-k; n++)
        {
            isn.i = n+minN;
            isn.x = isn.i*hx;

            alpha = thermalDiffusivity()*h;

            double g1 = alpha;
            double g2 = -2.0*alpha-1.0;
            double g3 = alpha;
            double fi = -u[m-1][n] - ht*f(isn,tn);

            g2 /= -g1;
            g3 /= -g1;
            fi /= +g1;
            g1  = 1.0;

            A[0][0] = A[0][1] + g2;
            A[0][1] = g3;
            b[0]    = b[0] - fi;
            \
            A[0][1] /= A[0][0];
            b[0]    /= A[0][0];
            A[0][0] = 1.0;

            ems[n-1][0] = A[0][1];
            ems[n-1][1] = b[0];
        }

        isn.i = maxN-1;
        isn.x = isn.i*hx;
        alpha = thermalDiffusivity()*h;
        A[1][0] = alpha;
        A[1][1] = -2.0*alpha - 1.0;
        b[1]    = -u[m-1][N-1] - alpha*u[m][N] - ht*f(isn,tn);

        LinearEquation::GaussianElimination(A, b, x);

        u[m][N-1] = x[1];
        u[m][N-2] = x[0];
        for (unsigned int n=N-(k+1); n>=1; n--)
        {
            u[m][n] = -ems[n-1][0]*u[m][n+1]+ems[n-1][1];
        }
    }

    ems.clear();
    x.clear();
    b.clear();
    A.clear();
}

#define SCHEME_2
#define __NORMALIZE__X
void IHeatEquationIBVP::calculateN4L2RD(DoubleMatrix &u) const
{
    /* get parameters */
    double ht = _timeDimension.step();
    unsigned int minM = _timeDimension.min();
    unsigned int maxM = _timeDimension.max();
    unsigned int M = maxM - minM;

    double hx = _spaceDimensionX.step();
    unsigned int minN = _spaceDimensionX.min();
    unsigned int maxN = _spaceDimensionX.max();
    unsigned int N = maxN - minN;

    const unsigned int k = 4;
    double h = ht/(24.0*hx*hx);

    /*****************************************/

    u.clear();
    u.resize(M+1, N+1);

    DoubleMatrix A(k, k, 0.0);
    DoubleVector b(k, 0.0);
    DoubleVector x(k, 0.0);
    DoubleMatrix ems(N-k, k+1);

    /* initial condition */
    SpaceNodePDE isn;
    for (unsigned int n=0; n<=N; n++)
    {
        isn.i = n+minN;
        isn.x = isn.i*hx;
        u[0][n] = initial(isn, InitialCondition::InitialValue);
    }

    TimeNodePDE tn;
    SpaceNodePDE lsn;
    SpaceNodePDE rsn;
    lsn.i = minN; lsn.x = minN*hx;
    rsn.i = maxN; rsn.x = maxN*hx;
    BoundaryConditionPDE condition;
    for (unsigned int m=1; m<=M; m++)
    {
        tn.i = m+minM;
        tn.t = tn.i*ht;

        /* border conditions */
        u[m][0] = boundary(lsn, tn, condition);
        u[m][N] = boundary(rsn, tn, condition);

        /* using 2nd scheme, at point n=1 Березин И.С., Жидков Н.П. - Методы вычислений (том 1) */
        isn.i = minN+1;
        isn.x = isn.i*hx;
        double alpha = thermalDiffusivity()*h;
        A[0][0] = -40.0*alpha - 1.0;
        A[0][1] = +12.0*alpha;
        A[0][2] = +8.0*alpha;
        A[0][3] = -2.0*alpha;
        b[0]    = -u[m-1][1] - (+22.0*alpha)*u[m][0] - ht*f(isn,tn);

        /* using 5th scheme, at point n=4 Березин И.С., Жидков Н.П. - Методы вычислений (том 1) */
        //isn.i = minN+4;
        //isn.x = isn.i*hx;
        //double alpha = a(isn,tn)*h;
        //A[0][0] = -112.0*alpha;
        //A[0][1] = +228.0*alpha;
        //A[0][2] = -208.0*alpha;
        //A[0][3] = +70.0*alpha - 1.0;
        //b[0]    = -u[m-1][4] - (+22.0*alpha)*u[m][0] - ht*f(isn,tn);

        //double x1 = (1)*hx;
        //double x2 = (2)*hx;
        //double x3 = (3)*hx;
        //double x4 = (4)*hx;
        //printf("%4d %18.10f %18.10f\n", 0, b[0], A[0][0]*x1*x1*x1*tn.t+A[0][1]*x2*x2*x2*tn.t+A[0][2]*x3*x3*x3*tn.t+A[0][3]*x4*x4*x4*tn.t);

#ifdef __NORMALIZE__
        A[0][1] /= A[0][0];
        A[0][2] /= A[0][0];
        A[0][3] /= A[0][0];
        b[0]    /= A[0][0];
        A[0][0] = 1.0;
#endif

        ems[0][0] = A[0][0];
        ems[0][1] = A[0][1];
        ems[0][2] = A[0][2];
        ems[0][3] = A[0][3];
        ems[0][4] = b[0];

#ifdef SCHEME_1
        unsigned int start = 1;
        unsigned int end   = N-k-1;
        unsigned int s     = 0;
#endif
#ifdef SCHEME_2
        unsigned int start = 2;
        unsigned int end   = N-k+0;
        unsigned int s     = 1;
#endif
#ifdef SCHEME_3
        unsigned int start = 3;
        unsigned int end   = N-k+1;
        unsigned int s     = 2;
#endif
#ifdef SCHEME_4
        unsigned int start = 4;
        unsigned int end   = N-k+2;
        unsigned int s     = 3;
#endif
#ifdef SCHEME_5
        unsigned int start = 5;
        unsigned int end   = N-k+3;
        unsigned int s     = 4;
#endif
        for (unsigned int n=start; n<=end; n++)
        {
            isn.i = n+minN;
            isn.x = isn.i*hx;

            alpha = thermalDiffusivity()*h;

#ifdef SCHEME_1
            /* using 1nd scheme, at point n=1 Березин И.С., Жидков Н.П. - Методы вычислений (том 1) */
            double g1 = +70.0*alpha-1.0;
            double g2 = -208.0*alpha;
            double g3 = +228.0*alpha;
            double g4 = -112.0*alpha;
            double g5 = +22.0*alpha;
            double fi = -u[m-1][n] - ht*f(isn,tn);
#endif
#ifdef SCHEME_2
            /* using 2nd scheme, at point n=2 Березин И.С., Жидков Н.П. - Методы вычислений (том 1) */
            double g1 = +22.0*alpha;
            double g2 = -40.0*alpha-1.0;
            double g3 = +12.0*alpha;
            double g4 = +8.0*alpha;
            double g5 = -2.0*alpha;
            double fi = -u[m-1][n] - ht*f(isn,tn);
#endif
#ifdef SCHEME_3
            /* using 3nd scheme, at point n=3 Березин И.С., Жидков Н.П. - Методы вычислений (том 1) */
            double g1 = -2.0*alpha;
            double g2 = +32.0*alpha;
            double g3 = -60.0*alpha-1.0;
            double g4 = +32.0*alpha;
            double g5 = -2.0*alpha;
            double fi = -u[m-1][n] - ht*f(isn,tn);
#endif
#ifdef SCHEME_4
            /* using 4nd scheme, at point n=4 Березин И.С., Жидков Н.П. - Методы вычислений (том 1) */
            double g1 = -2.0*alpha;
            double g2 = +8.0*alpha;
            double g3 = +12.0*alpha;
            double g4 = -40.0*alpha-1.0;
            double g5 = +22.0*alpha;
            double fi = -u[m-1][n] - ht*f(isn,tn);
#endif
#ifdef SCHEME_5
            /* using 5nd scheme, at point n=5 Березин И.С., Жидков Н.П. - Методы вычислений (том 1) */
            double g1 = +22.0*alpha;
            double g2 = -112.0*alpha;
            double g3 = +228.0*alpha;
            double g4 = -208.0*alpha;
            double g5 = +70.0*alpha-1.0;
            double fi = -u[m-1][n] - ht*f(isn,tn);
#endif

            g2 /= -g1;
            g3 /= -g1;
            g4 /= -g1;
            g5 /= -g1;
            fi /= +g1;
            g1 = 1.0;

            double A00 = A[0][0];
            A[0][0] = A[0][1] + g2*A00;
            A[0][1] = A[0][2] + g3*A00;
            A[0][2] = A[0][3] + g4*A00;
            A[0][3] = g5*A00;
            b[0]    = b[0] - fi*A00;

            //double x1 = (n+0)*hx;
            //double x2 = (n+1)*hx;
            //double x3 = (n+2)*hx;
            //double x4 = (n+3)*hx;
            //printf("%4d %18.10f %18.10f\n", n-s, b[0], A[0][0]*x1*x1*x1*tn.t+A[0][1]*x2*x2*x2*tn.t+A[0][2]*x3*x3*x3*tn.t+A[0][3]*x4*x4*x4*tn.t);

#ifdef __NORMALIZE__
            A[0][1] /= A[0][0];
            A[0][2] /= A[0][0];
            A[0][3] /= A[0][0];
            b[0]    /= A[0][0];
            A[0][0] = 1.0;
#endif

            // printf("0 %14.10f %14.10f\n", b[0], A[0][0]*x1*x1*t+A[0][1]*x2*x2*t+A[0][2]*x3*x3*t+A[0][3]*x4*x4*t);

            ems[n-s][0] = A[0][0];
            ems[n-s][1] = A[0][1];
            ems[n-s][2] = A[0][2];
            ems[n-s][3] = A[0][3];
            ems[n-s][4] = b[0];
        }

        /* using 2nd scheme, at point N-3 Березин И.С., Жидков Н.П. - Методы вычислений (том 1) */
        isn.i = maxN-3;
        isn.x = isn.i*hx;
        alpha = thermalDiffusivity()*h;
        A[1][0] = +22.0*alpha;
        A[1][1] = -40.0*alpha - 1.0;
        A[1][2] = +12.0*alpha;
        A[1][3] = +8.0*alpha;
        b[1]    = -u[m-1][N-3] - (-2.0*alpha)*u[m][N] - ht*f(isn,tn);

        /* using 3rd scheme, at point N-2 Березин И.С., Жидков Н.П. - Методы вычислений (том 1) */
        isn.i = maxN-2;
        isn.x = isn.i*hx;
        alpha = thermalDiffusivity()*h;
        A[2][0] = -2.0*alpha;
        A[2][1] = +32.0*alpha;
        A[2][2] = -60.0*alpha - 1.0;
        A[2][3] = +32.0*alpha;
        b[2]    = -u[m-1][N-2] - (-2.0*alpha)*u[m][N] - ht*f(isn,tn);

        /* using 4th scheme, at point N-1 Березин И.С., Жидков Н.П. - Методы вычислений (том 1) */
        isn.i = maxN-1;
        isn.x = isn.i*hx;
        A[3][0] = -2.0*alpha;
        A[3][1] = +8.0*alpha;
        A[3][2] = +12.0*alpha;
        A[3][3] = -40.0*alpha - 1.0;
        b[3]    = -u[m-1][N-1] - (+22.0*alpha)*u[m][N] - ht*f(isn,tn);

        //double xN4 = (N-4)*hx;
        //double xN3 = (N-3)*hx;
        //double xN2 = (N-2)*hx;
        //double xN1 = (N-1)*hx;

        //puts("------------------------");
        //printf("0 %18.10f %18.10f\n", b[0], A[0][0]*xN4*xN4*xN4*tn.t+A[0][1]*xN3*xN3*xN3*tn.t+A[0][2]*xN2*xN2*xN2*tn.t+A[0][3]*xN1*xN1*xN1*tn.t);
        //printf("1 %18.10f %18.10f\n", b[1], A[1][0]*xN4*xN4*xN4*tn.t+A[1][1]*xN3*xN3*xN3*tn.t+A[1][2]*xN2*xN2*xN2*tn.t+A[1][3]*xN1*xN1*xN1*tn.t);
        //printf("2 %18.10f %18.10f\n", b[2], A[2][0]*xN4*xN4*xN4*tn.t+A[2][1]*xN3*xN3*xN3*tn.t+A[2][2]*xN2*xN2*xN2*tn.t+A[2][3]*xN1*xN1*xN1*tn.t);
        //printf("3 %18.10f %18.10f\n", b[3], A[3][0]*xN4*xN4*xN4*tn.t+A[3][1]*xN3*xN3*xN3*tn.t+A[3][2]*xN2*xN2*xN2*tn.t+A[3][3]*xN1*xN1*xN1*tn.t);
        //puts("------------------------");

        LinearEquation::GaussianElimination(A, b, x);

        //printf("x %14.10f %14.10f %14.10f %14.10f\n", x[0], x[1], x[2], x[3]);

        u[m][N-1] = x[3];
        u[m][N-2] = x[2];
        u[m][N-3] = x[1];
        u[m][N-4] = x[0];

        for (unsigned int n=end-s; n>=start-s; n--) u[m][n] = -(ems[n-1][1]*u[m][n+1] + ems[n-1][2]*u[m][n+2] + ems[n-1][3]*u[m][n+3] - ems[n-1][4]) / ems[n-1][0];

        //IPrinter::printVector(14, 10, u.row(m));
        //break;
    }

    ems.clear();
    x.clear();
    b.clear();
    A.clear();
}

void IHeatEquationIBVP::calculateN4L2RDX(DoubleMatrix &u) const
{
    /* get parameters */
    double ht = _timeDimension.step();
    unsigned int minM = _timeDimension.min();
    unsigned int maxM = _timeDimension.max();
    unsigned int M = maxM - minM;

    double hx = _spaceDimensionX.step();
    unsigned int minN = _spaceDimensionX.min();
    unsigned int maxN = _spaceDimensionX.max();
    unsigned int N = maxN - minN;

    const unsigned int k = 4;
    double h = ht/(24.0*hx*hx);

    /*****************************************/

    u.clear();
    u.resize(M+1, N+1);

    DoubleMatrix A(k, k, 0.0);
    DoubleVector b(k, 0.0);
    DoubleVector x(k, 0.0);
    DoubleMatrix ems(N-k, k);

    /* initial condition */
    SpaceNodePDE isn;
    for (unsigned int n=0; n<=N; n++)
    {
        isn.i = n+minN;
        isn.x = isn.i*hx;
        u[0][n] = initial(isn, InitialCondition::InitialValue);
    }

    TimeNodePDE tn;
    SpaceNodePDE lsn;
    SpaceNodePDE rsn;
    lsn.i = minN; lsn.x = minN*hx;
    rsn.i = maxN; rsn.x = maxN*hx;
    BoundaryConditionPDE condition;
    for (unsigned int m=1; m<=M; m++)
    {
        tn.i = m+minM;
        tn.t = tn.i*ht;

        /* border conditions */
        u[m][0] = boundary(lsn, tn, condition);
        u[m][N] = boundary(rsn, tn, condition);

        /* using 2nd scheme, at point n=1 Березин И.С., Жидков Н.П. - Методы вычислений (том 1) */
        isn.i = minN+1;
        isn.x = isn.i*hx;
        double alpha = thermalDiffusivity()*h;
        A[0][0] = -40.0*alpha - 1.0;
        A[0][1] = +12.0*alpha;
        A[0][2] = +8.0*alpha;
        A[0][3] = -2.0*alpha;
        b[0]    = -u[m-1][1] - (+22.0*alpha)*u[m][0] - ht*f(isn,tn);

        /* using 5th scheme, at point n=4 Березин И.С., Жидков Н.П. - Методы вычислений (том 1) */
        //isn.i = minN+4;
        //isn.x = isn.i*hx;
        //double alpha = a(isn,tn)*h;
        //A[0][0] = -112.0*alpha;
        //A[0][1] = +228.0*alpha;
        //A[0][2] = -208.0*alpha;
        //A[0][3] = +70.0*alpha - 1.0;
        //b[0]    = -u[m-1][4] - (+22.0*alpha)*u[m][0] - ht*f(isn,tn);

        //double x1 = (1)*hx;
        //double x2 = (2)*hx;
        //double x3 = (3)*hx;
        //double x4 = (4)*hx;
        //printf("0 %14.10f %14.10f\n", b[0], A[0][0]*x1*x1*tn.t+A[0][1]*x2*x2*tn.t+A[0][2]*x3*x3*tn.t+A[0][3]*x4*x4*tn.t);

        A[0][1] /= A[0][0];
        A[0][2] /= A[0][0];
        A[0][3] /= A[0][0];
        b[0]    /= A[0][0];
        A[0][0] = 1.0;

        ems[0][0] = A[0][1];
        ems[0][1] = A[0][2];
        ems[0][2] = A[0][3];
        ems[0][3] = b[0];

#ifdef SCHEME_1
        unsigned int start = 1;
        unsigned int end   = N-k-1;
#endif
#ifdef SCHEME_2
        unsigned int start = 2;
        unsigned int end   = N-k+0;
#endif
#ifdef SCHEME_3
        unsigned int start = 3;
        unsigned int end   = N-k+1;
#endif
#ifdef SCHEME_4
        unsigned int start = 4;
        unsigned int end   = N-k+2;
#endif
#ifdef SCHEME_5
        unsigned int start = 5;
        unsigned int end   = N-k+3;
#endif
        for (unsigned int n=start; n<=end; n++)
        {
            isn.i = n+minN;
            isn.x = isn.i*hx;

            alpha = thermalDiffusivity()*h;

#ifdef SCHEME_1
            /* using 1nd scheme, at point n=1 Березин И.С., Жидков Н.П. - Методы вычислений (том 1) */
            double g1 = +70.0*alpha-1.0;
            double g2 = -208.0*alpha;
            double g3 = +228.0*alpha;
            double g4 = -112.0*alpha;
            double g5 = +22.0*alpha;
            double fi = -u[m-1][n] - ht*f(isn,tn);
#endif
#ifdef SCHEME_2
            /* using 2nd scheme, at point n=2 Березин И.С., Жидков Н.П. - Методы вычислений (том 1) */
            double g1 = +22.0*alpha;
            double g2 = -40.0*alpha-1.0;
            double g3 = +12.0*alpha;
            double g4 = +8.0*alpha;
            double g5 = -2.0*alpha;
            double fi = -u[m-1][n] - ht*f(isn,tn);
#endif
#ifdef SCHEME_3
            /* using 3nd scheme, at point n=3 Березин И.С., Жидков Н.П. - Методы вычислений (том 1) */
            double g1 = -2.0*alpha;
            double g2 = +32.0*alpha;
            double g3 = -60.0*alpha-1.0;
            double g4 = +32.0*alpha;
            double g5 = -2.0*alpha;
            double fi = -u[m-1][n] - ht*f(isn,tn);
#endif
#ifdef SCHEME_4
            /* using 4nd scheme, at point n=4 Березин И.С., Жидков Н.П. - Методы вычислений (том 1) */
            double g1 = -2.0*alpha;
            double g2 = +8.0*alpha;
            double g3 = +12.0*alpha;
            double g4 = -40.0*alpha-1.0;
            double g5 = +22.0*alpha;
            double fi = -u[m-1][n] - ht*f(isn,tn);
#endif
#ifdef SCHEME_5
            /* using 5nd scheme, at point n=5 Березин И.С., Жидков Н.П. - Методы вычислений (том 1) */
            double g1 = +22.0*alpha;
            double g2 = -112.0*alpha;
            double g3 = +228.0*alpha;
            double g4 = -208.0*alpha;
            double g5 = +70.0*alpha-1.0;
            double fi = -u[m-1][n] - ht*f(isn,tn);
#endif

            g2 /= -g1;
            g3 /= -g1;
            g4 /= -g1;
            g5 /= -g1;
            fi /= +g1;
            g1 = 1.0;

            A[0][0] = A[0][1] + g2;
            A[0][1] = A[0][2] + g3;
            A[0][2] = A[0][3] + g4;
            A[0][3] = g5;
            b[0]    = b[0] - fi;

            //double x1 = (n+1)*hx;
            //double x2 = (n+2)*hx;
            //double x3 = (n+3)*hx;
            //double x4 = (n+4)*hx;
            //printf("0 %14.10f %14.10f\n", b[0], A[0][0]*x1*x1*t+A[0][1]*x2*x2*t+A[0][2]*x3*x3*t+A[0][3]*x4*x4*t);
            \
            A[0][1] /= A[0][0];
            A[0][2] /= A[0][0];
            A[0][3] /= A[0][0];
            b[0]    /= A[0][0];
            A[0][0] = 1.0;

            // printf("0 %14.10f %14.10f\n", b[0], A[0][0]*x1*x1*t+A[0][1]*x2*x2*t+A[0][2]*x3*x3*t+A[0][3]*x4*x4*t);

#ifdef SCHEME_1
            ems[n+0][0] = A[0][1]; ems[n+0][1] = A[0][2]; ems[n+0][2] = A[0][3]; ems[n+0][3] = b[0];
#endif
#ifdef SCHEME_2
            ems[n-1][0] = A[0][1]; ems[n-1][1] = A[0][2]; ems[n-1][2] = A[0][3]; ems[n-1][3] = b[0];
#endif
#ifdef SCHEME_3
            ems[n-2][0] = A[0][1]; ems[n-2][1] = A[0][2]; ems[n-2][2] = A[0][3]; ems[n-2][3] = b[0];
#endif
#ifdef SCHEME_4
            ems[n-3][0] = A[0][1]; ems[n-3][1] = A[0][2]; ems[n-3][2] = A[0][3]; ems[n-3][3] = b[0];
#endif
#ifdef SCHEME_5
            ems[n-4][0] = A[0][1]; ems[n-4][1] = A[0][2]; ems[n-4][2] = A[0][3]; ems[n-4][3] = b[0];
#endif
        }

        /* using 2nd scheme, at point N-3 Березин И.С., Жидков Н.П. - Методы вычислений (том 1) */
        isn.i = maxN-3;
        isn.x = isn.i*hx;
        alpha = thermalDiffusivity()*h;
        A[1][0] = +22.0*alpha;
        A[1][1] = -40.0*alpha - 1.0;
        A[1][2] = +12.0*alpha;
        A[1][3] = +8.0*alpha;
        b[1]    = -u[m-1][N-3] - (-2.0*alpha)*u[m][N] - ht*f(isn,tn);

        /* using 3rd scheme, at point N-2 Березин И.С., Жидков Н.П. - Методы вычислений (том 1) */
        isn.i = maxN-2;
        isn.x = isn.i*hx;
        alpha = thermalDiffusivity()*h;
        A[2][0] = -2.0*alpha;
        A[2][1] = +32.0*alpha;
        A[2][2] = -60.0*alpha - 1.0;
        A[2][3] = +32.0*alpha;
        b[2]    = -u[m-1][N-2] - (-2.0*alpha)*u[m][N] - ht*f(isn,tn);

        /* using 4th scheme, at point N-1 Березин И.С., Жидков Н.П. - Методы вычислений (том 1) */
        isn.i = maxN-1;
        isn.x = isn.i*hx;
        A[3][0] = -2.0*alpha;
        A[3][1] = +8.0*alpha;
        A[3][2] = +12.0*alpha;
        A[3][3] = -40.0*alpha - 1.0;
        b[3]    = -u[m-1][N-1] - (+22.0*alpha)*u[m][N] - ht*f(isn,tn);

        //double xN4 = (N-4)*hx;
        //double xN3 = (N-3)*hx;
        //double xN2 = (N-2)*hx;
        //double xN1 = (N-1)*hx;

        //puts("------------------------");
        //printf("0 %14.10f %14.10f\n", b[0], A[0][0]*xN4*xN4*tn.t+A[0][1]*xN3*xN3*tn.t+A[0][2]*xN2*xN2*tn.t+A[0][3]*xN1*xN1*tn.t);
        //printf("1 %14.10f %14.10f\n", b[1], A[1][0]*xN4*xN4*tn.t+A[1][1]*xN3*xN3*tn.t+A[1][2]*xN2*xN2*tn.t+A[1][3]*xN1*xN1*tn.t);
        //printf("2 %14.10f %14.10f\n", b[2], A[2][0]*xN4*xN4*tn.t+A[2][1]*xN3*xN3*tn.t+A[2][2]*xN2*xN2*tn.t+A[2][3]*xN1*xN1*tn.t);
        //printf("3 %14.10f %14.10f\n", b[3], A[3][0]*xN4*xN4*tn.t+A[3][1]*xN3*xN3*tn.t+A[3][2]*xN2*xN2*tn.t+A[3][3]*xN1*xN1*tn.t);
        //puts("------------------------");

        LinearEquation::GaussianElimination(A, b, x);

        printf("x %14.10f %14.10f %14.10f %14.10f\n", x[0], x[1], x[2], x[3]);

        u[m][N-1] = x[3];
        u[m][N-2] = x[2];
        u[m][N-3] = x[1];
        u[m][N-4] = x[0];

#ifdef SCHEME_1
        for (unsigned int n=N-k-1; n>=1; n--) u[m][n+0] = -ems[n-1][0]*u[m][n+1] - ems[n-1][1]*u[m][n+2] - ems[n-1][2]*u[m][n+3] + ems[n-1][3];
#endif
#ifdef SCHEME_2
        for (unsigned int n=N-k+0; n>=2; n--) u[m][n-1] = -ems[n-2][0]*u[m][n+0] - ems[n-2][1]*u[m][n+1] - ems[n-2][2]*u[m][n+2] + ems[n-2][3];
#endif
#ifdef SCHEME_3
        for (unsigned int n=N-k+1; n>=3; n--) u[m][n-2] = -ems[n-3][0]*u[m][n-1] - ems[n-3][1]*u[m][n+0] - ems[n-3][2]*u[m][n+1] + ems[n-3][3];
#endif
#ifdef SCHEME_4
        for (unsigned int n=N-k+2; n>=4; n--) u[m][n-3] = -ems[n-4][0]*u[m][n-2] - ems[n-4][1]*u[m][n-1] - ems[n-4][2]*u[m][n+0] + ems[n-4][3];
#endif
#ifdef SCHEME_5
        for (unsigned int n=N-k+3; n>=5; n--) u[m][n-4] = -ems[n-5][0]*u[m][n-3] - ems[n-5][1]*u[m][n-2] - ems[n-5][2]*u[m][n-1] + ems[n-5][3];
#endif
        IPrinter::printVector(18, 10, u.row(m));
        break;
    }

    ems.clear();
    x.clear();
    b.clear();
    A.clear();
}

void IHeatEquationIBVP::calculateN6L2RD(DoubleMatrix &u) const
{
    /* get parameters */
    double ht = _timeDimension.step();
    unsigned int minM = _timeDimension.min();
    unsigned int maxM = _timeDimension.max();
    unsigned int M = maxM - minM;

    double hx = _spaceDimensionX.step();
    unsigned int minN = _spaceDimensionX.min();
    unsigned int maxN = _spaceDimensionX.max();
    unsigned int N = maxN - minN;

    const unsigned int k = 6;
    double h = ht/(24.0*hx*hx);

    /*****************************************/

    u.clear();
    u.resize(M+1, N+1);

    DoubleMatrix A(k, k, 0.0);
    DoubleVector b(k, 0.0);
    DoubleVector x(k, 0.0);
    DoubleMatrix ems(N-k, k);
    DoubleVector emk(N-k, k);

    /* initial condition */
    SpaceNodePDE isn;
    for (unsigned int n=0; n<=N; n++)
    {
        isn.i = n+minN;
        isn.x = isn.i*hx;
        u[0][n] = initial(isn, InitialCondition::InitialValue);
    }

    TimeNodePDE tn;
    SpaceNodePDE lsn;
    SpaceNodePDE rsn;
    lsn.i = minN; lsn.x = minN*hx;
    rsn.i = maxN; rsn.x = maxN*hx;
    BoundaryConditionPDE condition;
    for (unsigned int m=1; m<=M; m++)
    {
        tn.i = m+minM;
        tn.t = tn.i*ht;

        /* border conditions */
        u[m][0] = boundary(lsn, tn, condition);
        u[m][N] = boundary(rsn, tn, condition);

        /* using 2nd scheme, at point n=1 Березин И.С., Жидков Н.П. - Методы вычислений (том 1) */
        isn.i = minN+1;
        isn.x = isn.i*hx;
        double alpha = thermalDiffusivity()*h;
        A[0][0] = -40.0*alpha - 1.0;
        A[0][1] = +12.0*alpha;
        A[0][2] = +8.0*alpha;
        A[0][3] = -2.0*alpha;
        b[0]    = -u[m-1][1] - (+22.0*alpha)*u[m][0] - ht*f(isn,tn);

        /* using 5th scheme, at point n=4 Березин И.С., Жидков Н.П. - Методы вычислений (том 1) */
        //isn.i = minN+4;
        //isn.x = isn.i*hx;
        //double alpha = a(isn,tn)*h;
        //A[0][0] = -112.0*alpha;
        //A[0][1] = +228.0*alpha;
        //A[0][2] = -208.0*alpha;
        //A[0][3] = +70.0*alpha - 1.0;
        //b[0]    = -u[m-1][4] - (+22.0*alpha)*u[m][0] - ht*f(isn,tn);

        //double x1 = (1)*hx;
        //double x2 = (2)*hx;
        //double x3 = (3)*hx;
        //double x4 = (4)*hx;
        //printf("0 %14.10f %14.10f\n", b[0], A[0][0]*x1*x1*tn.t+A[0][1]*x2*x2*tn.t+A[0][2]*x3*x3*tn.t+A[0][3]*x4*x4*tn.t);

        //A[0][1] /= A[0][0];
        //A[0][2] /= A[0][0];
        //A[0][3] /= A[0][0];
        //b[0]    /= A[0][0];
        //A[0][0] = 1.0;

        ems[0][0] = A[0][1];
        ems[0][1] = A[0][2];
        ems[0][2] = A[0][3];
        ems[0][3] = b[0];
        emk[0] = A[0][0];

#define SCHEME_11

#ifdef SCHEME_1
        unsigned int start = 1;
        unsigned int end   = N-k-1;
#endif
#ifdef SCHEME_2
        unsigned int start = 2;
        unsigned int end   = N-k+0;
#endif
#ifdef SCHEME_3
        unsigned int start = 3;
        unsigned int end   = N-k+1;
#endif
#ifdef SCHEME_4
        unsigned int start = 4;
        unsigned int end   = N-k+2;
#endif
#ifdef SCHEME_5
        unsigned int start = 5;
        unsigned int end   = N-k+3;
#endif
        for (unsigned int n=start; n<=end; n++)
        {
            isn.i = n+minN;
            isn.x = isn.i*hx;

            alpha = thermalDiffusivity()*h;

#ifdef SCHEME_1
            /* using 1nd scheme, at point n=1 Березин И.С., Жидков Н.П. - Методы вычислений (том 1) */
            double g1 = +70.0*alpha-1.0;
            double g2 = -208.0*alpha;
            double g3 = +228.0*alpha;
            double g4 = -112.0*alpha;
            double g5 = +22.0*alpha;
            double fi = -u[m-1][n] - ht*f(isn,tn);
#endif
#ifdef SCHEME_2
            /* using 2nd scheme, at point n=2 Березин И.С., Жидков Н.П. - Методы вычислений (том 1) */
            double g1 = +22.0*alpha;
            double g2 = -40.0*alpha-1.0;
            double g3 = +12.0*alpha;
            double g4 = +8.0*alpha;
            double g5 = -2.0*alpha;
            double fi = -u[m-1][n] - ht*f(isn,tn);
#endif
#ifdef SCHEME_3
            /* using 3nd scheme, at point n=3 Березин И.С., Жидков Н.П. - Методы вычислений (том 1) */
            double g1 = -2.0*alpha;
            double g2 = +32.0*alpha;
            double g3 = -60.0*alpha-1.0;
            double g4 = +32.0*alpha;
            double g5 = -2.0*alpha;
            double fi = -u[m-1][n] - ht*f(isn,tn);
#endif
#ifdef SCHEME_4
            /* using 4nd scheme, at point n=4 Березин И.С., Жидков Н.П. - Методы вычислений (том 1) */
            double g1 = -2.0*alpha;
            double g2 = +8.0*alpha;
            double g3 = +12.0*alpha;
            double g4 = -40.0*alpha-1.0;
            double g5 = +22.0*alpha;
            double fi = -u[m-1][n] - ht*f(isn,tn);
#endif
#ifdef SCHEME_5
            /* using 5nd scheme, at point n=5 Березин И.С., Жидков Н.П. - Методы вычислений (том 1) */
            double g1 = +22.0*alpha;
            double g2 = -112.0*alpha;
            double g3 = +228.0*alpha;
            double g4 = -208.0*alpha;
            double g5 = +70.0*alpha-1.0;
            double fi = -u[m-1][n] - ht*f(isn,tn);
#endif

            g2 /= -g1;
            g3 /= -g1;
            g4 /= -g1;
            g5 /= -g1;
            fi /= +g1;
            g1 = 1.0;

            double a00 = A[0][0];
            A[0][0] = A[0][1] + g2*a00;
            A[0][1] = A[0][2] + g3*a00;
            A[0][2] = A[0][3] + g4*a00;
            A[0][3] = g5*a00;
            b[0]    = b[0] - fi*a00;

            //double x1 = (n+1)*hx;
            //double x2 = (n+2)*hx;
            //double x3 = (n+3)*hx;
            //double x4 = (n+4)*hx;
            //printf("0 %14.10f %14.10f\n", b[0], A[0][0]*x1*x1*t+A[0][1]*x2*x2*t+A[0][2]*x3*x3*t+A[0][3]*x4*x4*t);
            \
            //A[0][1] /= A[0][0];
            //A[0][2] /= A[0][0];
            //A[0][3] /= A[0][0];
            //b[0]    /= A[0][0];
            //A[0][0] = 1.0;

            // printf("0 %14.10f %14.10f\n", b[0], A[0][0]*x1*x1*t+A[0][1]*x2*x2*t+A[0][2]*x3*x3*t+A[0][3]*x4*x4*t);

#ifdef SCHEME_1
            ems[n+0][0] = A[0][1]; ems[n+0][1] = A[0][2]; ems[n+0][2] = A[0][3]; ems[n+0][3] = b[0]; emk[n+0] = A[0][0];
#endif
#ifdef SCHEME_2
            ems[n-1][0] = A[0][1]; ems[n-1][1] = A[0][2]; ems[n-1][2] = A[0][3]; ems[n-1][3] = b[0];
#endif
#ifdef SCHEME_3
            ems[n-2][0] = A[0][1]; ems[n-2][1] = A[0][2]; ems[n-2][2] = A[0][3]; ems[n-2][3] = b[0];
#endif
#ifdef SCHEME_4
            ems[n-3][0] = A[0][1]; ems[n-3][1] = A[0][2]; ems[n-3][2] = A[0][3]; ems[n-3][3] = b[0];
#endif
#ifdef SCHEME_5
            ems[n-4][0] = A[0][1]; ems[n-4][1] = A[0][2]; ems[n-4][2] = A[0][3]; ems[n-4][3] = b[0];
#endif
        }

        /* using 2nd scheme, at point N-3 Березин И.С., Жидков Н.П. - Методы вычислений (том 1) */
        isn.i = maxN-3;
        isn.x = isn.i*hx;
        alpha = thermalDiffusivity()*h;
        A[1][0] = +22.0*alpha;
        A[1][1] = -40.0*alpha - 1.0;
        A[1][2] = +12.0*alpha;
        A[1][3] = +8.0*alpha;
        b[1]    = -u[m-1][N-3] - (-2.0*alpha)*u[m][N] - ht*f(isn,tn);

        /* using 3rd scheme, at point N-2 Березин И.С., Жидков Н.П. - Методы вычислений (том 1) */
        isn.i = maxN-2;
        isn.x = isn.i*hx;
        alpha = thermalDiffusivity()*h;
        A[2][0] = -2.0*alpha;
        A[2][1] = +32.0*alpha;
        A[2][2] = -60.0*alpha - 1.0;
        A[2][3] = +32.0*alpha;
        b[2]    = -u[m-1][N-2] - (-2.0*alpha)*u[m][N] - ht*f(isn,tn);

        /* using 4th scheme, at point N-1 Березин И.С., Жидков Н.П. - Методы вычислений (том 1) */
        isn.i = maxN-1;
        isn.x = isn.i*hx;
        A[3][0] = -2.0*alpha;
        A[3][1] = +8.0*alpha;
        A[3][2] = +12.0*alpha;
        A[3][3] = -40.0*alpha - 1.0;
        b[3]    = -u[m-1][N-1] - (+22.0*alpha)*u[m][N] - ht*f(isn,tn);

        //double xN4 = (N-4)*hx;
        //double xN3 = (N-3)*hx;
        //double xN2 = (N-2)*hx;
        //double xN1 = (N-1)*hx;

        //puts("------------------------");
        //printf("0 %14.10f %14.10f\n", b[0], A[0][0]*xN4*xN4*tn.t+A[0][1]*xN3*xN3*tn.t+A[0][2]*xN2*xN2*tn.t+A[0][3]*xN1*xN1*tn.t);
        //printf("1 %14.10f %14.10f\n", b[1], A[1][0]*xN4*xN4*tn.t+A[1][1]*xN3*xN3*tn.t+A[1][2]*xN2*xN2*tn.t+A[1][3]*xN1*xN1*tn.t);
        //printf("2 %14.10f %14.10f\n", b[2], A[2][0]*xN4*xN4*tn.t+A[2][1]*xN3*xN3*tn.t+A[2][2]*xN2*xN2*tn.t+A[2][3]*xN1*xN1*tn.t);
        //printf("3 %14.10f %14.10f\n", b[3], A[3][0]*xN4*xN4*tn.t+A[3][1]*xN3*xN3*tn.t+A[3][2]*xN2*xN2*tn.t+A[3][3]*xN1*xN1*tn.t);
        //puts("------------------------");

        LinearEquation::GaussianElimination(A, b, x);

        //printf("x %14.10f %14.10f %14.10f %14.10f\n", x[0], x[1], x[2], x[3]);

        u[m][N-1] = x[3];
        u[m][N-2] = x[2];
        u[m][N-3] = x[1];
        u[m][N-4] = x[0];

#ifdef SCHEME_1
        for (unsigned int n=N-k-1; n>=1; n--) u[m][n+0] = (-ems[n-1][0]*u[m][n+1] - ems[n-1][1]*u[m][n+2] - ems[n-1][2]*u[m][n+3] + ems[n-1][3])/emk[n-1];
#endif
#ifdef SCHEME_2
        for (unsigned int n=N-k+0; n>=2; n--) u[m][n-1] = -ems[n-2][0]*u[m][n+0] - ems[n-2][1]*u[m][n+1] - ems[n-2][2]*u[m][n+2] + ems[n-2][3];
#endif
#ifdef SCHEME_3
        for (unsigned int n=N-k+1; n>=3; n--) u[m][n-2] = -ems[n-3][0]*u[m][n-1] - ems[n-3][1]*u[m][n+0] - ems[n-3][2]*u[m][n+1] + ems[n-3][3];
#endif
#ifdef SCHEME_4
        for (unsigned int n=N-k+2; n>=4; n--) u[m][n-3] = -ems[n-4][0]*u[m][n-2] - ems[n-4][1]*u[m][n-1] - ems[n-4][2]*u[m][n+0] + ems[n-4][3];
#endif
#ifdef SCHEME_5
        for (unsigned int n=N-k+3; n>=5; n--) u[m][n-4] = -ems[n-5][0]*u[m][n-3] - ems[n-5][1]*u[m][n-2] - ems[n-5][2]*u[m][n-1] + ems[n-5][3];
#endif
        //IPrinter::printVector(18, 10, u.row(m));
        //break;
    }

    ems.clear();
    x.clear();
    b.clear();
    A.clear();
}

void NewtonHeatEquation::calculateGM1(DoubleVector &u, SweepMethodDirection direction)
{
    typedef void (*t_algorithm)(const double*, const double*, const double*, const double*, double*, unsigned int);
    t_algorithm algorithm = &tomasAlgorithm;
    if (direction == ForwardSweep) algorithm = &tomasAlgorithmL2R;
    if (direction == BackwardSweep) algorithm = &tomasAlgorithmR2L;

    double ht = _timeDimension.step();
    unsigned int minM = _timeDimension.min();
    unsigned int maxM = _timeDimension.max();
    unsigned int M = maxM-minM;

    double hx = _spaceDimensionX.step();
    unsigned int minN = _spaceDimensionX.min();
    unsigned int maxN = _spaceDimensionX.max();
    unsigned int N = maxN-minN;

    u.clear();
    u.resize(N+1);

    double *ka = (double*) malloc(sizeof(double)*(N+1));
    double *kb = (double*) malloc(sizeof(double)*(N+1));
    double *kc = (double*) malloc(sizeof(double)*(N+1));
    double *kd = (double*) malloc(sizeof(double)*(N+1));
    double *rx = (double*) malloc(sizeof(double)*(N+1));

    /* initial condition */
    SpaceNodePDE isn;
    for (unsigned int n=0; n<=N; n++)
    {
        isn.i = n+minN;
        isn.x = isn.i*hx;
        u[n] = initial(isn, InitialCondition::InitialValue);
    }
    //layerInfo(u, 0);

    TimeNodePDE tn;
    for (unsigned int m=1; m<=M; m++)
    {
        tn.i = m+minM;
        tn.t = tn.i*ht;

        double aa = 0.0;

        isn.i = 0+minN;
        isn.x = isn.i*hx;
        aa = a(isn,tn);

        ka[0] = 0.0;
        kb[0] = 1.0+2.0*aa*(ht/(hx*hx))+2.0*lambda1*aa*(ht/hx) + lambda0*ht;
        kc[0] = -2.0*aa*(ht/(hx*hx));
        kd[0] = u[0] + lambda0*ht*theta0(tn) + 2.0*lambda1*aa*(ht/hx)*theta1(tn) + ht*f(isn,tn);
        for (unsigned int n=1; n<=N-1; n++)
        {
            isn.i = n+minN;
            isn.x = isn.i*hx;
            aa = a(isn,tn);

            ka[n] = -aa*(ht/(hx*hx));
            kb[n] = 1.0 + 2.0*aa*(ht/(hx*hx)) + lambda0*ht;
            kc[n] = -aa*(ht/(hx*hx));
            kd[n] = u[n] + lambda0*ht*theta0(tn) + ht*f(isn,tn);
        }
        isn.i = N+minN;
        isn.x = isn.i*hx;
        aa = a(isn,tn);

        ka[N] = -2.0*aa*(ht/(hx*hx));
        kb[N] = 1.0+2.0*aa*(ht/(hx*hx))+2.0*lambda2*aa*(ht/hx) + lambda0*ht;
        kc[N] = 0.0;
        kd[N] = u[N] + lambda0*ht*theta0(tn) + 2.0*lambda2*aa*(ht/hx)*theta2(tn) + ht*f(isn,tn);

        (*algorithm)(ka, kb, kc, kd, rx, N+1);

        for (unsigned int n=0; n<=N; n++) u[n] = rx[n];

        //layerInfo(u, m);
        break;
    }

    free(ka);
    free(kb);
    free(kc);
    free(kd);
    free(rx);
}

void NewtonHeatEquation::calculateGM2(DoubleVector &u, SweepMethodDirection direction)
{
    C_UNUSED(direction);
    //typedef void (*t_algorithm)(const double*, const double*, const double*, const double*, double*, unsigned int);
    //t_algorithm algorithm = &tomasAlgorithm;
    //if (direction == ForwardSweep) algorithm = &tomasAlgorithmL2R;
    //if (direction == BackwardSweep) algorithm = &tomasAlgorithmR2L;

    double ht = _timeDimension.step();
    unsigned int minM = _timeDimension.min();
    unsigned int maxM = _timeDimension.max();
    unsigned int M = maxM-minM;

    double hx = _spaceDimensionX.step();
    unsigned int minN = _spaceDimensionX.min();
    unsigned int maxN = _spaceDimensionX.max();
    unsigned int N = maxN-minN;

    u.clear();
    u.resize(N+1);

    double *ka = (double*) malloc(sizeof(double)*(N+1));
    double *kb = (double*) malloc(sizeof(double)*(N+1));
    double *kc = (double*) malloc(sizeof(double)*(N+1));
    double *kd = (double*) malloc(sizeof(double)*(N+1));
    double *rx = (double*) malloc(sizeof(double)*(N+1));

    /* initial condition */
    SpaceNodePDE isn;
    for (unsigned int n=0; n<=N; n++)
    {
        isn.i = n+minN;
        isn.x = isn.i*hx;
        u[n] = initial(isn, InitialCondition::InitialValue);
    }
    //layerInfo(u, 0);

    TimeNodePDE tn;
    for (unsigned int m=1; m<=M; m++)
    {
        tn.i = m+minM;
        tn.t = tn.i*ht;

        double aa = 0.0;

        isn.i = 0+minN;
        isn.x = isn.i*hx;
        aa = a(isn,tn);

        ka[0] = 0.0;
        kb[0] = 1.0+2.0*aa*(ht/(hx*hx))+2.0*lambda1*aa*(ht/hx) + lambda0*ht;
        kc[0] = -2.0*aa*(ht/(hx*hx));
        kd[0] = u[0] + lambda0*ht*theta0(tn) + 2.0*lambda1*aa*(ht/hx)*theta1(tn) + ht*f(isn,tn);
        for (unsigned int n=1; n<=N-1; n++)
        {
            isn.i = n+minN;
            isn.x = isn.i*hx;
            aa = a(isn,tn);

            ka[n] = -aa*(ht/(hx*hx));
            kb[n] = 1.0 + 2.0*aa*(ht/(hx*hx)) + lambda0*ht;
            kc[n] = -aa*(ht/(hx*hx));
            kd[n] = u[n] + lambda0*ht*theta0(tn) + ht*f(isn,tn);
        }
        isn.i = N+minN;
        isn.x = isn.i*hx;
        aa = a(isn,tn);

        ka[N] = -2.0*aa*(ht/(hx*hx));
        kb[N] = 1.0+2.0*aa*(ht/(hx*hx))+2.0*lambda2*aa*(ht/hx) + lambda0*ht;
        kc[N] = 0.0;
        kd[N] = u[N] + lambda0*ht*theta0(tn) + 2.0* lambda2*aa*(ht/hx)*theta2(tn) + ht*f(isn,tn);

        //(*algorithm)(ka, kb, kc, kd, rx, N+1);

        double *betta = (double *)malloc(sizeof(double)*(N+1));
        for (unsigned int i=0; i<=N; i++) betta[i] = 0.0;
        double eta = kd[0];
        betta[0] = kb[0];
        betta[1] = kc[0];

        //IPrinter::printVector(betta,N+1);

        double betta_norm = 0.0;
        for (unsigned int i=0; i<=N; i++)
        {
            betta_norm += betta[i]*betta[i];
        }
        betta_norm += eta*eta;
        betta_norm = sqrt(betta_norm);
        //printf("norm: %f\n", betta_norm);

        for (unsigned int i=0; i<=N; i++)
        {
            betta[i] = betta[i]/betta_norm;
        }
        eta /= betta_norm;
        //.0IPrinter::printVector(betta,N+1);

        betta_norm = 0.0;
        for (unsigned int i=0; i<=N; i++)
        {
            betta_norm += betta[i]*betta[i];
        }
        betta_norm = sqrt(betta_norm);
        //printf("norm: %f\n", betta_norm);

        for (unsigned int n=1; n<=N-1; n++)
        {
            betta[n+0] = -betta[n-1]*(kb[n]/ka[n]) + betta[n];
            betta[n+1] = -betta[n-1]*(kc[n]/ka[n]);
            eta = eta - betta[n-1]*(kd[n]/ka[n]);

            double betta_norm = 0.0;
            for (unsigned int i=0; i<=N; i++)
            {
                betta_norm += betta[i]*betta[i];
            }
            betta_norm += eta*eta;
            betta_norm = sqrt(betta_norm);
            for (unsigned int i=0; i<=N; i++)
            {
                betta[i] = betta[i]/betta_norm;
            }
            eta /= betta_norm;

        }

        DoubleMatrix M(2,2);
        DoubleVector A(2);
        DoubleVector x(2);

        M[0][0] = betta[N-1]; M[0][1] = betta[N]; A[0] = eta;
        M[1][0] = ka[N];      M[1][1] = kb[N];    A[1] = kd[N];

        LinearEquation::GaussianElimination(M,A,x);

        //for (unsigned int n=0; n<=N; n++) u[n] = rx[n];
        u[N-0] = x[1];
        u[N-1] = x[0];
        for (unsigned int n=N-1; n!=1; n--)
        {
            //betta[n+1] = +betta[n-1]*(kc[n]/ka[n]);
            betta[n] = betta[n] + betta[n-1]*(kb[n]/ka[n]);
            eta      = eta      + betta[n-1]*(kd[n]/ka[n]);
            u[n-1] = (eta-betta[n]*u[n])/betta[n-1];
        }

        //IPrinter::printVector(18,10,u);

        //layerInfo(u, m);
        break;
    }

    free(ka);
    free(kb);
    free(kc);
    free(kd);
    free(rx);
}

void NewtonHeatEquation::calculateGM3(DoubleVector &u, SweepMethodDirection direction)
{
    typedef void (*t_algorithm)(const double*, const double*, const double*, const double*, double*, unsigned int);
    t_algorithm algorithm = &tomasAlgorithm;
    if (direction == ForwardSweep) algorithm = &tomasAlgorithmL2R;
    if (direction == BackwardSweep) algorithm = &tomasAlgorithmR2L;

    double ht = _timeDimension.step();
    unsigned int minM = _timeDimension.min();
    unsigned int maxM = _timeDimension.max();
    unsigned int M = maxM-minM;

    double hx = _spaceDimensionX.step();
    unsigned int minN = _spaceDimensionX.min();
    unsigned int maxN = _spaceDimensionX.max();
    unsigned int N = maxN-minN;

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
        u[n] = initial(isn, InitialCondition::InitialValue);
    }
    //layerInfo(u, 0);

    TimeNodePDE tn;
    for (unsigned int m=1; m<=M; m++)
    {
        tn.i = m+minM;
        tn.t = tn.i*ht;

        double aa = 0.0;

        isn.i = 0+minN;
        isn.x = isn.i*hx;
        aa = a(isn,tn);
        double a0 = +1.0+lambda1*aa*(ht/hx) + lambda0*ht;
        double b0 = +aa*(ht/(hx*hx));
        double c0 = -aa*(ht/(hx*hx));
        double d0 = u[0] + lambda0*ht*theta0(tn) + lambda1*aa*(ht/hx)*theta1(tn) + ht*f(isn,tn);

        isn.i = 1+minN;
        isn.x = isn.i*hx;
        aa = a(isn,tn);
        double a1 = -aa*(ht/(hx*hx));
        double b1 = 1.0 + 2.0*aa*(ht/(hx*hx)) + lambda0*ht;
        double c1 = -aa*(ht/(hx*hx));
        double d1 = u[1] + lambda0*ht*theta0(tn) + ht*f(isn,tn);

        ka[0] = 0.0;
        kb[0] = b1-b0*(a1/a0);
        kc[0] = c1-c0*(a1/a0);
        kd[0] = d1-d0*(a1/a0);

        for (unsigned int n=2; n<=N-2; n++)
        {
            isn.i = n+minN;
            isn.x = isn.i*hx;
            aa = a(isn,tn);

            ka[n-1] = -aa*(ht/(hx*hx));
            kb[n-1] = 1.0 + 2.0*aa*(ht/(hx*hx)) + lambda0*ht;
            kc[n-1] = -aa*(ht/(hx*hx));
            kd[n-1] = u[n] + lambda0*ht*theta0(tn) + ht*f(isn,tn);
        }

        isn.i = (N-1)+minN;
        isn.x = isn.i*hx;
        aa = a(isn,tn);
        double aN1 = -aa*(ht/(hx*hx));
        double bN1 = 1.0 + 2.0*aa*(ht/(hx*hx)) + lambda0*ht;
        double cN1 = -aa*(ht/(hx*hx));
        double dN1 = u[N-1] + lambda0*ht*theta0(tn) + ht*f(isn,tn);

        isn.i = N+minN;
        isn.x = isn.i*hx;
        aa = a(isn,tn);
        double aN = -aa*(ht/(hx*hx));
        double bN = +aa*(ht/(hx*hx));
        double cN = +1.0+lambda2*aa*(ht/hx) + lambda0*ht;
        double dN = u[N] + lambda0*ht*theta0(tn) + lambda2*aa*(ht/hx)*theta2(tn) + ht*f(isn,tn);

        ka[N-2] = aN1-aN*(cN1/cN);
        kb[N-2] = bN1-bN*(cN1/cN);
        kc[N-2] = 0.0;
        kd[N-2] = dN1-dN*(cN1/cN);

        (*algorithm)(ka, kb, kc, kd, rx, N-1);

        for (unsigned int n=1; n<=N-1; n++)
            u[n] = rx[n-1];
        u[0] = -(b0/a0)*u[1]  -(c0/a0)*u[2]  +d0/a0;
        u[N] = -(bN/cN)*u[N-1]-(aN/cN)*u[N-2]+dN/cN;

        //layerInfo(u, m);
        break;
    }

    free(ka);
    free(kb);
    free(kc);
    free(kd);
    free(rx);
}

IHeatEquationFBVP::IHeatEquationFBVP(double thermalDiffusivity) : _thermalDiffusivity(thermalDiffusivity) {}

IHeatEquationFBVP::IHeatEquationFBVP(const IHeatEquationFBVP &ibvp)
{
    this->_timeDimension = ibvp._timeDimension;
    this->_spaceDimensionX = ibvp._spaceDimensionX;
    this->_spaceDimensionY = ibvp._spaceDimensionY;
    this->_spaceDimensionZ = ibvp._spaceDimensionZ;
    this->_thermalDiffusivity = ibvp._thermalDiffusivity;
}

IHeatEquationFBVP & IHeatEquationFBVP::operator =(const IHeatEquationFBVP &other)
{
    if (this == &other) { return *this; }

    this->_timeDimension = other._timeDimension;
    this->_spaceDimensionX = other._spaceDimensionX;
    this->_spaceDimensionY = other._spaceDimensionY;
    this->_spaceDimensionZ = other._spaceDimensionZ;
    this->_thermalDiffusivity = other._thermalDiffusivity;
    return *this;
}

IHeatEquationFBVP::~IHeatEquationFBVP() {}

double IHeatEquationFBVP::thermalDiffusivity() const { return _thermalDiffusivity; }

void IHeatEquationFBVP::setThermalDiffusivity(double thermalDiffusivity)
{
    this->_thermalDiffusivity = thermalDiffusivity;
}

double IHeatEquationFBVP::thermalConductivity() const { return _thermalConductivity; }

void IHeatEquationFBVP::setThermalConductivity(double convection)
{
    this->_thermalConductivity = convection;
}

double IHeatEquationFBVP::thermalConvection() const { return _thermalConvection; }

void IHeatEquationFBVP::setThermalConvection(double thermalConvection)
{
    this->_thermalConvection = thermalConvection;
}

void IHeatEquationFBVP::implicit_calculate_D1V1CN() const
{
    const unsigned int N = static_cast<unsigned int>(spaceDimensionX().size()) - 1;
    const unsigned int M = static_cast<unsigned int>(timeDimension().size()) - 1;

    const double hx = spaceDimensionX().step();
    const double ht = timeDimension().step();

    const double td = thermalDiffusivity();
    const double tc = thermalConductivity();
    const double cv = thermalConvection();
    const double w1 = weight();
    const double w2 = 1.0 - weight();

    const double p_td_ht__hxhx_w1 = +((td*ht)/(hx*hx))*w1 - ((cv*ht)/(2.0*hx))*w1;
    const double b_td_ht__hxhx_w1_tc_ht_w1 = +(1.0 - ((2.0*td*ht)/(hx*hx))*w1 + tc*ht*w1);

    const double ht_tc_w2 = +ht*tc*w2;
    const double m_td_ht__hxhx_w2 = -((td*ht)/(hx*hx))*w2;
    const double m_cv_ht__20hx_w2 = -((cv*ht)/(2.0*hx))*w2;
    const double p_td_ht__hx_w1 = -((td*ht)/hx)*w1;

    double *ax = static_cast<double*>(malloc(sizeof(double)*(N+1)));
    double *bx = static_cast<double*>(malloc(sizeof(double)*(N+1)));
    double *cx = static_cast<double*>(malloc(sizeof(double)*(N+1)));
    double *dx = static_cast<double*>(malloc(sizeof(double)*(N+1)));
    double *rx = static_cast<double*>(malloc(sizeof(double)*(N+1)));

    for (unsigned int n=0; n<=N; n++)
    {
        ax[n] = p_td_ht__hxhx_w1;
        bx[n] = b_td_ht__hxhx_w1_tc_ht_w1;
        cx[n] = p_td_ht__hxhx_w1;
    }
    ax[0] = 0.0; cx[N] = 0.0;

    DoubleVector u00(N+1);
    DoubleVector u10(N+1);

    /***********************************************************************************************/
    /***********************************************************************************************/

    TimeNodePDE tn0; tn0.i = M; tn0.t = tn0.i*ht;
    SpaceNodePDE sn;
    for (unsigned int n=0; n<=N; n++)
    {
        sn.i = static_cast<int>(n); sn.x = sn.i*hx;
        u00[n] = final(sn, FinalCondition::FinalValue);
    }
    layerInfo(u00, tn0);

    /***********************************************************************************************/
    /***********************************************************************************************/

    for (unsigned int ln=M-1, size_ln = static_cast<unsigned int>(0)-1; ln != size_ln; ln--)
    {
        TimeNodePDE tn00; tn00.i = ln+1; tn00.t = tn00.i*ht;
        TimeNodePDE tn10; tn10.i = ln;   tn10.t = tn10.i*ht;

        for (unsigned int n=1; n<=N-1; n++)
        {
            sn.i = static_cast<int>(n); sn.x = n*hx;
            dx[n] = 0.0;
            dx[n] += u00[n] * (1.0 - ht_tc_w2);
            dx[n] += (u00[n-1] - 2.0*u00[n] + u00[n+1])*m_td_ht__hxhx_w2;
            dx[n] += (u00[n+1] - u00[n-1])*m_cv_ht__20hx_w2;
            dx[n] -= ht*(w1*f(sn, tn10)+w2*f(sn, tn00));
        }

        unsigned int s=0, e=N;
        BoundaryConditionPDE condition; double alpha, beta, gamma, value;

        sn.i = static_cast<int>(0); sn.x = 0*hx;
        value = boundary(sn, tn10, condition);
        alpha = condition.alpha();
        beta  = condition.beta();
        gamma = condition.gamma();
        if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
        {
            s = 1;

            u10[0] = (gamma/alpha)*value;
            dx[1] -= u10[0]*p_td_ht__hxhx_w1;
            ax[1] = ax[0] = bx[0] = cx[0] = dx[0] = rx[0] = 0.0;
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
        {
            s = 0;

            ax[s]  = 0.0;
            bx[s]  = beta *(b_td_ht__hxhx_w1_tc_ht_w1);
            cx[s]  = beta *(2.0*p_td_ht__hxhx_w1);

            dx[s]  = u00[s]*(1.0 - ht_tc_w2);
            dx[s] += (2.0*u00[s]-5.0*u00[s+1]+4.0*u00[s+2]-u00[s+3])*m_td_ht__hxhx_w2;
            dx[s] -= ht*(w1*f(sn, tn10)+w2*f(sn, tn00));
            dx[s] *= beta;

            dx[s] += gamma*(-2.0*p_td_ht__hx_w1)*value;
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Robin)
        {
            s = 0;

            ax[s]  = 0.0;
            bx[s]  = beta *(b_td_ht__hxhx_w1_tc_ht_w1);
            cx[s]  = beta *(2.0*p_td_ht__hxhx_w1);

            dx[s]  = u00[s]*(1.0 - ht_tc_w2);
            dx[s] += (2.0*u00[s]-5.0*u00[s+1]+4.0*u00[s+2]-u00[s+3])*m_td_ht__hxhx_w2;
            dx[s] -= ht*(w1*f(sn, tn10)+w2*f(sn, tn00));
            dx[s] *= beta;

            bx[s] += alpha*(-2.0*p_td_ht__hx_w1);
            dx[s] += gamma*(-2.0*p_td_ht__hx_w1)*value;
        }

        sn.i = static_cast<int>(N); sn.x = N*hx;
        value = boundary(sn, tn10, condition);
        alpha = condition.alpha();
        beta  = condition.beta();
        gamma = condition.gamma();
        if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
        {
            e = N-1;
            u10[N] = (gamma/alpha)*value;
            dx[N-1] -= u10[N]*p_td_ht__hxhx_w1;
            cx[N-1] = ax[N] = bx[N] = cx[N] = dx[N] = rx[N] = 0.0;
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
        {
            e = N;

            ax[e]  = beta *(2.0*p_td_ht__hxhx_w1);
            bx[e]  = beta *(b_td_ht__hxhx_w1_tc_ht_w1);
            cx[e]  = 0.0;

            dx[e]  = u00[e]*(1.0 - ht_tc_w2);
            dx[e] += (-u00[e]+4.0*u00[e-1]-5.0*u00[e-2]+2.0*u00[e-3])*m_td_ht__hxhx_w2;
            dx[e] -= ht*(w1*f(sn, tn10)+w2*f(sn, tn00));
            dx[e] *= beta;

            dx[e] += gamma*(+2.0*p_td_ht__hx_w1)*value;
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Robin)
        {
            e = N;

            ax[e]  = beta *(2.0*p_td_ht__hxhx_w1);
            bx[e]  = beta *(b_td_ht__hxhx_w1_tc_ht_w1);
            cx[e]  = 0.0;

            dx[e]  = u00[e]*(1.0 - ht_tc_w2);
            dx[e] += (-u00[e]+4.0*u00[e-1]-5.0*u00[e-2]+2.0*u00[e-3])*m_td_ht__hxhx_w2;
            dx[e] -= ht*(w1*f(sn, tn10)+w2*f(sn, tn00));
            dx[e] *= beta;

            bx[e] += alpha*(+2.0*p_td_ht__hx_w1);
            dx[e] += gamma*(+2.0*p_td_ht__hx_w1)*value;
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

void IHeatEquationFBVP::gridMethod(DoubleVector &p, SweepMethodDirection direction) const
{
    typedef void (*t_algorithm)(const double*, const double*, const double*, const double*, double*, unsigned int);
    t_algorithm algorithm = &tomasAlgorithm;
    if (direction == ForwardSweep) algorithm = &tomasAlgorithmL2R;
    if (direction == BackwardSweep) algorithm = &tomasAlgorithmR2L;

    double ht = _timeDimension.step();
    unsigned int minM = _timeDimension.min();
    unsigned int maxM = _timeDimension.max();
    unsigned int M = maxM-minM;

    double hx = _spaceDimensionX.step();
    unsigned int minN = _spaceDimensionX.min();
    unsigned int maxN = _spaceDimensionX.max();
    unsigned int N = maxN-minN;

    double h = ht/(hx*hx);

    p.clear();
    p.resize(N+1);

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
        p[n] = final(isn, FinalCondition::FinalValue);
    }
    //layerInfo(p, M);

    SpaceNodePDE lsn;
    lsn.i = minN;
    lsn.x = minN*hx;

    SpaceNodePDE rsn;
    rsn.i = maxN;
    rsn.x = maxN*hx;

    TimeNodePDE tn;
    for (unsigned int m=M-1; m!=UINT32_MAX; m--)
    {
        tn.i = m+minM;
        tn.t = tn.i*ht;

        for (unsigned int n=1; n<=N-1; n++)
        {
            isn.i = n+minN;
            isn.x = isn.i*hx;

            double alpha = -thermalDiffusivity()*h;
            double betta = 1.0 - 2.0*alpha;

            ka[n-1] = alpha;
            kb[n-1] = betta;
            kc[n-1] = alpha;
            kd[n-1] = p[n] - ht * f(isn, tn);
        }

        ka[0]   = 0.0;
        kc[N-2] = 0.0;

        /* border conditions */
        //p[0] = boundary(lsn, tn);
        //p[N] = boundary(rsn, tn);

        kd[0]   -= thermalDiffusivity() * h * p[0];
        kd[N-2] -= thermalDiffusivity() * h * p[N];

        (*algorithm)(ka, kb, kc, kd, rx, N-1);

        for (unsigned int n=1; n<=N-1; n++) p[n] = rx[n-1];

        //layerInfo(p, m);
    }

    free(ka);
    free(kb);
    free(kc);
    free(kd);
    free(rx);
}

void IHeatEquationFBVP::gridMethod(DoubleMatrix &p, SweepMethodDirection direction) const
{
    typedef void (*t_algorithm)(const double*, const double*, const double*, const double*, double*, unsigned int);
    t_algorithm algorithm = &tomasAlgorithm;
    if (direction == ForwardSweep) algorithm = &tomasAlgorithmL2R;
    if (direction == BackwardSweep) algorithm = &tomasAlgorithmR2L;

    double ht = _timeDimension.step();
    unsigned int minM = _timeDimension.min();
    unsigned int maxM = _timeDimension.max();
    unsigned int M = maxM-minM;

    double hx = _spaceDimensionX.step();
    unsigned int minN = _spaceDimensionX.min();
    unsigned int maxN = _spaceDimensionX.max();
    unsigned int N = maxN-minN;

    double h = ht/(hx*hx);

    p.clear();
    p.resize(M+1, N+1);

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
        p[M][n] = final(isn, FinalCondition::FinalValue);
    }
    //layerInfo(p, M);

    SpaceNodePDE lsn;
    lsn.i = minN;
    lsn.x = minN*hx;

    SpaceNodePDE rsn;
    rsn.i = maxN;
    rsn.x = maxN*hx;

    TimeNodePDE tn;
    for (unsigned int m=M-1; m!=UINT32_MAX; m--)
    {
        tn.i = m+minM;
        tn.t = tn.i*ht;

        for (unsigned int n=1; n<=N-1; n++)
        {
            isn.i = n+minN;
            isn.x = isn.i*hx;

            double alpha = -thermalDiffusivity()*h;
            double betta = 1.0 - 2.0*alpha;

            ka[n-1] = alpha;
            kb[n-1] = betta;
            kc[n-1] = alpha;
            kd[n-1] = p[m+1][n] - ht * f(isn, tn);
        }

        ka[0]   = 0.0;
        kc[N-2] = 0.0;

        /* border conditions */
        //p[m][0] = boundary(lsn, tn);
        //p[m][N] = boundary(rsn, tn);

        //kd[0]   += a(lsn,tn) * h * p[m][0];
        //kd[N-2] += a(rsn,tn) * h * p[m][N];

        (*algorithm)(ka, kb, kc, kd, rx, N-1);

        for (unsigned int n=1; n<=N-1; n++) p[m][n] = rx[n-1];

        //layerInfo(p, m);
    }

    free(ka);
    free(kb);
    free(kc);
    free(kd);
    free(rx);
}

double IHeatEquationIBVP::weight() const { return 0.5; }

double IHeatEquationFBVP::weight() const { return 0.5; }

