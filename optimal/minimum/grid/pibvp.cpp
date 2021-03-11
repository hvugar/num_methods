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

void IHeatEquationIBVP::explicit_calculate_D1V1() const
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

    if (ht >= 0.5*(hx*hx)) throw std::runtime_error("Differential scheme not steady!");

    // common parameters
    const double k1 = (a*ht)/(hx*hx);
    const double k2 = (b*ht)/(2.0*hx);
    const double k3 = 1.0 + c*ht;

    DoubleVector u0(N+1, 0.0);
    DoubleVector u1(N+1, 0.0);

    /***********************************************************************************************/
    /***********************************************************************************************/

    TimeNodePDE tn;
    SpaceNodePDE sn;
    unsigned int i=0;
    int n = 0;
    BoundaryConditionPDE condition; double value;

    for (n=xmin, sn.i=n, sn.x=n*hx, i=0; n<=xmax; ++n, sn.i=n, sn.x=n*hx, ++i)
    {
        u0[i] = initial(sn, InitialCondition::InitialValue);
    }

    tn.i = 0; tn.t = tn.i*ht;
    layerInfo(u0, tn);

    /***********************************************************************************************/

    for (unsigned int ln=1; ln<=L; ln++)
    {
        tn.i = ln; tn.t = tn.i*ht;

        for (n=xmin+1, sn.i=n, sn.x=n*hx, i=1; n<=xmax-1; ++n, sn.i=n, sn.x=i*hx, ++i)
        {
            u1[i] = k3*u0[i] + k1*(u0[i+1] - 2.0*u0[i] + u0[i-1]) + k2*(u0[i+1] - u0[i-1]) + ht*f(sn, tn);
        }

        sn.i = xmin; sn.x = xmin*hx;
        value = boundary(sn, tn, condition);

        if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
        {
            u1[0] = value;
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
        {
            u1[0] = u1[1] - hx*value;
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Robin)
        {
            const double alpha = condition.alpha();
            const double beta  = condition.beta();
            /*gamma = condition.gamma();*/
            const double gamma = 1.0;

            u1[0] = (beta*u1[1] - hx*gamma*value)/(beta-hx*alpha);
        }

        sn.i = xmax; sn.x = xmax*hx;
        value = boundary(sn, tn, condition);

        if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
        {
            u1[N] = /*(gamma/alpha)**/value;
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
        {
            u1[N] = u1[N-1] + hx*value;
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Robin)
        {
            const double alpha = condition.alpha();
            const double beta  = condition.beta();
            /*gamma = condition.gamma();*/
            const double gamma = 1.0;

            u1[N] = (beta*u1[N-1] + hx*gamma*value)/(beta+hx*alpha);
        }

        layerInfo(u1, tn);

        u0 = u1;
    }

    u0.clear();
    u1.clear();
}

void IHeatEquationIBVP::implicit_calculate_D1V1() const
{
    const size_t N = static_cast<size_t>(spaceDimensionX().size()) - 1;
    const size_t L = static_cast<size_t>(timeDimension().size()) - 1;

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
    const double k11 = -((a*ht)/(hx*hx))*w1 + ((b*ht)/(2.0*hx))*w1; // i-1
    const double k12 = +1.0 + ((2.0*a*ht)/(hx*hx))*w1 - c*ht*w1;    // i
    const double k13 = -((a*ht)/(hx*hx))*w1 - ((b*ht)/(2.0*hx))*w1; // i+1
    const double k21 = +((a*ht)/(hx*hx))*w2 - ((b*ht)/(2.0*hx))*w2; // i-1
    const double k22 = +1.0 - ((2.0*a*ht)/(hx*hx))*w2 + c*ht*w2;    // i
    const double k23 = +((a*ht)/(hx*hx))*w2 + ((b*ht)/(2.0*hx))*w2; // i+1
    
    // left border condition parameters
    const double b11 = +1.0 + ((2.0*a*ht)/(hx*hx))*w1 - c*ht*w1;
    const double b12 = -((2.0*a*ht)/hx)*w1 + b*ht*w1;
    const double b13 = -((2.0*a*ht)/(hx*hx))*w1;
    const double b14 = +1.0 - ((2.0*a*ht)/(hx*hx))*w2 + c*ht*w2;
    const double b15 = +((2.0*a*ht)/hx)*w2 - b*ht*w2;
    const double b16 = +((2.0*a*ht)/(hx*hx))*w2;
    const double b19 = -((2.0*a*ht)/hx) + b*ht;
    //const double b17 = -((2.0*a*ht)/hx)*w1 + b*ht*w1;
    //const double b18 = -((2.0*a*ht)/hx)*w2 + b*ht*w2;

    // right border condition parameters
    const double b21 = -((2.0*a*ht)/(hx*hx))*w1;
    const double b22 = +((2.0*a*ht)/hx)*w1 + b*ht*w1;
    const double b23 = +1.0 + ((2.0*a*ht)/(hx*hx))*w1 - c*ht*w1;
    const double b24 = +((2.0*a*ht)/(hx*hx))*w2;
    const double b25 = -((2.0*a*ht)/hx)*w2 - b*ht*w2;
    const double b26 = +1.0 - ((2.0*a*ht)/(hx*hx))*w2 + c*ht*w2;
    const double b29 = +((2.0*a*ht)/hx) + b*ht;
    //const double b27 = +((2.0*a*ht)/hx)*w1 + b*ht*w1;
    //const double b28 = +((2.0*a*ht)/hx)*w2 + b*ht*w2;

    // common parameters
    //const double ht_w1 = +ht*w1;
    //const double ht_w2 = +ht*w2;

    double *ax = static_cast<double*>(malloc(sizeof(double)*(N+1)));
    double *bx = static_cast<double*>(malloc(sizeof(double)*(N+1)));
    double *cx = static_cast<double*>(malloc(sizeof(double)*(N+1)));
    double *dx = static_cast<double*>(malloc(sizeof(double)*(N+1)));
    double *rx = static_cast<double*>(malloc(sizeof(double)*(N+1)));

    for (size_t n=0; n<=N; n++)
    {
        ax[n] = k11;
        bx[n] = k12;
        cx[n] = k13;
    }
    ax[0] = 0.0; cx[N] = 0.0;

    double *u0 = static_cast<double*>(malloc(sizeof(double)*(N+1)));
    double *u1 = static_cast<double*>(malloc(sizeof(double)*(N+1)));

    /***********************************************************************************************/
    /***********************************************************************************************/

    SpaceNodePDE sn;
    unsigned int i = 0;
    unsigned int k = static_cast<unsigned int>(timeDimension().min());
    for (int n=xmin; n<=xmax; n++, i++)
    {
        sn.i = n; sn.x = n*hx;
        u0[i] = initial(sn, InitialCondition::InitialValue);
    }
    TimeNodePDE tn; tn.i = k; tn.t = tn.i*ht; k++;
    layerInfo(DoubleVector(u0, N+1), tn);

    /***********************************************************************************************/
    /***********************************************************************************************/

    BoundaryConditionPDE condition; double value;

    for (size_t ln=1; ln<=L; ln++, k++)
    {
        TimeNodePDE tn0; tn0.i = static_cast<unsigned int>(k-1); tn0.t = tn0.i*ht;
        TimeNodePDE tn1; tn1.i = static_cast<unsigned int>(k  ); tn1.t = tn1.i*ht;

        size_t i = 1;
        for (int n=xmin+1; n<=xmax-1; n++, i++)
        {
            sn.i = n; sn.x = n*hx;

#ifdef PARABOLIC_IBVP_H_D1V1_FX
            double fx = f(sn, tn1);
#else
            double fx = w1*f(sn, tn1) + w2*f(sn, tn0);
#endif
            dx[i] = k21*u0[i-1] + k22*u0[i] + k23*u0[i+1] + ht*fx;
        }

        size_t s=0, e=N;

        sn.i = xmin; sn.x = xmin*hx;
#ifdef PARABOLIC_IBVP_H_D1V1_BR_LEFT
        value = boundary(sn, tn1, condition);
#else
        value = w1*boundary(sn, tn1, condition) + w2*boundary(sn, tn0, condition);
#endif

        if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
        {
            //            s = 1;
            //            u1[0]  = value;
            //            dx[1] -= k11 * u1[0];
            //            ax[1]  = ax[0] = bx[0] = cx[0] = dx[0] = rx[0] = 0.0;

            s = 0;
            ax[0] = 0.0;
            bx[0] = 1.0;
            cx[0] = 0.0;
            dx[0] = value;
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
        {
            s = 0;
#ifdef PARABOLIC_IBVP_H_D1V1_BORDER_O2
#ifdef PARABOLIC_IBVP_H_D1V1_FX
            double fx = f(sn, tn1);
#else
            double fx = w1*f(sn, tn1) + w2*f(sn, tn0);
#endif
            ax[s]  = 0.0;
            bx[s]  = b11;
            cx[s]  = b13;
            dx[s]  = b14*u0[s] + b16*u0[s+1];
            dx[s] += b19*value;
            dx[s] += ht*fx;
#else
            ax[s]  = 0.0;
            bx[s]  = -1.0;
            cx[s]  = +1.0;
            dx[s]  = ht*value;
#endif
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Robin)
        {
            const double alpha = condition.alpha();
            const double beta  = condition.beta();

            s = 0;
#ifdef PARABOLIC_IBVP_H_D1V1_BORDER_O2
#ifdef PARABOLIC_IBVP_H_D1V1_FX
            double fx = f(sn, tn1);
#else
            double fx = w1*f(sn, tn1) + w2*f(sn, tn0);
#endif
            ax[s]  = 0.0;
            bx[s]  = beta*b11 + alpha*b12;
            cx[s]  = beta*b13;
            dx[s]  = beta*b14*u0[s] + alpha*b15*u0[s] + beta*b16*u0[s+1];
            dx[s] += b19*value;
            dx[s] += beta*ht*fx;
#else
            ax[s]  = 0.0;
            bx[s]  = alpha*ht - beta;
            cx[s]  = beta;
            dx[s]  = ht*value;
#endif
        }

        sn.i = xmax; sn.x = xmax*hx;
#ifdef PARABOLIC_IBVP_H_D1V1_BR_RIGHT
        value = boundary(sn, tn1, condition);
#else
        value = w1*boundary(sn, tn1, condition) + w2*boundary(sn, tn0, condition);
#endif

        if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
        {
            //            e = N-1;
            //            u1[N] = value;
            //            dx[N-1] -= k13 * u1[N];
            //            cx[N-1] = ax[N] = bx[N] = cx[N] = dx[N] = rx[N] = 0.0;

            e = N;
            ax[N] = 0.0;
            bx[N] = 1.0;
            cx[N] = 0.0;
            dx[N] = value;
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
        {
            e = N;
#ifdef PARABOLIC_IBVP_H_D1V1_BORDER_O2
#ifdef PARABOLIC_IBVP_H_D1V1_FX
            double fx = f(sn, tn1);
#else
            double fx = w1*f(sn, tn1) + w2*f(sn, tn0);
#endif
            ax[e]  = b21;
            bx[e]  = b23;
            cx[e]  = 0.0;
            dx[e]  = b24 * u0[e-1] + b26 * u0[e];
            dx[e] += b29 * value;
            dx[e] += ht*fx;
#else
            ax[e]  = -1.0;
            bx[e]  = +1.0;
            cx[e]  = 0.0;
            dx[e]  = ht*value;
#endif
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Robin)
        {
            const double alpha = condition.alpha();
            const double beta  = condition.beta();

            e = N;
#ifdef PARABOLIC_IBVP_H_D1V1_BORDER_O2
#ifdef PARABOLIC_IBVP_H_D1V1_FX
            double fx = f(sn, tn1);
#else
            double fx = w1*f(sn, tn1) + w2*f(sn, tn0);
#endif
            ax[e]  = beta*b21;
            bx[e]  = alpha*b22 + beta*b23;
            cx[e]  = 0.0;
            dx[e]  = beta*b24*u0[e-1] + alpha*b25*u0[e] + beta*b26*u0[e];
            dx[e]  += b29*value;
            dx[e]  += beta*ht*fx;
#else
            ax[e]  = -beta;
            bx[e]  = alpha*ht + beta;
            cx[e]  = 0.0;
            dx[e]  = ht*value;
#endif
        }

        tomasAlgorithmLeft2Right(ax+s, bx+s, cx+s, dx+s, rx+s, static_cast<unsigned int>(e-s+1));
        memcpy(u1+s, rx+s, sizeof(double)*(e-s+1));
        //for (unsigned int n=s; n<=e; n++) u1[n] = rx[n];

        layerInfo(DoubleVector(u1, N+1), tn1);

        double *tmp_u = u0; u0 = u1; u1 = tmp_u;
    }

    free(u0);
    free(u1);

    free(rx);
    free(dx);
    free(cx);
    free(bx);
    free(ax);
}

void IHeatEquationIBVP::explicit_calculate_D2V1() const
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

    if (ht >= 0.5/(1.0/(hx*hx)+1.0/(hy*hy)))
    {
        std::string message = "Differential scheme is not steady!";
        fputs(message.data(), stderr);
        throw std::runtime_error(message);
    }

    // common parameters
    const double k1 = +((a1*ht)/(hx*hx));
    const double k2 = +((a1*ht)/(hy*hy));
    const double ht_tc = ht*b1;

    DoubleMatrix u00(M+1, N+1);
    DoubleMatrix u10(M+1, N+1);

    /***********************************************************************************************/
    /***********************************************************************************************/

    TimeNodePDE tn00, tn05, tn10;
    SpaceNodePDE sn;
    unsigned int i=0, j=0, s=0, e=N;
    int n = 0, m = 0;
    BoundaryConditionPDE condition; double value/*, alpha, beta, gamma*/;

    for (m=ymin, sn.j=m, sn.y=m*hy; m<=ymax; ++m, sn.j=m, sn.y=m*hy, ++j)
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
        tn00.i = ln; tn00.t = tn00.i*ht;

        /**************************************************** border conditions ***************************************************/

        SpaceNodePDE sn0, sn1;
        BoundaryConditionPDE condition;

        sn0.i = static_cast<int>(0); sn0.x = 0*hx;
        sn1.i = static_cast<int>(N); sn1.x = N*hx;
        for (unsigned int m=0; m<=M; m++)
        {
            sn0.j = sn1.j = static_cast<int>(m); sn0.y = sn1.y = m*hy;
            u10[m][0] = boundary(sn0, tn00, condition);
            u10[m][N] = boundary(sn1, tn00, condition);
        }

        sn0.j = static_cast<int>(0); sn0.y = 0*hy;
        sn1.j = static_cast<int>(M); sn1.y = M*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn0.i = sn1.i = static_cast<int>(n); sn0.x = sn1.x = n*hx;
            u10[0][n] = boundary(sn0, tn00, condition);
            u10[M][n] = boundary(sn1, tn00, condition);
        }

        /**************************************************** border conditions ***************************************************/

        sn0.j = ymin; sn0.y = ymin*hy;
        sn1.j = ymax; sn1.y = ymax*hy;
        for (n=xmin+1, sn.i=n, sn.x=n*hx, i=1; n<=xmax-1; ++n, sn.i=n, sn.x=i*hx, ++i)
        {
            sn0.i = n; sn0.x = n*hx;
            value = boundary(sn0, tn00, condition);

            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                u10[0][i] = value;
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
                u10[0][i] = (3.5*u10[1][i] - 2.0*u10[2][i] + 0.5*u10[3][i] + hy*value)/(2.0);
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                const double alpha = condition.alpha();
                const double beta  = condition.beta();
                /*gamma = condition.gamma();*/
                const double gamma = 1.0;

                u10[0][i] = (3.5*u10[1][i] - 2.0*u10[2][i] + 0.5*u10[3][i] + hy*(gamma/beta)*value)/(2.0 + (alpha/beta)*hy);
            }

            sn1.i = n; sn1.x = n*hx;
            value = boundary(sn1, tn00, condition);

            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                u10[M][i] = value;
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
                u10[M][i] = (3.5*u10[M-1][i] - 2.0*u10[M-2][i] + 0.5*u10[M-3][i] + hy*value)/(2.0);
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                const double alpha = condition.alpha();
                const double beta  = condition.beta();
                //gamma = condition.gamma();
                const double gamma = 1.0;

                u10[M][i] = (3.5*u10[M-1][i] - 2.0*u10[M-2][i] + 0.5*u10[M-3][i] + hy*(gamma/beta)*value)/(2.0 + (alpha/beta)*hy);
            }
        }

        /**************************************************************************************************************************/
        for (m=ymin+1, sn.j=m, sn.y=m*hy, j=1; m<=ymax-1; ++m, sn.j=m, sn.y=m*hy, ++j)
        {
            sn.i = xmin; sn.x = xmin*hx;
            value = boundary(sn, tn05, condition);

            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                u10[j][0] = value;
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                //const double alpha = condition.alpha();
                //const double beta  = condition.beta();
                /*gamma = condition.gamma();*/
                //const double gamma = 1.0;
            }

            for (n=xmin+1, sn.i=m, sn.x=n*hx, i=1; n<=xmax-1; ++n, sn.i=n, sn.x=i*hx, ++i)
            {
                u10[j][i] = u00[j][i] + k1*(u00[j][i-1] - 2.0*u00[j][i] + u00[j][i+1])
                        + k2*(u00[j-1][i] - 2.0*u00[j][i] + u00[j+1][i]) - ht_tc*u00[j][i] + ht*f(sn, tn00);
            }

            sn.i = xmax; sn.x = xmax*hx;
            value = boundary(sn, tn05, condition);

            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                u10[j][N] = /*(gamma/alpha)**/value;
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                //const double alpha = condition.alpha();
                //const double beta  = condition.beta();
                /*gamma = condition.gamma();*/
                //const double gamma = 1.0;

            }
        }
        layerInfo(u10, tn00);
        /**************************************************************************************************************************/

        u00 = u10;
    }

    u00.clear();
    u10.clear();
}

void IHeatEquationIBVP::implicit_calculate_D2V1() const
{
    const size_t N = static_cast<size_t>(spaceDimensionX().size()) - 1;
    const size_t M = static_cast<size_t>(spaceDimensionY().size()) - 1;
    const size_t L = static_cast<size_t>(timeDimension().size()) - 1;

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

    //double l1 = hx*fabs(xmax-xmin);
    //double l2 = hy*fabs(ymax-ymin);
    //double ht_max = (1.0/M_PI) * sqrt((hx*hx)/fabs(a1)+(hy*hy)/fabs(a2)) * (1.0/sqrt( (1.0/(l1*l1)) + (1.0/(l2*l2)) ));
    //if (ht > ht_max) { throw std::runtime_error("Differential scheme is conditionally steady."); }

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

    // left border condition parameters x = xmin
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

    for (size_t n=0; n<=N; n++)
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

    for (size_t m=0; m<=M; m++)
    {
        ay[m] = k21;
        by[m] = k22;
        cy[m] = k23;
    }
    ay[0] = 0.0; cy[M] = 0.0;

    double **u00 = static_cast<double**>(malloc(sizeof(double*)*(M+1)));
    double **u05 = static_cast<double**>(malloc(sizeof(double*)*(M+1)));
    double **u10 = static_cast<double**>(malloc(sizeof(double*)*(M+1)));
    for (size_t i=0; i<=M; i++)
    {
        u00[i] = static_cast<double*>(malloc(sizeof(double)*(N+1)));
        u05[i] = static_cast<double*>(malloc(sizeof(double)*(N+1)));
        u10[i] = static_cast<double*>(malloc(sizeof(double)*(N+1)));
    }

    /***********************************************************************************************/
    /***********************************************************************************************/

    TimeNodePDE tn00, tn10;
    SpaceNodePDE sn;
    size_t i=0, j=0, s=0, e=N;
    int n = 0, m = 0;
    unsigned int k = static_cast<unsigned int>(timeDimension().min());

    for (m=ymin, sn.j=m, sn.y=m*hy, j=0; m<=ymax; ++m, sn.j=m, sn.y=m*hy, ++j)
    {
        for (n=xmin, sn.i=n, sn.x=n*hx, i=0; n<=xmax; ++n, sn.i=n, sn.x=n*hx, ++i)
        {
            u00[j][i] = initial(sn, InitialCondition::InitialValue);
        }
    }

    tn00.i = k; tn00.t = tn00.i*ht; k++;
    layerInfo(DoubleMatrix(u00, M+1, N+1), tn00);

    /***********************************************************************************************/
    /***********************************************************************************************/

    BoundaryConditionPDE condition; double value;

    for (size_t ln=1; ln<=L; ln++, k++)
    {
        tn00.i = static_cast<unsigned int>(k-1); tn00.t = tn00.i*ht;
        tn10.i = static_cast<unsigned int>(k  ); tn10.t = tn10.i*ht;

        /**************************************************** x direction apprx ***************************************************/

        for (m=ymin+1, sn.j=m, sn.y=m*hy, j=1; m<=ymax-1; ++m, sn.j=m, sn.y=m*hy, ++j)
        {
            for (n=xmin+1, sn.i=n, sn.x=n*hx, i=1; n<=xmax-1; ++n, sn.i=n, sn.x=n*hx, ++i)
            {
#ifdef PARABOLIC_IBVP_H_D2V1_FX_X
                double fx = f(sn, tn00);
#else
                double fx = 0.5 * ( f(sn, tn00) + f(sn, tn10) );
#endif
                dx[i] = k14*u00[j-1][i] + k15*u00[j][i] + k16*u00[j+1][i] + ht05*fx;
            }

            sn.i = xmin; sn.x = xmin*hx;
#ifdef PARABOLIC_IBVP_H_D2V1_BR_X
            value = boundary(sn, tn00, condition);
#else
            value = 0.5 * ( boundary(sn, tn00, condition) + boundary(sn, tn10, condition) );
#endif

            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                s = 1;
                u05[j][0] = value;
                dx[1] -= k11 * u05[j][0];
                ax[1] = ax[0] = bx[0] = cx[0] = dx[0] = rx[0] = 0.0;
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
                s = 0;
#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
#ifdef PARABOLIC_IBVP_H_D2V1_FX_X
                double fx = f(sn, tn00);
#else
                double fx = 0.5 * ( f(sn, tn00) + f(sn, tn10) );
#endif
                ax[0]  = 0.0;
                bx[0]  = b111;
                cx[0]  = b113;
                dx[0]  = b114 * u00[j][0];
                dx[0] += b115 * ((u00[j+1][0]-2.0*u00[j][0]+u00[j-1][0])/(hy*hy));
                dx[0] += b116 * ((u00[j+1][0]-u00[j-1][0])/(2.0*hy));
                dx[0] += ht05 * fx;
                dx[0] += b117 * value;
#else
                ax[0]  = +0.0;
                bx[0]  = +1.0;
                cx[0]  = -1.0;
                dx[0]  = +hx * value;
#endif
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                const double alpha = condition.alpha();
                const double beta  = condition.beta();

                s = 0;
#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
#ifdef PARABOLIC_IBVP_H_D2V1_FX_X
                double fx = f(sn, tn00);
#else
                double fx = 0.5 * ( f(sn, tn00) + f(sn, tn10) );
#endif
                ax[0]  = 0.0;
                bx[0]  = beta * b111 + alpha * b112;
                cx[0]  = beta * b113;
                dx[0]  = beta * b114 * u00[j][0];
                dx[0] += beta * b115 * ((u00[j+1][0]-2.0*u00[j][0]+u00[j-1][0])/(hy*hy));
                dx[0] += beta * b116 * ((u00[j+1][0]-u00[j-1][0])/(2.0*hy));
                dx[0] += beta * ht05 * fx;
                dx[0] += b117 * value;
#else
                ax[0]  = +0.0;
                bx[0]  = +beta + alpha * hx;
                cx[0]  = -beta;
                dx[0]  = +hx * value;
#endif
            }

            sn.i = xmax; sn.x = xmax*hx;
#ifdef PARABOLIC_IBVP_H_D2V1_BR_X
            value = boundary(sn, tn00, condition);
#else
            value = 0.5 * ( boundary(sn, tn00, condition) + boundary(sn, tn10, condition) );
#endif
            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                e = N-1;
                u05[j][N] = value;
                dx[e] -= k13 * u05[j][N];
                cx[e] = ax[N] = bx[N] = cx[N] = dx[N] = rx[N] = 0.0;
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
                e = N;
#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
#ifdef PARABOLIC_IBVP_H_D2V1_FX_X
                double fx = f(sn, tn00);
#else
                double fx = 0.5 * ( f(sn, tn00) + f(sn, tn10) );
#endif
                ax[e]  = b121;
                bx[e]  = b123;
                cx[e]  = 0.0;
                dx[e]  = b124 * u00[j][e];
                dx[e] += b125 * ((u00[j+1][e]-2.0*u00[j][e]+u00[j-1][e])/(hy*hy));
                dx[e] += b126 * ((u00[j+1][e]-u00[j-1][e])/(2.0*hy));
                dx[e] += ht05 * fx;
                dx[e] += b127 * value;
#else
                ax[e]  = -1.0;
                bx[e]  = +1.0;
                cx[e]  = +0.0;
                dx[e]  = hx * value;
#endif
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                const double alpha = condition.alpha();
                const double beta  = condition.beta();

                e = N;
#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
#ifdef PARABOLIC_IBVP_H_D2V1_FX_X
                double fx = f(sn, tn00);
#else
                double fx = 0.5 * ( f(sn, tn00) + f(sn, tn10) );
#endif
                ax[e]  = beta * b121;
                bx[e]  = beta * b123 + alpha * b122;
                cx[e]  = 0.0;
                dx[e]  = beta * b124 * u00[j][e];
                dx[e] += beta * b125 * ((u00[j+1][e]-2.0*u00[j][e]+u00[j-1][e])/(hy*hy));
                dx[e] += beta * b126 * ((u00[j+1][e]-u00[j-1][e])/(2.0*hy));
                dx[e] += beta * ht05 * fx;
                dx[e] += b127 * value;
#else
                ax[e]  = -beta;
                bx[e]  = +beta + alpha * hx;
                cx[e]  = +0.0;
                dx[e]  = +hx * value;
#endif
            }

            tomasAlgorithmLeft2Right(ax+s, bx+s, cx+s, dx+s, rx+s, static_cast<unsigned int>(e-s+1));
            //for (size_t i=s; i<=e; i++) u05[j][i] = rx[i];
            memcpy(u05[j]+s, rx+s, sizeof (double*)*(e-s+1));
        }

        for (n=xmin, sn.i=n, sn.x=n*hx, i=0; n<=xmax; ++n, sn.i=n, sn.x=n*hx, ++i)
        {
            sn.j = ymin; sn.y = ymin*hy;
#ifdef PARABOLIC_IBVP_H_D2V1_BR_X
            value = boundary(sn, tn00, condition);
#else
            value = 0.5 * ( boundary(sn, tn00, condition) + boundary(sn, tn10, condition) );
#endif
            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                u05[0][i] = value;
            }
            if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
                u05[0][i] = (3.5*u05[1][i] - 2.0*u05[2][i] + 0.5*u05[3][i] + hy*value)/(2.0);
#else
                u05[0][i] = u05[1][i] + hy*value;
#endif
            }
            if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                const double alpha = condition.alpha();
                const double beta  = condition.beta();

#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
                //u05[0][i] = (3.5*u05[1][i] - 2.0*u05[2][i] + 0.5*u05[3][i] + hy*(1.0/beta)*value)/(2.0 + (alpha/beta)*hy);
                u05[0][i] = (4.0*beta*u05[1][i] - beta*u05[2][i] + 2.0*hy*value)/(2.0*alpha*hy+3.0*beta);
#else
                u05[0][i] = (beta*u05[1][i] + hy*value)/(beta + alpha*hy);
#endif
            }

            sn.j = ymax; sn.y = ymax*hy;
#ifdef PARABOLIC_IBVP_H_D2V1_BR_X
            value = boundary(sn, tn00, condition);
#else
            value = 0.5 * ( boundary(sn, tn00, condition) + boundary(sn, tn10, condition) );
#endif
            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                u05[M][i] = value;
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
                u05[M][i] = (3.5*u05[M-1][i] - 2.0*u05[M-2][i] + 0.5*u05[M-3][i] + hy*value)/(2.0);
#else
                u05[M][i] = u05[M-1][i] + hy*value;
#endif
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                const double alpha = condition.alpha();
                const double beta  = condition.beta();

#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
                //u05[M][i] = (3.5*u05[M-1][i] - 2.0*u05[M-2][i] + 0.5*u05[M-3][i] + hy*(1.0/beta)*value)/(2.0 + (alpha/beta)*hy);
                u05[M][i] = (4.0*beta*u05[M-1][i] - beta*u05[M-2][i] + 2.0*hy*value)/(2.0*alpha*hy+3.0*beta);
#else
                u05[M][i] = (beta*u05[M-1][i] + hy*value)/(beta + alpha*hy);
#endif
            }
        }

        /**************************************************** x direction apprx ***************************************************/

        /**************************************************** y direction apprx ***************************************************/

        for (n=xmin+1, sn.i=n, sn.x=n*hx, i=1; n<=xmax-1; ++n, sn.i=n, sn.x=n*hx, ++i)
        {
            for (m=ymin+1, sn.j=m, sn.y=m*hy, j=1; m<=ymax-1; ++m, sn.j=m, sn.y=m*hy, ++j)
            {
#ifdef PARABOLIC_IBVP_H_D2V1_FX_Y
                double fx = f(sn, tn10);
#else
                double fx = 0.5 * ( f(sn, tn00)+f(sn, tn10) );
#endif
                dy[j] = k24*u05[j][i-1] + k25*u05[j][i] + k26*u05[j][i+1] + ht05*fx;
            }

            sn.j = ymin; sn.y = ymin*hy;
#ifdef PARABOLIC_IBVP_H_D2V1_BR_Y
            value = boundary(sn, tn10, condition);
#else
            value = 0.5 * ( boundary(sn, tn00, condition) + boundary(sn, tn10, condition) );
#endif
            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                s = 1;
                u10[0][i] = value;
                dy[1] -= k21 * u10[0][i];
                ay[1] = ay[0] = by[0] = cy[0] = dy[0] = ry[0] = 0.0;
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
                s = 0;
#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
#ifdef PARABOLIC_IBVP_H_D2V1_FX_Y
                double fx = f(sn, tn10);
#else
                double fx = 0.5 * ( f(sn, tn00)+f(sn, tn10) );
#endif
                ay[0]  = 0.0;
                by[0]  = b211;
                cy[0]  = b213;
                dy[0]  = b214 * u05[0][i];
                dy[0] += b215 * ((u05[0][i+1]-2.0*u05[0][i]+u05[0][i-1])/(hx*hx));
                dy[0] += b216 * ((u05[0][i+1]-u05[0][i-1])/(2.0*hx));
                dy[0] += ht05 * fx;
                dy[0] += b217 * value;
#else
                ay[0]  = +0.0;
                by[0]  = +1.0;
                cy[0]  = -1.0;
                dy[0]  = +hx * value;
#endif
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                const double alpha = condition.alpha();
                const double beta  = condition.beta();

                s = 0;
#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
#ifdef PARABOLIC_IBVP_H_D2V1_BR_Y
                double fx = f(sn, tn10);
#else
                double fx = 0.5 * ( f(sn, tn00)+f(sn, tn10) );
#endif
                ay[s]  = 0.0;
                by[s]  = beta * b211 + alpha * b212;
                cy[s]  = beta * b213;
                dy[s]  = beta * b214 * u05[s][i];
                dy[s] += beta * b215 * ((u05[s][i+1]-2.0*u05[s][i]+u05[s][i-1])/(hx*hx));
                dy[s] += beta * b216 * ((u05[s][i+1]-u05[s][i-1])/(2.0*hx));
                dy[s] += beta * ht05 * fx;
                dy[s] += b217 * value;
#else
                ay[s]  = +0.0;
                by[s]  = +beta + alpha * hx;
                cy[s]  = -beta;
                dy[s]  = +hx * value;
#endif
            }

            sn.j = ymax; sn.y = ymax*hy;
#ifdef PARABOLIC_IBVP_H_D2V1_BR_Y
            value = boundary(sn, tn10, condition);
#else
            value = 0.5 * ( boundary(sn, tn00, condition) + boundary(sn, tn10, condition) );
#endif
            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                e = M-1;
                u10[M][i] = value;
                dy[e] -= k23 * u10[M][i];
                cy[e] = ay[M] = by[M] = cy[M] = dy[M] = ry[M] = 0.0;
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
                e = M;
#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
#ifdef PARABOLIC_IBVP_H_D2V1_FX_Y
                double fx = f(sn, tn10);
#else
                double fx = 0.5 * ( f(sn, tn00)+f(sn, tn10) );
#endif
                ay[e]  = b221;
                by[e]  = b223;
                cy[e]  = 0.0;
                dy[e]  = b224 * u05[e][i];
                dy[e] += b225 * ((u05[e][i+1]-2.0*u05[e][i]+u05[e][i-1])/(hx*hx));
                dy[e] += b226 * ((u05[e][i+1]-u05[e][i-1])/(2.0*hx));
                dy[e] += ht05 * fx;
                dy[e] += b227 * value;
#else
                ay[e]  = -1.0;
                by[e]  = +1.0;
                cy[e]  = +0.0;
                dy[e]  = +hx * value;
#endif
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                const double alpha = condition.alpha();
                const double beta  = condition.beta();

                e = M;
#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
#ifdef PARABOLIC_IBVP_H_D2V1_FX_Y
                double fx = f(sn, tn10);
#else
                double fx = 0.5 * ( f(sn, tn00)+f(sn, tn10) );
#endif
                ay[e]  = beta * b221;
                by[e]  = beta * b223 + alpha * b222;
                cy[e]  = 0.0;
                dy[e]  = beta * b224 * u05[e][i];
                dy[e] += beta * b225 * ((u05[e][i+1]-2.0*u05[e][i]+u05[e][i-1])/(hx*hx));
                dy[e] += beta * b226 * ((u05[e][i+1]-u05[e][i-1])/(2.0*hx));
                dy[e] += beta * ht05 * fx;
                dy[e] += b227 * value;
#else
                ay[e]  = -beta;
                by[e]  = +beta + alpha * hx;
                cy[e]  = +0.0;
                dy[e]  = +hx * value;
#endif
            }

            tomasAlgorithmLeft2Right(ay+s, by+s, cy+s, dy+s, ry+s, static_cast<size_t>(e-s+1));
            for (size_t j=s; j<=e; j++) u10[j][i] = ry[j];
            //memcpy(u10[i]+s, ry+s, sizeof(double*)*(e-s+1));
        }

        for (m=ymin, sn.j=m, sn.y=m*hy, j=0; m<=ymax; ++m, sn.j=m, sn.y=m*hy, ++j)
        {
            sn.i = xmin; sn.x = xmin*hx;
#ifdef PARABOLIC_IBVP_H_D2V1_BR_Y
            value = boundary(sn, tn10, condition);
#else
            value = 0.5 * ( boundary(sn, tn00, condition) + boundary(sn, tn10, condition) );
#endif
            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                u10[j][0] = value;
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
                u10[j][0] = (3.5*u10[j][1] - 2.0*u10[j][2] + 0.5*u10[j][3] + hx*value)/(2.0);
#else
                u10[j][0] = u10[j][1] + hx*value;
#endif
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                const double alpha = condition.alpha();
                const double beta  = condition.beta();

#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
                //u10[j][0] = (3.5*u10[j][1] - 2.0*u10[j][2] + 0.5*u10[j][3] + hx*(1.0/beta)*value)/(2.0 + (alpha/beta)*hx);
                u10[j][0] = (4.0*beta*u10[j][1] - beta*u10[j][2] + 2.0*hx*value)/(2.0*alpha*hx+3.0*beta);
#else
                u10[j][0] = (beta*u10[j][1] + hx*value)/(beta + alpha*hx);
#endif
            }

            sn.i = xmax; sn.x = xmax*hx;
#ifdef PARABOLIC_IBVP_H_D2V1_BR_Y
            value = boundary(sn, tn10, condition);
#else
            value = 0.5 * ( boundary(sn, tn00, condition) + boundary(sn, tn10, condition) );
#endif
            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                u10[j][N] = value;
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
                u10[j][N] = (3.5*u10[j][N-1] - 2.0*u10[j][N-2] + 0.5*u10[j][N-3] + hx*value)/(2.0);
#else
                u10[j][N] = u10[j][N-1] + hx*value;
#endif
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                const double alpha = condition.alpha();
                const double beta  = condition.beta();

#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
                //u10[j][N] = (3.5*u10[j][N-1] - 2.0*u10[j][N-2] + 0.5*u10[j][N-3] + hx*(1.0/beta)*value)/(2.0 + (alpha/beta)*hx);
                u10[j][N] = (4.0*beta*u10[j][N-1] - beta*u10[j][N-2] + 2.0*hx*value)/(2.0*alpha*hx+3.0*beta);
#else
                u10[j][N] = (beta*u10[j][N-1] + hx*value)/(beta + alpha*hx);
#endif
            }
        }

        layerInfo(DoubleMatrix(u10, M+1, N+1), tn10);

        /**************************************************** y direction apprx ***************************************************/
        double **_tmp = u00; u00 = u10; u10 = _tmp;
    }

    for (size_t i=0; i<=M; i++) { free(u00[i]); free(u05[i]); free(u10[i]); }
    free(u00); free(u05); free(u10);

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

void IHeatEquationFBVP::explicit_calculate_D1V1() const
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

    if (ht >= 0.5*(hx*hx)) throw std::runtime_error("Differential scheme not steady!");

    // common parameters
    const double k1 = -(a*ht)/(hx*hx);
    const double k2 = -(b*ht)/(2.0*hx);
    const double k3 = 1.0 - c*ht;

    DoubleVector u0(N+1, 0.0);
    DoubleVector u1(N+1, 0.0);

    /***********************************************************************************************/
    /***********************************************************************************************/

    TimeNodePDE tn;
    SpaceNodePDE sn;
    unsigned int i=0;
    int n = 0;
    BoundaryConditionPDE condition; double value, alpha, beta, gamma;

    for (n=xmin, sn.i=n, sn.x=n*hx, i=0; n<=xmax; ++n, sn.i=n, sn.x=n*hx, ++i)
    {
        u0[i] = final(sn, FinalCondition::FinalValue);
    }

    tn.i = 0; tn.t = tn.i*ht;
    layerInfo(u0, tn);

    /***********************************************************************************************/

    for (unsigned int ln=L-1, size_ln = static_cast<unsigned int>(0)-1; ln != size_ln; ln--)
    {
        tn.i = ln; tn.t = tn.i*ht;

        for (n=xmin+1, sn.i=n, sn.x=n*hx, i=1; n<=xmax-1; ++n, sn.i=n, sn.x=i*hx, ++i)
        {
            u1[i] = k3*u0[i] + k1*(u0[i+1] - 2.0*u0[i] + u0[i-1]) + k2*(u0[i+1] - u0[i-1]) - ht*f(sn, tn);
        }

        sn.i = xmin; sn.x = xmin*hx;
        value = boundary(sn, tn, condition);

        if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
        {
            u1[0] = /*(gamma/alpha)**/value;
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
        {
            u1[0] = u1[1] - hx*/*(gamma/beta)**/value;
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Robin)
        {
            const double alpha = condition.alpha();
            const double beta  = condition.beta();
            //const double gamma = condition.gamma();
            const double gamma = 1.0;

            u1[0] = (beta*u1[1] - hx*gamma*value)/(beta-hx*alpha);
        }

        sn.i = xmax; sn.x = xmax*hx;
        value = boundary(sn, tn, condition);

        if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
        {
            u1[N] = /*(gamma/alpha)**/value;
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
        {
            u1[N] = u1[N-1] + hx*/*(gamma/beta)**/value;
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Robin)
        {
            const double alpha = condition.alpha();
            const double beta  = condition.beta();
            //const double gamma = condition.gamma();
            const double gamma = 1.0;

            u1[N] = (beta*u1[N-1] + hx*gamma*value)/(beta+hx*alpha);
        }

        layerInfo(u1, tn);

        u0 = u1;
    }

    u0.clear();
    u1.clear();
}

void IHeatEquationFBVP::implicit_calculate_D1V1() const
{
    const size_t N = static_cast<size_t>(spaceDimensionX().size()) - 1;
    const size_t M = static_cast<size_t>(timeDimension().size()) - 1;

    const int xmin = spaceDimensionX().min();
    const int xmax = spaceDimensionX().max();
    const int tmin = timeDimension().min();
    //const int tmax = timeDimension().max();

    const double hx = spaceDimensionX().step();
    const double ht = timeDimension().step();

    const double a = thermalDiffusivity();
    const double b = thermalConductivity();
    const double c = thermalConvection();
    const double w1 = weight();
    const double w2 = 1.0 - w1;

    // equation parameters
    const double k11 = -((a*ht)/(hx*hx))*w1 + ((b*ht)/(2.0*hx))*w1; // i-1
    const double k12 = -1.0 + ((2.0*a*ht)/(hx*hx))*w1 - c*ht*w1;    // i
    const double k13 = -((a*ht)/(hx*hx))*w1 - ((b*ht)/(2.0*hx))*w1; // i+1
    const double k21 = +((a*ht)/(hx*hx))*w2 - ((b*ht)/(2.0*hx))*w2; // i-1
    const double k22 = -1.0 - ((2.0*a*ht)/(hx*hx))*w2 + c*ht*w2;    // i
    const double k23 = +((a*ht)/(hx*hx))*w2 + ((b*ht)/(2.0*hx))*w2; // i+1

    // left border condition parameters
    const double b11 = -1.0 + ((2.0*a*ht)/(hx*hx))*w1 - c*ht*w1;
    const double b12 = -((2.0*a*ht)/hx)*w1 + b*ht*w1;
    const double b13 = -((2.0*a*ht)/(hx*hx))*w1;
    const double b14 = -1.0 - ((2.0*a*ht)/(hx*hx))*w2 + c*ht*w2;
    const double b15 = +((2.0*a*ht)/hx)*w2 - b*ht*w2;
    const double b16 = +((2.0*a*ht)/(hx*hx))*w2;
    const double b19 = -((2.0*a*ht)/hx) + b*ht;
    //const double b17 = -((2.0*a*ht)/hx)*w1 + b*ht*w1;
    //const double b18 = -((2.0*a*ht)/hx)*w2 + b*ht*w2;

    // right border condition parameters
    const double b21 = -((2.0*a*ht)/(hx*hx))*w1;
    const double b22 = +((2.0*a*ht)/hx)*w1 + b*ht*w1;
    const double b23 = -1.0 + ((2.0*a*ht)/(hx*hx))*w1 - c*ht*w1;
    const double b24 = +((2.0*a*ht)/(hx*hx))*w2;
    const double b25 = -((2.0*a*ht)/hx)*w2 - b*ht*w2;
    const double b26 = -1.0 - ((2.0*a*ht)/(hx*hx))*w2 + c*ht*w2;
    const double b29 = +((2.0*a*ht)/hx) + b*ht;
    const double b27 = +((2.0*a*ht)/hx)*w1 + b*ht*w1;
    const double b28 = +((2.0*a*ht)/hx)*w2 + b*ht*w2;

    // common parameters
    //const double ht_w1 = +ht*w1;
    //const double ht_w2 = +ht*w2;

    double *ax = static_cast<double*>(malloc(sizeof(double)*(N+1)));
    double *bx = static_cast<double*>(malloc(sizeof(double)*(N+1)));
    double *cx = static_cast<double*>(malloc(sizeof(double)*(N+1)));
    double *dx = static_cast<double*>(malloc(sizeof(double)*(N+1)));
    double *rx = static_cast<double*>(malloc(sizeof(double)*(N+1)));

    for (size_t n=0; n<=N; n++)
    {
        ax[n] = k11;
        bx[n] = k12;
        cx[n] = k13;
    }
    ax[0] = 0.0; cx[N] = 0.0;

    double *p0 = static_cast<double*>(malloc(sizeof(double)*(N+1)));
    double *p1 = static_cast<double*>(malloc(sizeof(double)*(N+1)));

    /***********************************************************************************************/
    /***********************************************************************************************/

    SpaceNodePDE sn;
    unsigned int i = 0;
    unsigned int k = static_cast<unsigned int>(timeDimension().max());
    for (int n=xmin; n<=xmax; n++, i++)
    {
        sn.i = n; sn.x = sn.i*hx;
        p0[i] = final(sn, FinalCondition::FinalValue);
    }
    TimeNodePDE tn; tn.i = k; tn.t = tn.i*ht; k--;
    layerInfo(DoubleVector(p0, N+1), tn);

    /***********************************************************************************************/
    /***********************************************************************************************/

    BoundaryConditionPDE condition; double value;

    for (size_t ln=M-1, size_ln = static_cast<size_t>(0)-1; ln != size_ln; ln--)
    {
        TimeNodePDE tn0; tn0.i = static_cast<unsigned int>(k+1); tn0.t = tn0.i*ht;
        TimeNodePDE tn1; tn1.i = static_cast<unsigned int>(k  ); tn1.t = tn1.i*ht;

        size_t i = 1;
        for (int n=xmin+1; n<=xmax-1; n++, i++)
        {
            sn.i = n; sn.x = n*hx;

#ifdef PARABOLIC_IBVP_H_D1V1_FX
            double fx = f(sn, tn1);
#else
            double fx = w1*f(sn, tn1) + w2*f(sn, tn0);
#endif
            dx[i] = k21*p0[i-1] + k22*p0[i] + k23*p0[i+1] + ht*fx;
        }

        size_t s=0, e=N;

        sn.i = xmin; sn.x = xmin*hx;
#ifdef PARABOLIC_IBVP_H_D1V1_BR
        value = boundary(sn, tn1, condition);
#else
        value = w1*boundary(sn, tn1, condition) + w2*boundary(sn, tn0, condition);
#endif

        if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
        {
            s = 1;
            p1[0] = value;
            dx[1] -= k11 * p1[0];
            ax[1] = ax[0] = bx[0] = cx[0] = dx[0] = rx[0] = 0.0;
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
        {
            s = 0;
#ifdef PARABOLIC_IBVP_H_D1V1_BORDER_O2
#ifdef PARABOLIC_IBVP_H_D1V1_FX
            double fx = f(sn, tn1);
#else
            double fx = w1*f(sn, tn1) + w2*f(sn, tn0);
#endif
            ax[s]  = 0.0;
            bx[s]  = b11;
            cx[s]  = b13;
            dx[s]  = b14*p0[s] + b16*p0[s+1];
            dx[s] += b19*value;
            dx[s] += ht*fx;
#endif
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Robin)
        {
            const double alpha = condition.alpha();
            const double beta  = condition.beta();

            s = 0;
#ifdef PARABOLIC_IBVP_H_D1V1_BORDER_O2
#ifdef PARABOLIC_IBVP_H_D1V1_FX
            double fx = f(sn, tn1);
#else
            double fx = w1*f(sn, tn1) + w2*f(sn, tn0);
#endif
            ax[s]  = 0.0;
            bx[s]  = beta * b11 + alpha * b12;
            cx[s]  = beta * b13;
            dx[s]  = beta * b14*p0[s] + alpha * b15*p0[s] + beta * b16*p0[s+1];
            dx[s] += b19*value;
            dx[s] += beta*ht*fx;
#endif
        }

        sn.i = xmax; sn.x = xmax*hx;
#ifdef PARABOLIC_IBVP_H_D1V1_BR
        value = boundary(sn, tn1, condition);
#else
        value = w1*boundary(sn, tn1, condition) + w2*boundary(sn, tn0, condition);
#endif

        if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
        {
            e = N-1;
            p1[N] = value;
            dx[N-1] -= k13 * p1[N];
            cx[N-1] = ax[N] = bx[N] = cx[N] = dx[N] = rx[N] = 0.0;
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
        {
            e = N;
#ifdef PARABOLIC_IBVP_H_D1V1_BORDER_O2
#ifdef PARABOLIC_IBVP_H_D1V1_FX
            double fx = f(sn, tn1);
#else
            double fx = w1*f(sn, tn1) + w2*f(sn, tn0);
#endif
            ax[e]  = b21;
            bx[e]  = b23;
            cx[e]  = 0.0;
            dx[e]  = b24*p0[e-1] + b26*p0[e];
            dx[e] += b29*value;
            dx[e] += ht*fx;
#endif
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Robin)
        {
            const double alpha = condition.alpha();
            const double beta  = condition.beta();

            e = N;
#ifdef PARABOLIC_IBVP_H_D1V1_BORDER_O2
#ifdef PARABOLIC_IBVP_H_D1V1_FX
            double fx = f(sn, tn1);
#else
            double fx = w1*f(sn, tn1) + w2*f(sn, tn0);
#endif
            ax[e]  = beta*b21;
            bx[e]  = alpha*b22 + beta*b23;
            cx[e]  = 0.0;
            dx[e]  = beta*b24*p0[e-1] + alpha*b25*p0[e] + beta*b26*p0[e];
            dx[e] += b29*value;
            dx[e] += beta*ht*fx;
#endif
        }

        tomasAlgorithmLeft2Right(ax+s, bx+s, cx+s, dx+s, rx+s, static_cast<unsigned int>(e-s+1));
        memcpy(p1+s, rx+s, sizeof(double)*(e-s+1));
        //for (unsigned int n=s; n<=e; n++) u1[n] = rx[n];

        layerInfo(DoubleVector(p1, N+1), tn1);

        double *tmp_p = p0; p0 = p1; p1 = tmp_p;
    }

    free(p0);
    free(p1);

    free(rx);
    free(dx);
    free(cx);
    free(bx);
    free(ax);
}

void IHeatEquationFBVP::implicit_calculate_D2V1() const
{
    const size_t N = static_cast<size_t>(spaceDimensionX().size()) - 1;
    const size_t M = static_cast<size_t>(spaceDimensionY().size()) - 1;
    const size_t L = static_cast<size_t>(timeDimension().size()) - 1;

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

    //    double l1 = hx*fabs(xmax-xmin);
    //    double l2 = hx*fabs(xmax-xmin);
    //    double ht_max = (1.0/M_PI) * sqrt((a1*a1)/(hx*hx)+(a2*a2)/(hy*hy)) / sqrt((1.0/(l1*l1))+(1.0/(l2*l2)));
    //    if (ht > ht_max) { throw std::runtime_error("Differential scheme is conditionally steady."); }

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

    for (size_t n=0; n<=N; n++)
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

    for (size_t m=0; m<=M; m++)
    {
        ay[m] = k21;
        by[m] = k22;
        cy[m] = k23;
    }
    ay[0] = 0.0; cy[M] = 0.0;

    double **p00 = static_cast<double**>(malloc(sizeof(double*)*(M+1)));
    double **p05 = static_cast<double**>(malloc(sizeof(double*)*(M+1)));
    double **p10 = static_cast<double**>(malloc(sizeof(double*)*(M+1)));
    for (size_t i=0; i<=M; i++)
    {
        p00[i] = static_cast<double*>(malloc(sizeof(double)*(N+1)));
        p05[i] = static_cast<double*>(malloc(sizeof(double)*(N+1)));
        p10[i] = static_cast<double*>(malloc(sizeof(double)*(N+1)));
    }

    /***********************************************************************************************/
    /***********************************************************************************************/

    TimeNodePDE tn00, tn10;
    SpaceNodePDE sn;
    size_t i=0, j=0, s=0, e=N;
    int n = 0, m = 0;

    for (m=ymin, sn.j=m, sn.y=m*hy, j=0; m<=ymax; ++m, sn.j=m, sn.y=m*hy, ++j)
    {
        for (n=xmin, sn.i=n, sn.x=n*hx, i=0; n<=xmax; ++n, sn.i=n, sn.x=n*hx, ++i)
        {
            p00[j][i] = final(sn, FinalCondition::FinalValue);
        }
    }

    tn00.i = static_cast<unsigned int>(L); tn00.t = tn00.i*ht;
    layerInfo(DoubleMatrix(p00, M+1, N+1), tn00);

    /***********************************************************************************************/
    BoundaryConditionPDE condition; double value;

    for (size_t ln=L-1, size_ln = static_cast<size_t>(0)-1; ln != size_ln; ln--)
    {
        tn00.i = static_cast<unsigned int>(ln+1); tn00.t = tn00.i*ht;
        tn10.i = static_cast<unsigned int>(ln-0); tn10.t = tn10.i*ht;

        /**************************************************** x direction apprx ***************************************************/

        for (m=ymin+1, sn.j=m, sn.y=m*hy, j=1; m<=ymax-1; ++m, sn.j=m, sn.y=m*hy, ++j)
        {
            for (n=xmin+1, sn.i=n, sn.x=n*hx, i=1; n<=xmax-1; ++n, sn.i=n, sn.x=n*hx, ++i)
            {
#ifdef PARABOLIC_IBVP_H_D2V1_FX_X
                double fx = f(sn, tn00);
#else
                double fx = 0.5 * ( f(sn, tn00) + f(sn, tn10) );
#endif
                dx[i] = k14*p00[j-1][i] + k15*p00[j][i] + k16*p00[j+1][i] + ht05*fx;
            }

            sn.i = xmin; sn.x = xmin*hx;
#ifdef PARABOLIC_IBVP_H_D2V1_BR_X
            value = boundary(sn, tn00, condition);
#else
            value = 0.5*( boundary(sn, tn00, condition) + boundary(sn, tn10, condition) );
#endif
            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                s = 1;
                p05[j][0] = value;
                dx[1] -= k11 * p05[j][0];
                ax[1] = ax[0] = bx[0] = cx[0] = dx[0] = rx[0] = 0.0;
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
                s = 0;
#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
#ifdef PARABOLIC_IBVP_H_D2V1_FX_X
                double fx = f(sn, tn00);
#else
                double fx = 0.5 * ( f(sn, tn00) + f(sn, tn10) );
#endif
                ax[0]  = 0.0;
                bx[0]  = b111;
                cx[0]  = b113;
                dx[0]  = b114 * p00[j][0];
                dx[0] += b115 * ((p00[j+1][0]-2.0*p00[j][0]+p00[j-1][0])/(hy*hy));
                dx[0] += b116 * ((p00[j+1][0]-p00[j-1][0])/(2.0*hy));
                dx[0] += ht05 * fx;
                dx[0] += b117 * value;
#else
                ax[0]  = +0.0;
                bx[0]  = +1.0;
                cx[0]  = -1.0;
                dx[0]  = +hx * value;
#endif
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                const double alpha = condition.alpha();
                const double beta  = condition.beta();

                s = 0;
#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
#ifdef PARABOLIC_IBVP_H_D2V1_FX_X
                double fx = f(sn, tn00);
#else
                double fx = 0.5 * ( f(sn, tn00) + f(sn, tn10) );
#endif
                ax[0]  = 0.0;
                bx[0]  = beta * b111 + alpha * b112;
                cx[0]  = beta * b113;
                dx[0]  = beta * b114 * p00[j][0];
                dx[0] += beta * b115 * ((p00[j+1][0]-2.0*p00[j][0]+p00[j-1][0])/(hy*hy));
                dx[0] += beta * b116 * ((p00[j+1][0]-p00[j-1][0])/(2.0*hy));
                dx[0] += beta * ht05 * fx;
                dx[0] += b117 * value;
#else
                ax[0]  = +0.0;
                bx[0]  = +beta + alpha * hx;
                cx[0]  = -beta;
                dx[0]  = +hx * value;
#endif
            }

            sn.i = xmax; sn.x = xmax*hx;
#ifdef PARABOLIC_IBVP_H_D2V1_BR_X
            value = boundary(sn, tn00, condition);
#else
            value = 0.5 * ( boundary(sn, tn00, condition) + boundary(sn, tn10, condition) );
#endif
            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                e = N-1;
                p05[j][N] = value;
                dx[e] -= k13 * p05[j][N];
                cx[e] = ax[N] = bx[N] = cx[N] = dx[N] = rx[N] = 0.0;
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
                e = N;
#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
#ifdef PARABOLIC_IBVP_H_D2V1_FX_X
                double fx = f(sn, tn00);
#else
                double fx = 0.5 * ( f(sn, tn00) + f(sn, tn10) );
#endif
                ax[e]  = b121;
                bx[e]  = b123;
                cx[e]  = 0.0;
                dx[e]  = b124 * p00[j][e];
                dx[e] += b125 * ((p00[j+1][e]-2.0*p00[j][e]+p00[j-1][e])/(hy*hy));
                dx[e] += b126 * ((p00[j+1][e]-p00[j-1][e])/(2.0*hy));
                dx[e] += ht05 * fx;
                dx[e] += b127 * value;
#else
                ax[e]  = -1.0;
                bx[e]  = +1.0;
                cx[e]  = +0.0;
                dx[e]  = hx * value;
#endif
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                const double alpha = condition.alpha();
                const double beta  = condition.beta();

                e = N;
#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
#ifdef PARABOLIC_IBVP_H_D2V1_FX_X
                double fx = f(sn, tn00);
#else
                double fx = 0.5 * ( f(sn, tn00) + f(sn, tn10) );
#endif
                ax[e]  = beta * b121;
                bx[e]  = beta * b123 + alpha * b122;
                cx[e]  = 0.0;
                dx[e]  = beta * b124 * p00[j][e];
                dx[e] += beta * b125 * ((p00[j+1][e]-2.0*p00[j][e]+p00[j-1][e])/(hy*hy));
                dx[e] += beta * b126 * ((p00[j+1][e]-p00[j-1][e])/(2.0*hy));
                dx[e] += beta * ht05 * fx;
                dx[e] += b127 * value;
#else
                ax[e]  = -beta;
                bx[e]  = +beta + alpha * hx;
                cx[e]  = +0.0;
                dx[e]  = +hx * value;
#endif
            }

            tomasAlgorithmLeft2Right(ax+s, bx+s, cx+s, dx+s, rx+s, static_cast<unsigned int>(e-s+1));
            //for (unsigned int i=s; i<=e; i++) p05[j][i] = rx[i];
            memcpy(p05[j]+s, rx+s, sizeof (double*)*(e-s+1));
        }

        for (n=xmin, sn.i=n, sn.x=n*hx, i=0; n<=xmax; ++n, sn.i=n, sn.x=n*hx, ++i)
        {
            sn.j = ymin; sn.y = ymin*hy;
#ifdef PARABOLIC_IBVP_H_D2V1_BR_X
            value = boundary(sn, tn00, condition);
#else
            value = 0.5 * ( boundary(sn, tn00, condition) + boundary(sn, tn10, condition) );
#endif
            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                p05[0][i] = value;
            }
            if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
                p05[0][i] = (3.5*p05[1][i] - 2.0*p05[2][i] + 0.5*p05[3][i] + hy*value)/(2.0);
#else
                p05[0][i] = p05[1][i] + hy*value;
#endif
            }
            if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                const double alpha = condition.alpha();
                const double beta  = condition.beta();

#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
                p05[0][i] = (3.5*p05[1][i] - 2.0*p05[2][i] + 0.5*p05[3][i] + hy*(1.0/beta)*value)/(2.0 + (alpha/beta)*hy);
#else
                p05[0][i] = (beta*p05[1][i] + hy*value)/(beta + alpha*hy);
#endif
            }

            sn.j = ymax; sn.y = ymax*hy;
#ifdef PARABOLIC_IBVP_H_D2V1_BR_X
            value = boundary(sn, tn00, condition);
#else
            value = 0.5 * ( boundary(sn, tn00, condition) + boundary(sn, tn10, condition) );
#endif
            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                p05[M][i] = value;
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
                p05[M][i] = (3.5*p05[M-1][i] - 2.0*p05[M-2][i] + 0.5*p05[M-3][i] + hy*value)/(2.0);
#else
                p05[M][i] = p05[M-1][i] + hy*value;
#endif
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                const double alpha = condition.alpha();
                const double beta  = condition.beta();

#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
                p05[M][i] = (3.5*p05[M-1][i] - 2.0*p05[M-2][i] + 0.5*p05[M-3][i] + hy*(1.0/beta)*value)/(2.0 + (alpha/beta)*hy);
#else
                p05[M][i] = (beta*p05[M-1][i] + hy*value)/(beta + alpha*hy);
#endif
            }
        }

        /**************************************************** x direction apprx ***************************************************/

        /**************************************************** y direction apprx ***************************************************/

        for (n=xmin+1, sn.i=n, sn.x=n*hx, i=1; n<=xmax-1; ++n, sn.i=n, sn.x=n*hx, ++i)
        {
            for (m=ymin+1, sn.j=m, sn.y=m*hy, j=1; m<=ymax-1; ++m, sn.j=m, sn.y=m*hy, ++j)
            {
#ifdef PARABOLIC_IBVP_H_D2V1_FX_Y
                double fx = f(sn, tn10);
#else
                double fx = 0.5 * ( f(sn, tn00)+f(sn, tn10) );
#endif
                dy[j] = k24*p05[j][i-1] + k25*p05[j][i] + k26*p05[j][i+1] + ht05*fx;
            }

            sn.j = ymin; sn.y = ymin*hy;
#ifdef PARABOLIC_IBVP_H_D2V1_BR_Y
            value = boundary(sn, tn10, condition);
#else
            value = 0.5 * ( boundary(sn, tn00, condition) + boundary(sn, tn10, condition) );
#endif
            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                s = 1;
                p10[0][i] = value;
                dy[1] -= k21 * p10[0][i];
                ay[1] = ay[0] = by[0] = cy[0] = dy[0] = ry[0] = 0.0;
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
                s = 0;
#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
#ifdef PARABOLIC_IBVP_H_D2V1_FX_Y
                double fx = f(sn, tn10);
#else
                double fx = 0.5 * ( f(sn, tn00)+f(sn, tn10) );
#endif
                ay[0]  = 0.0;
                by[0]  = b211;
                cy[0]  = b213;
                dy[0]  = b214 * p05[0][i];
                dy[0] += b215 * ((p05[0][i+1]-2.0*p05[0][i]+p05[0][i-1])/(hx*hx));
                dy[0] += b216 * ((p05[0][i+1]-p05[0][i-1])/(2.0*hx));
                dy[0] += ht05 * fx;
                dy[0] += b217 * value;
#else
                ay[0]  = +0.0;
                by[0]  = +1.0;
                cy[0]  = -1.0;
                dy[0]  = +hx * value;
#endif
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                const double alpha = condition.alpha();
                const double beta  = condition.beta();

                s = 0;
#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
#ifdef PARABOLIC_IBVP_H_D2V1_BR_Y
                double fx = f(sn, tn10);
#else
                double fx = 0.5 * ( f(sn, tn00)+f(sn, tn10) );
#endif
                ay[s]  = 0.0;
                by[s]  = beta * b211 + alpha * b212;
                cy[s]  = beta * b213;
                dy[s]  = beta * b214 * p05[s][i];
                dy[s] += beta * b215 * ((p05[s][i+1]-2.0*p05[s][i]+p05[s][i-1])/(hx*hx));
                dy[s] += beta * b216 * ((p05[s][i+1]-p05[s][i-1])/(2.0*hx));
                dy[s] += beta * ht05 * fx;
                dy[s] += b217 * value;
#else
                ay[s]  = +0.0;
                by[s]  = +beta + alpha * hx;
                cy[s]  = -beta;
                dy[s]  = +hx * value;
#endif
            }

            sn.j = ymax; sn.y = ymax*hy;
#ifdef PARABOLIC_IBVP_H_D2V1_BR_Y
            value = boundary(sn, tn10, condition);
#else
            value = 0.5 * ( boundary(sn, tn00, condition) + boundary(sn, tn10, condition) );
#endif
            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                e = M-1;
                p10[M][i] = value;
                dy[e] -= k23 * p10[M][i];
                cy[e] = ay[M] = by[M] = cy[M] = dy[M] = ry[M] = 0.0;
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
                e = M;
#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
#ifdef PARABOLIC_IBVP_H_D2V1_FX_Y
                double fx = f(sn, tn10);
#else
                double fx = 0.5 * ( f(sn, tn00)+f(sn, tn10) );
#endif
                ay[e]  = b221;
                by[e]  = b223;
                cy[e]  = 0.0;
                dy[e]  = b224 * p05[e][i];
                dy[e] += b225 * ((p05[e][i+1]-2.0*p05[e][i]+p05[e][i-1])/(hx*hx));
                dy[e] += b226 * ((p05[e][i+1]-p05[e][i-1])/(2.0*hx));
                dy[e] += ht05 * fx;
                dy[e] += b227 * value;
#else
                ay[e]  = -1.0;
                by[e]  = +1.0;
                cy[e]  = +0.0;
                dy[e]  = +hx * value;
#endif
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                const double alpha = condition.alpha();
                const double beta  = condition.beta();

                e = M;
#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
#ifdef PARABOLIC_IBVP_H_D2V1_FX_Y
                double fx = f(sn, tn10);
#else
                double fx = 0.5 * ( f(sn, tn00)+f(sn, tn10) );
#endif
                ay[e]  = beta * b221;
                by[e]  = beta * b223 + alpha * b222;
                cy[e]  = 0.0;
                dy[e]  = beta * b224 * p05[e][i];
                dy[e] += beta * b225 * ((p05[e][i+1]-2.0*p05[e][i]+p05[e][i-1])/(hx*hx));
                dy[e] += beta * b226 * ((p05[e][i+1]-p05[e][i-1])/(2.0*hx));
                dy[e] += beta * ht05 * fx;
                dy[e] += b227 * value;
#else
                ay[e]  = -beta;
                by[e]  = +beta + alpha * hx;
                cy[e]  = 0.0;
                dy[e]  = +hx * value;
#endif
            }

            tomasAlgorithmLeft2Right(ay+s, by+s, cy+s, dy+s, ry+s, static_cast<unsigned int>(e-s+1));
            for (size_t j=s; j<=e; j++) p10[j][i] = ry[j];
            //memcpy(p10[i]+s, ry+s, sizeof(double*)*(e-s+1));
        }

        for (m=ymin, sn.j=m, sn.y=m*hy, j=0; m<=ymax; ++m, sn.j=m, sn.y=m*hy, ++j)
        {
            sn.i = xmin; sn.x = xmin*hx;
#ifdef PARABOLIC_IBVP_H_D2V1_BR_Y
            value = boundary(sn, tn10, condition);
#else
            value = 0.5 * ( boundary(sn, tn00, condition) + boundary(sn, tn10, condition) );
#endif
            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                p10[j][0] = value;
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
                p10[j][0] = (3.5*p10[j][1] - 2.0*p10[j][2] + 0.5*p10[j][3] + hx*value)/(2.0);
#else
                p10[j][0] = p10[j][1] + hx*value;
#endif
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                const double alpha = condition.alpha();
                const double beta  = condition.beta();

#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
                p10[j][0] = (3.5*p10[j][1] - 2.0*p10[j][2] + 0.5*p10[j][3] + hx*(1.0/beta)*value)/(2.0 + (alpha/beta)*hx);
#else
                p10[j][0] = (beta*p10[j][1] + hx*value)/(beta + alpha*hx);
#endif
            }

            sn.i = xmax; sn.x = xmax*hx;
#ifdef PARABOLIC_IBVP_H_D2V1_BR_Y
            value = boundary(sn, tn10, condition);
#else
            value = 0.5 * ( boundary(sn, tn00, condition) + boundary(sn, tn10, condition) );
#endif
            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                p10[j][N] = value;
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
                p10[j][N] = (3.5*p10[j][N-1] - 2.0*p10[j][N-2] + 0.5*p10[j][N-3] + hx*value)/(2.0);
#else
                p10[j][N] = p10[j][N-1] + hx*value;
#endif
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                const double alpha = condition.alpha();
                const double beta  = condition.beta();

#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
                p10[j][N] = (3.5*p10[j][N-1] - 2.0*p10[j][N-2] + 0.5*p10[j][N-3] + hx*(1.0/beta)*value)/(2.0 + (alpha/beta)*hx);
#else
                p10[j][N] = (beta*p10[j][N-1] + hx*value)/(beta + alpha*hx);
#endif
            }
        }

        layerInfo(DoubleMatrix(p10, M+1, N+1), tn10);

        /**************************************************** y direction apprx ***************************************************/
        double **_tmp = p00; p00 = p10; p10 = _tmp;
    }

    for (unsigned int i=0; i<=M; i++) { free(p00[i]); free(p05[i]); free(p10[i]); }
    free(p00); free(p05); free(p10);

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

//------------------------------------------------------------------------------------------------------------------------------------------------------//

ILoadedHeatEquationIBVP::ILoadedHeatEquationIBVP(double thermalDiffusivity, double thermalConductivity, double thermalConvection) :
    IHeatEquationIBVP(thermalDiffusivity, thermalConductivity, thermalConvection) {}

ILoadedHeatEquationIBVP::~ILoadedHeatEquationIBVP() {}

ILoadedHeatEquationIBVP::ILoadedHeatEquationIBVP(const ILoadedHeatEquationIBVP &) {}

ILoadedHeatEquationIBVP& ILoadedHeatEquationIBVP::operator=(const ILoadedHeatEquationIBVP &)
{
    return *this;
}

void ILoadedHeatEquationIBVP::explicit_calculate_D1V1() const {}

void ILoadedHeatEquationIBVP::implicit_calculate_D1V1() const
{
    const size_t N = static_cast<size_t>(spaceDimensionX().size()) - 1;
    const size_t L = static_cast<size_t>(timeDimension().size()) - 1;

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
    const double k11 = -((a*ht)/(hx*hx))*w1 + ((b*ht)/(2.0*hx))*w1; // i-1
    const double k12 = +1.0 + ((2.0*a*ht)/(hx*hx))*w1 - c*ht*w1;    // i
    const double k13 = -((a*ht)/(hx*hx))*w1 - ((b*ht)/(2.0*hx))*w1; // i+1
    const double k21 = +((a*ht)/(hx*hx))*w2 - ((b*ht)/(2.0*hx))*w2; // i-1
    const double k22 = +1.0 - ((2.0*a*ht)/(hx*hx))*w2 + c*ht*w2;    // i
    const double k23 = +((a*ht)/(hx*hx))*w2 + ((b*ht)/(2.0*hx))*w2; // i+1

    // left border condition parameters
    const double b11 = +1.0 + ((2.0*a*ht)/(hx*hx))*w1 - c*ht*w1;
    const double b12 = -((2.0*a*ht)/hx)*w1 + b*ht*w1;
    const double b13 = -((2.0*a*ht)/(hx*hx))*w1;
    const double b14 = +1.0 - ((2.0*a*ht)/(hx*hx))*w2 + c*ht*w2;
    const double b15 = +((2.0*a*ht)/hx)*w2 - b*ht*w2;
    const double b16 = +((2.0*a*ht)/(hx*hx))*w2;
    const double b19 = -((2.0*a*ht)/hx) + b*ht;
    //const double b17 = -((2.0*a*ht)/hx)*w1 + b*ht*w1;
    //const double b18 = -((2.0*a*ht)/hx)*w2 + b*ht*w2;

    // right border condition parameters
    const double b21 = -((2.0*a*ht)/(hx*hx))*w1;
    const double b22 = +((2.0*a*ht)/hx)*w1 + b*ht*w1;
    const double b23 = +1.0 + ((2.0*a*ht)/(hx*hx))*w1 - c*ht*w1;
    const double b24 = +((2.0*a*ht)/(hx*hx))*w2;
    const double b25 = -((2.0*a*ht)/hx)*w2 - b*ht*w2;
    const double b26 = +1.0 - ((2.0*a*ht)/(hx*hx))*w2 + c*ht*w2;
    const double b29 = +((2.0*a*ht)/hx) + b*ht;
    //const double b27 = +((2.0*a*ht)/hx)*w1 + b*ht*w1;
    //const double b28 = +((2.0*a*ht)/hx)*w2 + b*ht*w2;

    // common parameters
    //const double ht_w1 = +ht*w1;
    //const double ht_w2 = +ht*w2;

    double *ax = static_cast<double*>(malloc(sizeof(double)*(N+1)));
    double *bx = static_cast<double*>(malloc(sizeof(double)*(N+1)));
    double *cx = static_cast<double*>(malloc(sizeof(double)*(N+1)));
    double *dx = static_cast<double*>(malloc(sizeof(double)*(N+1)));
    double *rx = static_cast<double*>(malloc(sizeof(double)*(N+1)));

    /////////////// TO-DO

    const size_t lps = loadedPoints().size();
    size_t* minX = new size_t[lps];
    size_t* maxX = new size_t[lps];

    std::vector<size_t> loaded_indecies;
    const size_t k1 = 0;
    for (size_t i=0; i<lps; i++)
    {
        const SpacePoint &p = _loadedPoints[i];
        const size_t rx = static_cast<size_t>(round(p.x/hx));
        minX[i] = rx - k1;
        maxX[i] = rx + k1;
        for (size_t n=minX[i]; n<=maxX[i]; n++)
        {
            if ( std::find(loaded_indecies.begin(), loaded_indecies.end(), n) == loaded_indecies.end() )
            {
                loaded_indecies.push_back(n);
            }
        }
    }
    sort(loaded_indecies.begin(), loaded_indecies.end());

    double **wx = static_cast<double**>(malloc(sizeof(double*)*(N+1)));
    for (size_t i=0; i<=N; i++)
    {
        wx[i] = static_cast<double*>(malloc(sizeof(double)*(N+1)));
        for (size_t j=0; j<=N; j++) wx[i][j] = 0.0;
    }

    //size_t row_size = loaded_indecies.size();

    const double sigma = hx;
    double loadedPart = 0.0;
    for (size_t i=0; i<lps; i++)
    {
        const LoadedSpacePoint &lp = _loadedPoints[i];

        SpacePoint p;
        double uv = 0.0;
        for (size_t n=minX[i]; n<=maxX[i]; n++)
        {
            p.x = n*hx;
            double w = DeltaFunction::gaussian(p, lp, sigma);
            uv += w;
        }
        uv *= hx;
        loadedPart += lp.d * uv;
    }

    /////////////// TO-DO

    for (size_t n=0; n<=N; n++)
    {
        ax[n] = k11;
        bx[n] = k12;
        cx[n] = k13;
    }
    ax[0] = 0.0; cx[N] = 0.0;

    double *u0 = static_cast<double*>(malloc(sizeof(double)*(N+1)));
    double *u1 = static_cast<double*>(malloc(sizeof(double)*(N+1)));

    /***********************************************************************************************/
    /***********************************************************************************************/

    SpaceNodePDE sn;
    unsigned int i = 0;
    unsigned int k = static_cast<unsigned int>(timeDimension().min());
    for (int n=xmin; n<=xmax; n++, i++)
    {
        sn.i = n; sn.x = n*hx;
        u0[i] = initial(sn, InitialCondition::InitialValue);
    }
    TimeNodePDE tn; tn.i = k; tn.t = tn.i*ht; k++;
    layerInfo(DoubleVector(u0, N+1), tn);

    /***********************************************************************************************/
    /***********************************************************************************************/

    BoundaryConditionPDE condition; double value;

    for (size_t ln=1; ln<=L; ln++, k++)
    {
        TimeNodePDE tn0; tn0.i = static_cast<unsigned int>(k-1); tn0.t = tn0.i*ht;
        TimeNodePDE tn1; tn1.i = static_cast<unsigned int>(k  ); tn1.t = tn1.i*ht;

        size_t i = 1;
        for (int n=xmin+1; n<=xmax-1; n++, i++)
        {
            sn.i = n; sn.x = n*hx;

#ifdef PARABOLIC_IBVP_H_D1V1_FX
            double fx = f(sn, tn1);
#else
            double fx = w1*f(sn, tn1) + w2*f(sn, tn0);
#endif
            dx[i] = k21*u0[i-1] + k22*u0[i] + k23*u0[i+1] + ht*fx;
        }

        size_t s=0, e=N;

        sn.i = xmin; sn.x = xmin*hx;
#ifdef PARABOLIC_IBVP_H_D1V1_BR_LEFT
        value = boundary(sn, tn1, condition);
#else
        value = w1*boundary(sn, tn1, condition) + w2*boundary(sn, tn0, condition);
#endif

        if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
        {
            //            s = 1;
            //            u1[0]  = value;
            //            dx[1] -= k11 * u1[0];
            //            ax[1]  = ax[0] = bx[0] = cx[0] = dx[0] = rx[0] = 0.0;

            s = 0;
            ax[0] = 0.0;
            bx[0] = 1.0;
            cx[0] = 0.0;
            dx[0] = value;
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
        {
            s = 0;
#ifdef PARABOLIC_IBVP_H_D1V1_BORDER_O2
#ifdef PARABOLIC_IBVP_H_D1V1_FX
            double fx = f(sn, tn1);
#else
            double fx = w1*f(sn, tn1) + w2*f(sn, tn0);
#endif
            ax[s]  = 0.0;
            bx[s]  = b11;
            cx[s]  = b13;
            dx[s]  = b14*u0[s] + b16*u0[s+1];
            dx[s] += b19*value;
            dx[s] += ht*fx;
#else
            ax[s]  = 0.0;
            bx[s]  = -1.0;
            cx[s]  = +1.0;
            dx[s]  = ht*value;
#endif
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Robin)
        {
            const double alpha = condition.alpha();
            const double beta  = condition.beta();

            s = 0;
#ifdef PARABOLIC_IBVP_H_D1V1_BORDER_O2
#ifdef PARABOLIC_IBVP_H_D1V1_FX
            double fx = f(sn, tn1);
#else
            double fx = w1*f(sn, tn1) + w2*f(sn, tn0);
#endif
            ax[s]  = 0.0;
            bx[s]  = beta*b11 + alpha*b12;
            cx[s]  = beta*b13;
            dx[s]  = beta*b14*u0[s] + alpha*b15*u0[s] + beta*b16*u0[s+1];
            dx[s] += b19*value;
            dx[s] += beta*ht*fx;
#else
            ax[s]  = 0.0;
            bx[s]  = alpha*ht - beta;
            cx[s]  = beta;
            dx[s]  = ht*value;
#endif
        }

        sn.i = xmax; sn.x = xmax*hx;
#ifdef PARABOLIC_IBVP_H_D1V1_BR_RIGHT
        value = boundary(sn, tn1, condition);
#else
        value = w1*boundary(sn, tn1, condition) + w2*boundary(sn, tn0, condition);
#endif

        if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
        {
            //            e = N-1;
            //            u1[N] = value;
            //            dx[N-1] -= k13 * u1[N];
            //            cx[N-1] = ax[N] = bx[N] = cx[N] = dx[N] = rx[N] = 0.0;

            e = N;
            ax[N] = 0.0;
            bx[N] = 1.0;
            cx[N] = 0.0;
            dx[N] = value;
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
        {
            e = N;
#ifdef PARABOLIC_IBVP_H_D1V1_BORDER_O2
#ifdef PARABOLIC_IBVP_H_D1V1_FX
            double fx = f(sn, tn1);
#else
            double fx = w1*f(sn, tn1) + w2*f(sn, tn0);
#endif
            ax[e]  = b21;
            bx[e]  = b23;
            cx[e]  = 0.0;
            dx[e]  = b24 * u0[e-1] + b26 * u0[e];
            dx[e] += b29 * value;
            dx[e] += ht*fx;
#else
            ax[e]  = -1.0;
            bx[e]  = +1.0;
            cx[e]  = 0.0;
            dx[e]  = ht*value;
#endif
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Robin)
        {
            const double alpha = condition.alpha();
            const double beta  = condition.beta();

            e = N;
#ifdef PARABOLIC_IBVP_H_D1V1_BORDER_O2
#ifdef PARABOLIC_IBVP_H_D1V1_FX
            double fx = f(sn, tn1);
#else
            double fx = w1*f(sn, tn1) + w2*f(sn, tn0);
#endif
            ax[e]  = beta*b21;
            bx[e]  = alpha*b22 + beta*b23;
            cx[e]  = 0.0;
            dx[e]  = beta*b24*u0[e-1] + alpha*b25*u0[e] + beta*b26*u0[e];
            dx[e]  += b29*value;
            dx[e]  += beta*ht*fx;
#else
            ax[e]  = -beta;
            bx[e]  = alpha*ht + beta;
            cx[e]  = 0.0;
            dx[e]  = ht*value;
#endif
        }

        tomasAlgorithmLeft2Right(ax+s, bx+s, cx+s, dx+s, rx+s, static_cast<unsigned int>(e-s+1));
        memcpy(u1+s, rx+s, sizeof(double)*(e-s+1));
        //for (unsigned int n=s; n<=e; n++) u1[n] = rx[n];

        layerInfo(DoubleVector(u1, N+1), tn1);

        double *tmp_u = u0; u0 = u1; u1 = tmp_u;
    }

    free(u0);
    free(u1);

    free(rx);
    free(dx);
    free(cx);
    free(bx);
    free(ax);

}

void ILoadedHeatEquationIBVP::explicit_calculate_D2V1() const {}

void ILoadedHeatEquationIBVP::implicit_calculate_D2V1() const
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

    // left border condition parameters x = xmin
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

    double **u00 = static_cast<double**>(malloc(sizeof(double*)*(M+1)));
    double **u05 = static_cast<double**>(malloc(sizeof(double*)*(M+1)));
    double **u10 = static_cast<double**>(malloc(sizeof(double*)*(M+1)));
    for (unsigned int i=0; i<=M; i++)
    {
        u00[i] = static_cast<double*>(malloc(sizeof(double)*(N+1)));
        u05[i] = static_cast<double*>(malloc(sizeof(double)*(N+1)));
        u10[i] = static_cast<double*>(malloc(sizeof(double)*(N+1)));
    }

    /***********************************************************************************************/
    /***********************************************************************************************/

    TimeNodePDE tn00, tn10;
    SpaceNodePDE sn;
    unsigned int i=0, j=0, s=0, e=N;
    int n = 0, m = 0;

    for (m=ymin, sn.j=m, sn.y=m*hy, j=0; m<=ymax; ++m, sn.j=m, sn.y=m*hy, ++j)
    {
        for (n=xmin, sn.i=n, sn.x=n*hx, i=0; n<=xmax; ++n, sn.i=n, sn.x=n*hx, ++i)
        {
            u00[j][i] = initial(sn, InitialCondition::InitialValue);
        }
    }

    tn00.i = 0; tn00.t = tn00.i*ht;
    layerInfo(DoubleMatrix(u00, M+1, N+1), tn00);

    /***********************************************************************************************/
    BoundaryConditionPDE condition; double value/*, alpha, beta, gamma*/;

    for (unsigned int ln=1; ln<=L; ln++)
    {
        tn00.i = ln-1; tn00.t = tn00.i*ht;
        tn10.i = ln-0; tn10.t = tn10.i*ht;

        /**************************************************************************************************************************/

        const size_t lps = loadedPoints().size();
        size_t* minX = new size_t[lps];
        size_t* maxX = new size_t[lps];
        size_t* minY = new size_t[lps];
        size_t* maxY = new size_t[lps];

        std::vector<size_t> loaded_indecies_x;
        std::vector<size_t> loaded_indecies_y;
        const size_t k = 0;
        for (size_t i=0; i<lps; i++)
        {
            const SpacePoint &p = _loadedPoints[i];
            const size_t rx = static_cast<size_t>(round(p.x/hx));
            const size_t ry = static_cast<size_t>(round(p.y/hy));

            minX[i] = rx - k;
            maxX[i] = rx + k;
            minY[i] = ry - k;
            maxY[i] = ry + k;

            for (size_t n=minX[i]; n<=maxX[i]; n++) { if (std::find(loaded_indecies_x.begin(), loaded_indecies_x.end(), n) == loaded_indecies_x.end()) loaded_indecies_x.push_back(n); }
            for (size_t m=minY[i]; m<=maxY[i]; m++) { if (std::find(loaded_indecies_y.begin(), loaded_indecies_y.end(), m) == loaded_indecies_y.end()) loaded_indecies_y.push_back(m); }
        }
        sort(loaded_indecies_x.begin(), loaded_indecies_x.end());
        sort(loaded_indecies_y.begin(), loaded_indecies_y.end());

        size_t row_size = loaded_indecies_y.size();
        double **w1= static_cast<double**> ( malloc(sizeof(double*)*row_size) );
        for (unsigned int row=0; row < row_size; row++) w1[row] = static_cast<double*> ( malloc(sizeof(double)*row_size) );

        //for (size_t i=0; i<loaded_indecies_x.size(); i++) printf(">> x: %zu\n", loaded_indecies_x[i]); puts("");
        //for (size_t i=0; i<loaded_indecies_y.size(); i++) printf(">> y: %zu\n", loaded_indecies_y[i]); puts("");
        //printf("%zu %d %d\n", loadedPoints().size(), loaded_indecies_x.size(), loaded_indecies_y.size());

        /**************************************************************************************************************************/

        const SpacePoint sigma(hx, hy);
        double loadedPart = 0.0;
        for (size_t i=0; i<lps; i++)
        {
            const LoadedSpacePoint &lp = _loadedPoints[i];

            SpacePoint p;
            double uv = 0.0;
            for (size_t m=minY[i]; m<=maxY[i]; m++)
            {
                p.y = m*hy;
                for (size_t n=minX[i]; n<=maxX[i]; n++)
                {
                    p.x = n*hx;
                    double w = DeltaFunction::gaussian(p, lp, sigma);
                    uv += w;// * u00[m][n];
                }
            }
            uv *= (hx*hy);
            loadedPart += lp.d * uv;
        }

        double loadedPart00 = (0.2+0.2+tn00.t) + (0.8+0.8+tn00.t);
        double loadedPart10 = (0.2+0.2+tn10.t) + (0.8+0.8+tn10.t);

        std::cout << loadedPart << " " << loadedPart00 << " " << loadedPart10 << std::endl;

        delete [] minX;
        delete [] maxX;
        delete [] minY;
        delete [] maxY;

        double *aa = static_cast<double*>(malloc(sizeof(double)*((M+1)*row_size)));
        double *bb = static_cast<double*>(malloc(sizeof(double)*((M+1)*row_size)));
        double *cc = static_cast<double*>(malloc(sizeof(double)*((M+1)*row_size)));
        double *dd = static_cast<double*>(malloc(sizeof(double)*((M+1)*row_size)));
        double *rr = static_cast<double*>(malloc(sizeof(double)*((M+1)*row_size)));

        for (size_t s=0; s<row_size; s++)
        {
            size_t offset = s*(M+1);
            for (size_t m=1; m<M; m++)
            {
                aa[offset+m] = k21;
                bb[offset+m] = k22;
                cc[offset+m] = k23;
            }
        }

        loaded_indecies_x.clear();
        loaded_indecies_y.clear();

        /**************************************************************************************************************************/

        /**************************************************** x direction apprx ***************************************************/

        for (m=ymin+1, sn.j=m, sn.y=m*hy, j=1; m<=ymax-1; ++m, sn.j=m, sn.y=m*hy, ++j)
        {
            for (n=xmin+1, sn.i=n, sn.x=n*hx, i=1; n<=xmax-1; ++n, sn.i=n, sn.x=n*hx, ++i)
            {
#ifdef PARABOLIC_IBVP_H_D2V1_FX_X
                double fx = f(sn, tn00);
#else
                double fx = 0.5 * ( f(sn, tn00) + f(sn, tn10) );
#endif
                dx[i] = k14*u00[j-1][i] + k15*u00[j][i] + k16*u00[j+1][i] + ht05*fx;// + ht05*loadedPart00;
            }

            sn.i = xmin; sn.x = xmin*hx;
#ifdef PARABOLIC_IBVP_H_D2V1_BR_X
            value = boundary(sn, tn00, condition);
#else
            value = 0.5 * ( boundary(sn, tn00, condition) + boundary(sn, tn10, condition) );
#endif
            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                s = 1;
                u05[j][0] = value;
                dx[1] -= k11 * u05[j][0];
                ax[1] = ax[0] = bx[0] = cx[0] = dx[0] = rx[0] = 0.0;
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
                s = 0;
#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
#ifdef PARABOLIC_IBVP_H_D2V1_FX_X
                double fx = f(sn, tn00);
#else
                double fx = 0.5 * ( f(sn, tn00) + f(sn, tn10) );
#endif
                ax[0]  = 0.0;
                bx[0]  = b111;
                cx[0]  = b113;
                dx[0]  = b114 * u00[j][0];
                dx[0] += b115 * ((u00[j+1][0]-2.0*u00[j][0]+u00[j-1][0])/(hy*hy));
                dx[0] += b116 * ((u00[j+1][0]-u00[j-1][0])/(2.0*hy));
                dx[0] += ht05 * fx;
                dx[0] += b117 * value;
#else
                ax[0]  = +0.0;
                bx[0]  = +1.0;
                cx[0]  = -1.0;
                dx[0]  = +hx * value;
#endif
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                const double alpha = condition.alpha();
                const double beta  = condition.beta();

                s = 0;
#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
#ifdef PARABOLIC_IBVP_H_D2V1_FX_X
                double fx = f(sn, tn00);
#else
                double fx = 0.5 * ( f(sn, tn00) + f(sn, tn10) );
#endif
                ax[0]  = 0.0;
                bx[0]  = beta * b111 + alpha * b112;
                cx[0]  = beta * b113;
                dx[0]  = beta * b114 * u00[j][0];
                dx[0] += beta * b115 * ((u00[j+1][0]-2.0*u00[j][0]+u00[j-1][0])/(hy*hy));
                dx[0] += beta * b116 * ((u00[j+1][0]-u00[j-1][0])/(2.0*hy));
                dx[0] += beta * ht05 * fx;
                dx[0] += b117 * value;
#else
                ax[0]  = +0.0;
                bx[0]  = +beta + alpha * hx;
                cx[0]  = -beta;
                dx[0]  = +hx * value;
#endif
            }

            sn.i = xmax; sn.x = xmax*hx;
#ifdef PARABOLIC_IBVP_H_D2V1_BR_X
            value = boundary(sn, tn00, condition);
#else
            value = 0.5 * ( boundary(sn, tn00, condition) + boundary(sn, tn10, condition) );
#endif
            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                e = N-1;
                u05[j][N] = value;
                dx[e] -= k13 * u05[j][N];
                cx[e] = ax[N] = bx[N] = cx[N] = dx[N] = rx[N] = 0.0;
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
                e = N;
#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
#ifdef PARABOLIC_IBVP_H_D2V1_FX_X
                double fx = f(sn, tn00);
#else
                double fx = 0.5 * ( f(sn, tn00) + f(sn, tn10) );
#endif
                ax[e]  = b121;
                bx[e]  = b123;
                cx[e]  = 0.0;
                dx[e]  = b124 * u00[j][e];
                dx[e] += b125 * ((u00[j+1][e]-2.0*u00[j][e]+u00[j-1][e])/(hy*hy));
                dx[e] += b126 * ((u00[j+1][e]-u00[j-1][e])/(2.0*hy));
                dx[e] += ht05 * fx;
                dx[e] += b127 * value;
#else
                ax[e]  = -1.0;
                bx[e]  = +1.0;
                cx[e]  = +0.0;
                dx[e]  = hx * value;
#endif
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                const double alpha = condition.alpha();
                const double beta  = condition.beta();

                e = N;
#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
#ifdef PARABOLIC_IBVP_H_D2V1_FX_X
                double fx = f(sn, tn00);
#else
                double fx = 0.5 * ( f(sn, tn00) + f(sn, tn10) );
#endif
                ax[e]  = beta * b121;
                bx[e]  = beta * b123 + alpha * b122;
                cx[e]  = 0.0;
                dx[e]  = beta * b124 * u00[j][e];
                dx[e] += beta * b125 * ((u00[j+1][e]-2.0*u00[j][e]+u00[j-1][e])/(hy*hy));
                dx[e] += beta * b126 * ((u00[j+1][e]-u00[j-1][e])/(2.0*hy));
                dx[e] += beta * ht05 * fx;
                dx[e] += b127 * value;
#else
                ax[e]  = -beta;
                bx[e]  = +beta + alpha * hx;
                cx[e]  = +0.0;
                dx[e]  = +hx * value;
#endif
            }

            tomasAlgorithmLeft2Right(ax+s, bx+s, cx+s, dx+s, rx+s, e-s+1);
            //for (size_t i=s; i<=e; i++) u05[j][i] = rx[i];
            memcpy(u05[j]+s, rx+s, sizeof (double*)*(e-s+1));
        }

        for (n=xmin, sn.i=n, sn.x=n*hx, i=0; n<=xmax; ++n, sn.i=n, sn.x=n*hx, ++i)
        {
            sn.j = ymin; sn.y = ymin*hy;
#ifdef PARABOLIC_IBVP_H_D2V1_BR_X
            value = boundary(sn, tn00, condition);
#else
            value = 0.5 * ( boundary(sn, tn00, condition) + boundary(sn, tn10, condition) );
#endif
            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                u05[0][i] = value;
            }
            if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
                u05[0][i] = (3.5*u05[1][i] - 2.0*u05[2][i] + 0.5*u05[3][i] + hy*value)/(2.0);
#else
                u05[0][i] = u05[1][i] + hy*value;
#endif
            }
            if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                const double alpha = condition.alpha();
                const double beta  = condition.beta();

#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
                //u05[0][i] = (3.5*u05[1][i] - 2.0*u05[2][i] + 0.5*u05[3][i] + hy*(1.0/beta)*value)/(2.0 + (alpha/beta)*hy);
                u05[0][i] = (4.0*beta*u05[1][i] - beta*u05[2][i] + 2.0*hy*value)/(2.0*alpha*hy+3.0*beta);
#else
                u05[0][i] = (beta*u05[1][i] + hy*value)/(beta + alpha*hy);
#endif
            }

            sn.j = ymax; sn.y = ymax*hy;
#ifdef PARABOLIC_IBVP_H_D2V1_BR_X
            value = boundary(sn, tn00, condition);
#else
            value = 0.5 * ( boundary(sn, tn00, condition) + boundary(sn, tn10, condition) );
#endif
            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                u05[M][i] = value;
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
                u05[M][i] = (3.5*u05[M-1][i] - 2.0*u05[M-2][i] + 0.5*u05[M-3][i] + hy*value)/(2.0);
#else
                u05[M][i] = u05[M-1][i] + hy*value;
#endif
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                const double alpha = condition.alpha();
                const double beta  = condition.beta();

#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
                //u05[M][i] = (3.5*u05[M-1][i] - 2.0*u05[M-2][i] + 0.5*u05[M-3][i] + hy*(1.0/beta)*value)/(2.0 + (alpha/beta)*hy);
                u05[M][i] = (4.0*beta*u05[M-1][i] - beta*u05[M-2][i] + 2.0*hy*value)/(2.0*alpha*hy+3.0*beta);
#else
                u05[M][i] = (beta*u05[M-1][i] + hy*value)/(beta + alpha*hy);
#endif
            }
        }

        /**************************************************** x direction apprx ***************************************************/


        /**************************************************** y direction apprx ***************************************************/

        for (n=xmin+1, sn.i=n, sn.x=n*hx, i=1; n<=xmax-1; ++n, sn.i=n, sn.x=n*hx, ++i)
        {
            for (m=ymin+1, sn.j=m, sn.y=m*hy, j=1; m<=ymax-1; ++m, sn.j=m, sn.y=m*hy, ++j)
            {
#ifdef PARABOLIC_IBVP_H_D2V1_FX_Y
                double fx = f(sn, tn10);
#else
                double fx = 0.5 * ( f(sn, tn00)+f(sn, tn10) );
#endif
                dy[j] = k24*u05[j][i-1] + k25*u05[j][i] + k26*u05[j][i+1] + ht05*fx + ht*loadedPart;
            }

            sn.j = ymin; sn.y = ymin*hy;
#ifdef PARABOLIC_IBVP_H_D2V1_BR_Y
            value = boundary(sn, tn10, condition);
#else
            value = 0.5 * ( boundary(sn, tn00, condition) + boundary(sn, tn10, condition) );
#endif
            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                s = 1;
                u10[0][i] = value;
                dy[1] -= k21 * u10[0][i];
                ay[1] = ay[0] = by[0] = cy[0] = dy[0] = ry[0] = 0.0;
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
                s = 0;
#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
#ifdef PARABOLIC_IBVP_H_D2V1_FX_Y
                double fx = f(sn, tn10);
#else
                double fx = 0.5 * ( f(sn, tn00)+f(sn, tn10) );
#endif
                ay[0]  = 0.0;
                by[0]  = b211;
                cy[0]  = b213;
                dy[0]  = b214 * u05[0][i];
                dy[0] += b215 * ((u05[0][i+1]-2.0*u05[0][i]+u05[0][i-1])/(hx*hx));
                dy[0] += b216 * ((u05[0][i+1]-u05[0][i-1])/(2.0*hx));
                dy[0] += ht05 * fx;
                dy[0] += b217 * value;
#else
                ay[0]  = +0.0;
                by[0]  = +1.0;
                cy[0]  = -1.0;
                dy[0]  = +hx * value;
#endif
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                const double alpha = condition.alpha();
                const double beta  = condition.beta();

                s = 0;
#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
#ifdef PARABOLIC_IBVP_H_D2V1_BR_Y
                double fx = f(sn, tn10);
#else
                double fx = 0.5 * ( f(sn, tn00)+f(sn, tn10) );
#endif
                ay[s]  = 0.0;
                by[s]  = beta * b211 + alpha * b212;
                cy[s]  = beta * b213;
                dy[s]  = beta * b214 * u05[s][i];
                dy[s] += beta * b215 * ((u05[s][i+1]-2.0*u05[s][i]+u05[s][i-1])/(hx*hx));
                dy[s] += beta * b216 * ((u05[s][i+1]-u05[s][i-1])/(2.0*hx));
                dy[s] += beta * ht05 * fx;
                dy[s] += b217 * value;
#else
                ay[s]  = +0.0;
                by[s]  = +beta + alpha * hx;
                cy[s]  = -beta;
                dy[s]  = +hx * value;
#endif
            }

            sn.j = ymax; sn.y = ymax*hy;
#ifdef PARABOLIC_IBVP_H_D2V1_BR_Y
            value = boundary(sn, tn10, condition);
#else
            value = 0.5 * ( boundary(sn, tn00, condition) + boundary(sn, tn10, condition) );
#endif
            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                e = M-1;
                u10[M][i] = value;
                dy[e] -= k23 * u10[M][i];
                cy[e] = ay[M] = by[M] = cy[M] = dy[M] = ry[M] = 0.0;
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
                e = M;
#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
#ifdef PARABOLIC_IBVP_H_D2V1_FX_Y
                double fx = f(sn, tn10);
#else
                double fx = 0.5 * ( f(sn, tn00)+f(sn, tn10) );
#endif
                ay[e]  = b221;
                by[e]  = b223;
                cy[e]  = 0.0;
                dy[e]  = b224 * u05[e][i];
                dy[e] += b225 * ((u05[e][i+1]-2.0*u05[e][i]+u05[e][i-1])/(hx*hx));
                dy[e] += b226 * ((u05[e][i+1]-u05[e][i-1])/(2.0*hx));
                dy[e] += ht05 * fx;
                dy[e] += b227 * value;
#else
                ay[e]  = -1.0;
                by[e]  = +1.0;
                cy[e]  = +0.0;
                dy[e]  = +hx * value;
#endif
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                const double alpha = condition.alpha();
                const double beta  = condition.beta();

                e = M;
#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
#ifdef PARABOLIC_IBVP_H_D2V1_FX_Y
                double fx = f(sn, tn10);
#else
                double fx = 0.5 * ( f(sn, tn00)+f(sn, tn10) );
#endif
                ay[e]  = beta * b221;
                by[e]  = beta * b223 + alpha * b222;
                cy[e]  = 0.0;
                dy[e]  = beta * b224 * u05[e][i];
                dy[e] += beta * b225 * ((u05[e][i+1]-2.0*u05[e][i]+u05[e][i-1])/(hx*hx));
                dy[e] += beta * b226 * ((u05[e][i+1]-u05[e][i-1])/(2.0*hx));
                dy[e] += beta * ht05 * fx;
                dy[e] += b227 * value;
#else
                ay[e]  = -beta;
                by[e]  = +beta + alpha * hx;
                cy[e]  = +0.0;
                dy[e]  = +hx * value;
#endif
            }

            tomasAlgorithmLeft2Right(ay+s, by+s, cy+s, dy+s, ry+s, e-s+1);
            for (unsigned int j=s; j<=e; j++) u10[j][i] = ry[j];
            //memcpy(u10[i]+s, ry+s, sizeof(double*)*(e-s+1));
        }

        for (m=ymin, sn.j=m, sn.y=m*hy, j=0; m<=ymax; ++m, sn.j=m, sn.y=m*hy, ++j)
        {
            sn.i = xmin; sn.x = xmin*hx;
            sn.i = xmin; sn.x = xmin*hx;
#ifdef PARABOLIC_IBVP_H_D2V1_BR_Y
            value = boundary(sn, tn10, condition);
#else
            value = 0.5 * ( boundary(sn, tn00, condition) + boundary(sn, tn10, condition) );
#endif
            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                u10[j][0] = value;
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
                u10[j][0] = (3.5*u10[j][1] - 2.0*u10[j][2] + 0.5*u10[j][3] + hx*value)/(2.0);
#else
                u10[j][0] = u10[j][1] + hx*value;
#endif
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                const double alpha = condition.alpha();
                const double beta  = condition.beta();

#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
                //u10[j][0] = (3.5*u10[j][1] - 2.0*u10[j][2] + 0.5*u10[j][3] + hx*(1.0/beta)*value)/(2.0 + (alpha/beta)*hx);
                u10[j][0] = (4.0*beta*u10[j][1] - beta*u10[j][2] + 2.0*hx*value)/(2.0*alpha*hx+3.0*beta);
#else
                u10[j][0] = (beta*u10[j][1] + hx*value)/(beta + alpha*hx);
#endif
            }

            sn.i = xmax; sn.x = xmax*hx;
#ifdef PARABOLIC_IBVP_H_D2V1_BR_Y
            value = boundary(sn, tn10, condition);
#else
            value = 0.5 * ( boundary(sn, tn00, condition) + boundary(sn, tn10, condition) );
#endif
            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                u10[j][N] = value;
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
                u10[j][N] = (3.5*u10[j][N-1] - 2.0*u10[j][N-2] + 0.5*u10[j][N-3] + hx*value)/(2.0);
#else
                u10[j][N] = u10[j][N-1] + hx*value;
#endif
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                const double alpha = condition.alpha();
                const double beta  = condition.beta();

#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
                //u10[j][N] = (3.5*u10[j][N-1] - 2.0*u10[j][N-2] + 0.5*u10[j][N-3] + hx*(1.0/beta)*value)/(2.0 + (alpha/beta)*hx);
                u10[j][N] = (4.0*beta*u10[j][N-1] - beta*u10[j][N-2] + 2.0*hx*value)/(2.0*alpha*hx+3.0*beta);
#else
                u10[j][N] = (beta*u10[j][N-1] + hx*value)/(beta + alpha*hx);
#endif
            }
        }

        layerInfo(DoubleMatrix(u10, M+1, N+1), tn10);

        /**************************************************** y direction apprx ***************************************************/
        double **_tmp = u00; u00 = u10; u10 = _tmp;
    }

    for (unsigned int i=0; i<=M; i++) { free(u00[i]); free(u05[i]); free(u10[i]); }
    free(u00); free(u05); free(u10);

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

void ILoadedHeatEquationIBVP::setLoadedPoints(const std::vector<LoadedSpacePoint> &loadedPoints)
{
    this->_loadedPoints = loadedPoints;
}

void ILoadedHeatEquationIBVP::setLoadedPoints0(const std::vector<LoadedSpacePoint> &loadedPoints)
{
    this->_loadedPoints = loadedPoints;
}

void ILoadedHeatEquationIBVP::setLoadedPoints1(const std::vector<LoadedSpacePoint> &loadedPoints)
{
    this->_loadedPoints = loadedPoints;
}

const std::vector<LoadedSpacePoint> ILoadedHeatEquationIBVP::loadedPoints() const
{
    return _loadedPoints;
}

//------------------------------------------------------------------------------------------------------------------------------------------------------//

ILoadedHeatEquationFBVP::ILoadedHeatEquationFBVP(double thermalDiffusivity, double thermalConductivity, double thermalConvection) :
    IHeatEquationFBVP(thermalDiffusivity, thermalConductivity, thermalConvection) {}

ILoadedHeatEquationFBVP::~ILoadedHeatEquationFBVP() {}

ILoadedHeatEquationFBVP::ILoadedHeatEquationFBVP(const ILoadedHeatEquationFBVP &) {}

ILoadedHeatEquationFBVP& ILoadedHeatEquationFBVP::operator=(const ILoadedHeatEquationFBVP &)
{
    return *this;
}

void ILoadedHeatEquationFBVP::explicit_calculate_D1V1() const {}

void ILoadedHeatEquationFBVP::implicit_calculate_D1V1() const {}

void ILoadedHeatEquationFBVP::explicit_calculate_D2V1() const {}

void ILoadedHeatEquationFBVP::implicit_calculate_D2V1() const
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

    double **p00 = static_cast<double**>(malloc(sizeof(double*)*(M+1)));
    double **p05 = static_cast<double**>(malloc(sizeof(double*)*(M+1)));
    double **p10 = static_cast<double**>(malloc(sizeof(double*)*(M+1)));
    for (unsigned int i=0; i<=M; i++)
    {
        p00[i] = static_cast<double*>(malloc(sizeof(double)*(N+1)));
        p05[i] = static_cast<double*>(malloc(sizeof(double)*(N+1)));
        p10[i] = static_cast<double*>(malloc(sizeof(double)*(N+1)));
    }

    /***********************************************************************************************/
    /***********************************************************************************************/

    TimeNodePDE tn00, tn10;
    SpaceNodePDE sn;
    unsigned int i=0, j=0, s=0, e=N;
    int n = 0, m = 0;

    for (m=ymin, sn.j=m, sn.y=m*hy, j=0; m<=ymax; ++m, sn.j=m, sn.y=m*hy, ++j)
    {
        for (n=xmin, sn.i=n, sn.x=n*hx, i=0; n<=xmax; ++n, sn.i=n, sn.x=n*hx, ++i)
        {
            p00[j][i] = final(sn, FinalCondition::FinalValue);
        }
    }

    tn00.i = L; tn00.t = tn00.i*ht;
    layerInfo(DoubleMatrix(p00, M+1, N+1), tn00);

    /***********************************************************************************************/
    BoundaryConditionPDE condition; double value/*, alpha, beta, gamma*/;

    for (unsigned int ln=L-1, size_ln = static_cast<unsigned int>(0)-1; ln != size_ln; ln--)
    {
        tn00.i = ln+1; tn00.t = tn00.i*ht;
        tn10.i = ln-0; tn10.t = tn10.i*ht;

        /**************************************************************************************************************************/

        const size_t lps = loadedPoints().size();
        size_t* minX = new size_t[lps];
        size_t* maxX = new size_t[lps];
        size_t* minY = new size_t[lps];
        size_t* maxY = new size_t[lps];

        std::vector<size_t> loaded_indecies_x;
        std::vector<size_t> loaded_indecies_y;
        const size_t k = 4;
        for (size_t i=0; i<lps; i++)
        {
            const SpacePoint &p = _loadedPoints[i];
            const size_t rx = static_cast<size_t>(round(p.x/hx));
            const size_t ry = static_cast<size_t>(round(p.y/hy));

            minX[i] = rx - k;
            maxX[i] = rx + k;
            minY[i] = ry - k;
            maxY[i] = ry + k;

            for (size_t n=minX[i]; n<=maxX[i]; n++) { if (std::find(loaded_indecies_x.begin(), loaded_indecies_x.end(), n) == loaded_indecies_x.end()) loaded_indecies_x.push_back(n); }
            for (size_t m=minY[i]; m<=maxY[i]; m++) { if (std::find(loaded_indecies_y.begin(), loaded_indecies_y.end(), m) == loaded_indecies_y.end()) loaded_indecies_y.push_back(m); }
        }
        sort(loaded_indecies_x.begin(), loaded_indecies_x.end());
        sort(loaded_indecies_y.begin(), loaded_indecies_y.end());

        size_t row_size = loaded_indecies_y.size();
        double **w1= static_cast<double**> ( malloc(sizeof(double*)*row_size) );
        for (unsigned int row=0; row < row_size; row++) w1[row] = static_cast<double*> ( malloc(sizeof(double)*row_size) );

        for (size_t i=0; i<loaded_indecies_x.size(); i++) printf(">> x: %zu\n", loaded_indecies_x[i]); puts("");
        for (size_t i=0; i<loaded_indecies_y.size(); i++) printf(">> y: %zu\n", loaded_indecies_y[i]); puts("");
        printf("%zu %d %d\n", loadedPoints().size(), loaded_indecies_x.size(), loaded_indecies_y.size());

        /**************************************************************************************************************************/

        const SpacePoint sigma(hx, hy);
        double loadedPart = 0.0;
        for (size_t i=0; i<lps; i++)
        {
            const LoadedSpacePoint &lp = _loadedPoints[i];

            SpacePoint p;
            double uv = 0.0;
            for (size_t m=minY[i]; m<=maxY[i]; m++)
            {
                p.y = m*hy;
                for (size_t n=minX[i]; n<=maxX[i]; n++)
                {
                    p.x = n*hx;
                    double w = DeltaFunction::gaussian(p, lp, sigma);
                    uv += w * p00[m][n];
                }
            }
            uv *= (hx*hy);
            loadedPart += lp.d * uv;
        }

        delete [] minX;
        delete [] maxX;
        delete [] minY;
        delete [] maxY;

        loaded_indecies_x.clear();
        loaded_indecies_y.clear();

        /**************************************************************************************************************************/

        /**************************************************** x direction apprx ***************************************************/

        for (m=ymin+1, sn.j=m, sn.y=m*hy, j=1; m<=ymax-1; ++m, sn.j=m, sn.y=m*hy, ++j)
        {
            for (n=xmin+1, sn.i=n, sn.x=n*hx, i=1; n<=xmax-1; ++n, sn.i=n, sn.x=n*hx, ++i)
            {
                //double fx = f(sn, tn00)
                double fx = 0.5 * ( f(sn, tn00) + f(sn, tn10) );
                dx[i] = k14*p00[j-1][i] + k15*p00[j][i] + k16*p00[j+1][i] + ht05*fx + ht05*loadedPart;
            }

            sn.i = xmin; sn.x = xmin*hx;
            value = 0.5*( boundary(sn, tn00, condition) + boundary(sn, tn10, condition) );

            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                s = 1;
                p05[j][0] = /*(gamma/alpha)**/value;
                dx[1] -= k11 * p05[j][0];
                ax[1] = ax[0] = bx[0] = cx[0] = dx[0] = rx[0] = 0.0;
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
                s = 0;
#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
                ax[s]  = 0.0;
                bx[s]  = b111;
                cx[s]  = b113;
                dx[s]  = b114 * p00[j][s];
                dx[s] += b115 * ((p00[j+1][s]-2.0*p00[j][s]+p00[j-1][s])/(hy*hy));
                dx[s] += b116 * ((p00[j+1][s]-p00[j-1][s])/(2.0*hy));
                dx[s] += /*beta  **/ ht05 * f(sn, tn00);
                dx[s] += /*gamma **/ b117 * value;
#else
                ax[s]  = 0.0;
                bx[s]  = +/*beta*/1.0 + /*alpha **/ hx;
                cx[s]  = -/*beta*/1.0;
                dx[s]  = +/*gamma **/ hx * value;
#endif
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                const double alpha = condition.alpha();
                const double beta  = condition.beta();
                //const double gamma = condition.gamma();
                const double gamma = 1.0;

                s = 0;
#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
                ax[s]  = 0.0;
                bx[s]  = beta  * b111 + alpha * b112;
                cx[s]  = beta  * b113;
                dx[s]  = beta  * b114 * p00[j][s];
                dx[s] += beta  * b115 * ((p00[j+1][s]-2.0*p00[j][s]+p00[j-1][s])/(hy*hy));
                dx[s] += beta  * b116 * ((p00[j+1][s]-p00[j-1][s])/(2.0*hy));
                dx[s] += beta  * ht05 * f(sn, tn00);
                dx[s] += gamma * b117 * value;
#else
                ax[s]  = 0.0;
                bx[s]  = +beta + alpha * hx;
                cx[s]  = -beta;
                dx[s]  = +gamma * hx * value;
#endif
            }

            sn.i = xmax; sn.x = xmax*hx;
            value = 0.5 * ( boundary(sn, tn00, condition) + boundary(sn, tn10, condition) );

            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                e = N-1;
                p05[j][N] = /*(gamma/alpha)**/value;
                dx[N-1] -= k13 * p05[j][N];
                cx[N-1] = ax[N] = bx[N] = cx[N] = dx[N] = rx[N] = 0.0;
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
                e = N;
#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
                ax[e]  = /*beta  **/ b121;
                bx[e]  = /*beta  **/ b123;
                cx[e]  = 0.0;
                dx[e]  = /*beta  **/ b124 * p00[j][e];
                dx[e] += /*beta  **/ b125 * ((p00[j+1][e]-2.0*p00[j][e]+p00[j-1][e])/(hy*hy));
                dx[e] += /*beta  **/ b126 * ((p00[j+1][e]-p00[j-1][e])/(2.0*hy));
                dx[e] += /*beta  **/ ht05 * f(sn, tn00);
                dx[e] += /*gamma **/ b127 * value;
#else
                ax[e]  = -/*beta*/1.0;
                bx[e]  = +/*beta*/1.0 + /*alpha **/ hx;
                cx[e]  = 0.0;
                dx[e]  = /*gamma **/ hx * value;
#endif
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                const double alpha = condition.alpha();
                const double beta  = condition.beta();
                //const double gamma = condition.gamma();
                const double gamma = 1.0;

                e = N;
#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
                ax[e]  = beta  * b121;
                bx[e]  = beta  * b123 + alpha * b122;
                cx[e]  = 0.0;
                dx[e]  = beta  * b124 * p00[j][e];
                dx[e] += beta  * b125 * ((p00[j+1][e]-2.0*p00[j][e]+p00[j-1][e])/(hy*hy));
                dx[e] += beta  * b126 * ((p00[j+1][e]-p00[j-1][e])/(2.0*hy));
                dx[e] += beta  * ht05 * f(sn, tn00);
                dx[e] += gamma * b127 * value;
#else
                ax[e]  = -beta;
                bx[e]  = +beta + alpha * hx;
                cx[e]  = 0.0;
                dx[e]  = gamma * hx * value;
#endif
            }

            tomasAlgorithmLeft2Right(ax+s, bx+s, cx+s, dx+s, rx+s, e-s+1);
            //for (unsigned int i=s; i<=e; i++) p05[j][i] = rx[i];
            memcpy(p05[j]+s, rx+s, sizeof (double*)*(e-s+1));
        }

        for (n=xmin, sn.i=n, sn.x=n*hx, i=0; n<=xmax; ++n, sn.i=n, sn.x=n*hx, ++i)
        {
            sn.j = ymin; sn.y = ymin*hy;
            value = 0.5 * ( boundary(sn, tn00, condition) + boundary(sn, tn10, condition) );

            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                p05[0][i] = /*(gamma/alpha)**/value;
            }
            if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
                p05[0][i] = (3.5*p05[1][i] - 2.0*p05[2][i] + 0.5*p05[3][i] + hy/**(gamma/beta)*/*value)/(2.0);
#else
                p05[0][i] = p05[1][i] + hy/**(gamma/beta)*/*value;
#endif

            }
            if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                const double alpha = condition.alpha();
                const double beta  = condition.beta();
                //const double gamma = condition.gamma();
                const double gamma = 1.0;

#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
                p05[0][i] = (3.5*p05[1][i] - 2.0*p05[2][i] + 0.5*p05[3][i] + hy*(gamma/beta)*value)/(2.0 + (alpha/beta)*hy);
#else
                p05[0][i] = (p05[1][i] + hy*(gamma/beta)*value)/(1.0 + hy*(alpha/beta));
#endif
            }

            sn.j = ymax; sn.y = ymax*hy;
            value = 0.5 * ( boundary(sn, tn00, condition) + boundary(sn, tn10, condition) );

            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                p05[M][i] = /*(gamma/alpha)**/value;
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
                p05[M][i] = (3.5*p05[M-1][i] - 2.0*p05[M-2][i] + 0.5*p05[M-3][i] + hy/**(gamma/beta)*/*value)/(2.0 + /*(alpha/beta)**/hy);
#else
                p05[M][i] = p05[M-1][i] + hy/**(gamma/beta)*/*value;
#endif
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                const double alpha = condition.alpha();
                const double beta  = condition.beta();
                //const double gamma = condition.gamma();
                const double gamma = 1.0;

#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
                p05[M][i] = (3.5*p05[M-1][i] - 2.0*p05[M-2][i] + 0.5*p05[M-3][i] + hy*(gamma/beta)*value)/(2.0 + (alpha/beta)*hy);
#else
                p05[M][i] = (p05[M-1][i] + hy*(gamma/beta)*value)/(1.0 + hy*(alpha/beta));
#endif
            }
        }

        /**************************************************** x direction apprx ***************************************************/

        /**************************************************** y direction apprx ***************************************************/

        for (n=xmin+1, sn.i=n, sn.x=n*hx, i=1; n<=xmax-1; ++n, sn.i=n, sn.x=n*hx, ++i)
        {
            for (m=ymin+1, sn.j=m, sn.y=m*hy, j=1; m<=ymax-1; ++m, sn.j=m, sn.y=m*hy, ++j)
            {
                //double fx = f(sn, tn10);
                double fx = 0.5 * ( f(sn, tn00)+f(sn, tn10) );
                dy[j] = k24*p05[j][i-1] + k25*p05[j][i] + k26*p05[j][i+1] + ht05*fx + ht05*loadedPart;
            }

            sn.j = ymin; sn.y = ymin*hy;
            value = boundary(sn, tn10, condition);

            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                s = 1;
                p10[0][i] = /*(gamma/alpha)**/value;
                dy[1] -= k21 * p10[0][i];
                ay[1] = ay[0] = by[0] = cy[0] = dy[0] = ry[0] = 0.0;
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
                s = 0;
#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
                ay[s]  = 0.0;
                by[s]  = /*beta  **/ b211;
                cy[s]  = /*beta  **/ b213;
                dy[s]  = /*beta  **/ b214 * p05[s][i];
                dy[s] += /*beta  **/ b215 * ((p05[s][i+1]-2.0*p05[s][i]+p05[s][i-1])/(hx*hx));
                dy[s] += /*beta  **/ b216 * ((p05[s][i+1]-p05[s][i-1])/(2.0*hx));
                dy[s] += /*beta  **/ ht05 * f(sn, tn10);
                dy[s] += /*gamma **/ b217 * value;
#else
                ay[s]  = 0.0;
                by[s]  = +/*beta*/1.0 + /*alpha **/ hx;
                cy[s]  = -/*beta*/1.0;
                dy[s]  = +/*gamma **/ hx * value;
#endif
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                const double alpha = condition.alpha();
                const double beta  = condition.beta();
                //const double gamma = condition.gamma();
                const double gamma = 1.0;

                s = 0;
#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
                ay[s]  = 0.0;
                by[s]  = beta  * b211 + alpha * b212;
                cy[s]  = beta  * b213;
                dy[s]  = beta  * b214 * p05[s][i];
                dy[s] += beta  * b215 * ((p05[s][i+1]-2.0*p05[s][i]+p05[s][i-1])/(hx*hx));
                dy[s] += beta  * b216 * ((p05[s][i+1]-p05[s][i-1])/(2.0*hx));
                dy[s] += beta  * ht05 * f(sn, tn10);
                dy[s] += gamma * b217 * value;
#else
                ay[s]  = 0.0;
                by[s]  = +beta + alpha * hx;
                cy[s]  = -beta;
                dy[s]  = +gamma * hx * value;
#endif
            }

            sn.j = ymax; sn.y = ymax*hy;
            value = boundary(sn, tn10, condition);

            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                e = M-1;
                p10[M][i] = /*(gamma/alpha)**/value;
                dy[M-1] -= k23 * p10[M][i];
                cy[M-1] = ay[M] = by[M] = cy[M] = dy[M] = ry[M] = 0.0;
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
                e = M;
#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
                ay[e]  = /*beta  **/ b221;
                by[e]  = /*beta  **/ b223;
                cy[e]  = 0.0;
                dy[e]  = /*beta  **/ b224 * p05[e][i];
                dy[e] += /*beta  **/ b225 * ((p05[e][i+1]-2.0*p05[e][i]+p05[e][i-1])/(hx*hx));
                dy[e] += /*beta  **/ b226 * ((p05[e][i+1]-p05[e][i-1])/(2.0*hx));
                dy[e] += /*beta  **/ ht05 * f(sn, tn10);
                dy[e] += /*gamma **/ b227 * value;
#else
                ay[e]  = -/*beta*/1.0;
                by[e]  = +/*beta*/1.0 + /*alpha **/ hx;
                cy[e]  = 0.0;
                dy[e]  = +/*gamma **/ hx * value;
#endif
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                const double alpha = condition.alpha();
                const double beta  = condition.beta();
                //const double gamma = condition.gamma();
                const double gamma = 1.0;

                e = M;
#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
                ay[e]  = beta  * b221;
                by[e]  = beta  * b223 + alpha * b222;
                cy[e]  = 0.0;
                dy[e]  = beta  * b224 * p05[e][i];
                dy[e] += beta  * b225 * ((p05[e][i+1]-2.0*p05[e][i]+p05[e][i-1])/(hx*hx));
                dy[e] += beta  * b226 * ((p05[e][i+1]-p05[e][i-1])/(2.0*hx));
                dy[e] += beta  * ht05 * f(sn, tn10);
                dy[e] += gamma * b227 * value;
#else
                ay[e]  = -beta;
                by[e]  = +beta + alpha * hx;
                cy[e]  = 0.0;
                dy[e]  = +gamma * hx * value;
#endif
            }

            tomasAlgorithmLeft2Right(ay+s, by+s, cy+s, dy+s, ry+s, e-s+1);
            for (unsigned int j=s; j<=e; j++) p10[j][i] = ry[j];
            //memcpy(p10[i]+s, ry+s, sizeof(double*)*(e-s+1));
        }

        for (m=ymin, sn.j=m, sn.y=m*hy, j=0; m<=ymax; ++m, sn.j=m, sn.y=m*hy, ++j)
        {
            sn.i = xmin; sn.x = xmin*hx;
            value = boundary(sn, tn10, condition);

            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                p10[j][0] = /*(gamma/alpha)**/value;
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
                p10[j][0] = (3.5*p10[j][1] - 2.0*p10[j][2] + 0.5*p10[j][3] + hx/**(gamma/beta)*/*value)/(2.0);
#else
                p10[j][0] = p10[j][1] + hx/**(gamma/beta)*/*value;
#endif
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                const double alpha = condition.alpha();
                const double beta  = condition.beta();
                //const double gamma = condition.gamma();
                const double gamma = 1.0;

#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
                p10[j][0] = (3.5*p10[j][1] - 2.0*p10[j][2] + 0.5*p10[j][3] + hx*(gamma/beta)*value)/(2.0 + (alpha/beta)*hx);
#else
                p10[j][0] = (p10[j][1] + hx*(gamma/beta)*value)/(1.0 + hx*(alpha/beta));
#endif
            }

            sn.i = xmax; sn.x = xmax*hx;
            value = boundary(sn, tn10, condition);

            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                p10[j][N] = /*(gamma/alpha)**/value;
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
                p10[j][N] = (3.5*p10[j][N-1] - 2.0*p10[j][N-2] + 0.5*p10[j][N-3] + hx/**(gamma/beta)*/*value)/(2.0);
#else
                p10[j][N] = p10[j][N-1] + hx/**(gamma/beta)*/*value;
#endif
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                const double alpha = condition.alpha();
                const double beta  = condition.beta();
                //const double gamma = condition.gamma();
                const double gamma = 1.0;

#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
                p10[j][N] = (3.5*p10[j][N-1] - 2.0*p10[j][N-2] + 0.5*p10[j][N-3] + hx*(gamma/beta)*value)/(2.0 + (alpha/beta)*hx);
#else
                p10[j][N] = (p10[j][N-1] + hx*(gamma/beta)*value)/(1.0 + hx*(alpha/beta));
#endif
            }
        }

        layerInfo(DoubleMatrix(p10, M+1, N+1), tn10);

        /**************************************************** y direction apprx ***************************************************/
        double **_tmp = p00; p00 = p10; p10 = _tmp;
    }

    for (unsigned int i=0; i<=M; i++) { free(p00[i]); free(p05[i]); free(p10[i]); }
    free(p00); free(p05); free(p10);

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

void ILoadedHeatEquationFBVP::implicit_calculate_D2V2() const
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

    double **p00 = static_cast<double**>(malloc(sizeof(double*)*(M+1)));
    double **p05 = static_cast<double**>(malloc(sizeof(double*)*(M+1)));
    double **p10 = static_cast<double**>(malloc(sizeof(double*)*(M+1)));
    for (unsigned int i=0; i<=M; i++)
    {
        p00[i] = static_cast<double*>(malloc(sizeof(double)*(N+1)));
        p05[i] = static_cast<double*>(malloc(sizeof(double)*(N+1)));
        p10[i] = static_cast<double*>(malloc(sizeof(double)*(N+1)));
    }

    /***********************************************************************************************/
    /***********************************************************************************************/

    TimeNodePDE tn00, tn10;
    SpaceNodePDE sn;
    unsigned int i=0, j=0, s=0, e=N;
    int n = 0, m = 0;

    for (m=ymin, sn.j=m, sn.y=m*hy, j=0; m<=ymax; ++m, sn.j=m, sn.y=m*hy, ++j)
    {
        for (n=xmin, sn.i=n, sn.x=n*hx, i=0; n<=xmax; ++n, sn.i=n, sn.x=n*hx, ++i)
        {
            p00[j][i] = final(sn, FinalCondition::FinalValue);
        }
    }

    tn00.i = L; tn00.t = tn00.i*ht;
    layerInfo(DoubleMatrix(p00, M+1, N+1), tn00);

    /***********************************************************************************************/
    BoundaryConditionPDE condition; double value/*, alpha, beta, gamma*/;

    for (unsigned int ln=L-1, size_ln = static_cast<unsigned int>(0)-1; ln != size_ln; ln--)
    {
        tn00.i = ln+1; tn00.t = tn00.i*ht;
        tn10.i = ln-0; tn10.t = tn10.i*ht;

        /**************************************************** x direction apprx ***************************************************/

        for (m=ymin+1, sn.j=m, sn.y=m*hy, j=1; m<=ymax-1; ++m, sn.j=m, sn.y=m*hy, ++j)
        {
            for (n=xmin+1, sn.i=n, sn.x=n*hx, i=1; n<=xmax-1; ++n, sn.i=n, sn.x=n*hx, ++i)
            {
                double fx = f(sn, tn00);
                //double fx = 0.5 * ( f(sn, tn00) + f(sn, tn10) );
                dx[i] = k14*p00[j-1][i] + k15*p00[j][i] + k16*p00[j+1][i] + ht05*fx;
            }

            sn.i = xmin; sn.x = xmin*hx;
            value = 0.5*( boundary(sn, tn00, condition) + boundary(sn, tn10, condition) );

            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                s = 1;
                p05[j][0] = /*(gamma/alpha)**/value;
                dx[1] -= k11 * p05[j][0];
                ax[1] = ax[0] = bx[0] = cx[0] = dx[0] = rx[0] = 0.0;
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
                s = 0;
#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
                double fx = f(sn, tn00);
                //double fx = 0.5 * ( f(sn, tn00) + f(sn, tn10) );

                ax[s]  = 0.0;
                bx[s]  = /*beta  **/ b111;
                cx[s]  = /*beta  **/ b113;
                dx[s]  = /*beta  **/ b114 * p00[j][s];
                dx[s] += /*beta  **/ b115 * ((p00[j+1][s]-2.0*p00[j][s]+p00[j-1][s])/(hy*hy));
                dx[s] += /*beta  **/ b116 * ((p00[j+1][s]-p00[j-1][s])/(2.0*hy));
                dx[s] += /*beta  **/ ht05 * fx;
                dx[s] += /*gamma **/ b117 * value;
#else
                ax[s]  = 0.0;
                bx[s]  = +/*beta*/1.0 + /*alpha **/ hx;
                cx[s]  = -/*beta*/1.0;
                dx[s]  = +/*gamma **/ hx * value;
#endif
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                const double alpha = condition.alpha();
                const double beta  = condition.beta();
                //const double gamma = condition.gamma();
                const double gamma = 1.0;

                s = 0;
#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
                double fx = f(sn, tn00);
                //double fx = 0.5 * ( f(sn, tn00) + f(sn, tn10) );

                ax[s]  = 0.0;
                bx[s]  = beta  * b111 + alpha * b112;
                cx[s]  = beta  * b113;
                dx[s]  = beta  * b114 * p00[j][s];
                dx[s] += beta  * b115 * ((p00[j+1][s]-2.0*p00[j][s]+p00[j-1][s])/(hy*hy));
                dx[s] += beta  * b116 * ((p00[j+1][s]-p00[j-1][s])/(2.0*hy));
                dx[s] += beta  * ht05 * fx;
                dx[s] += gamma * b117 * value;
#else
                ax[s]  = 0.0;
                bx[s]  = +beta + alpha * hx;
                cx[s]  = -beta;
                dx[s]  = +gamma * hx * value;
#endif
            }

            sn.i = xmax; sn.x = xmax*hx;
            value = 0.5 * ( boundary(sn, tn00, condition) + boundary(sn, tn10, condition) );

            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                e = N-1;
                p05[j][N] = /*(gamma/alpha)**/value;
                dx[N-1] -= k13 * p05[j][N];
                cx[N-1] = ax[N] = bx[N] = cx[N] = dx[N] = rx[N] = 0.0;
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
                e = N;
#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
                double fx = f(sn, tn00);
                //double fx = 0.5 * ( f(sn, tn00) + f(sn, tn10) );

                ax[e]  = /*beta  **/ b121;
                bx[e]  = /*beta  **/ b123;
                cx[e]  = 0.0;
                dx[e]  = /*beta  **/ b124 * p00[j][e];
                dx[e] += /*beta  **/ b125 * ((p00[j+1][e]-2.0*p00[j][e]+p00[j-1][e])/(hy*hy));
                dx[e] += /*beta  **/ b126 * ((p00[j+1][e]-p00[j-1][e])/(2.0*hy));
                dx[e] += /*beta  **/ ht05 * fx;
                dx[e] += /*gamma **/ b127 * value;
#else
                ax[e]  = -/*beta*/1.0;
                bx[e]  = +/*beta*/1.0 + /*alpha **/ hx;
                cx[e]  = 0.0;
                dx[e]  = /*gamma **/ hx * value;
#endif
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                const double alpha = condition.alpha();
                const double beta  = condition.beta();
                //const double gamma = condition.gamma();
                const double gamma = 1.0;

                e = N;
#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
                double fx = f(sn, tn00);
                //double fx = 0.5 * ( f(sn, tn00) + f(sn, tn10) );

                ax[e]  = beta  * b121;
                bx[e]  = beta  * b123 + alpha * b122;
                cx[e]  = 0.0;
                dx[e]  = beta  * b124 * p00[j][e];
                dx[e] += beta  * b125 * ((p00[j+1][e]-2.0*p00[j][e]+p00[j-1][e])/(hy*hy));
                dx[e] += beta  * b126 * ((p00[j+1][e]-p00[j-1][e])/(2.0*hy));
                dx[e] += beta  * ht05 * fx;
                dx[e] += gamma * b127 * value;
#else
                ax[e]  = -beta;
                bx[e]  = +beta + alpha * hx;
                cx[e]  = 0.0;
                dx[e]  = gamma * hx * value;
#endif
            }

            tomasAlgorithmLeft2Right(ax+s, bx+s, cx+s, dx+s, rx+s, e-s+1);
            //for (unsigned int i=s; i<=e; i++) p05[j][i] = rx[i];
            memcpy(p05[j]+s, rx+s, sizeof (double*)*(e-s+1));
        }

        for (n=xmin, sn.i=n, sn.x=n*hx, i=0; n<=xmax; ++n, sn.i=n, sn.x=n*hx, ++i)
        {
            sn.j = ymin; sn.y = ymin*hy;
            value = 0.5 * ( boundary(sn, tn00, condition) + boundary(sn, tn10, condition) );

            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                p05[0][i] = /*(gamma/alpha)**/value;
            }
            if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
                p05[0][i] = (3.5*p05[1][i] - 2.0*p05[2][i] + 0.5*p05[3][i] + hy/**(gamma/beta)*/*value)/(2.0);
#else
                p05[0][i] = p05[1][i] + hy/**(gamma/beta)*/*value;
#endif
            }
            if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                const double alpha = condition.alpha();
                const double beta  = condition.beta();
                //const double gamma = condition.gamma();
                const double gamma = 1.0;

#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
                p05[0][i] = (3.5*p05[1][i] - 2.0*p05[2][i] + 0.5*p05[3][i] + hy*(gamma/beta)*value)/(2.0 + (alpha/beta)*hy);
#else
                p05[0][i] = (p05[1][i] + hy*(gamma/beta)*value)/(1.0 + hy*(alpha/beta));
#endif
            }

            sn.j = ymax; sn.y = ymax*hy;
            value = 0.5 * ( boundary(sn, tn00, condition) + boundary(sn, tn10, condition) );

            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                p05[M][i] = /*(gamma/alpha)**/value;
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
                p05[M][i] = (3.5*p05[M-1][i] - 2.0*p05[M-2][i] + 0.5*p05[M-3][i] + hy/**(gamma/beta)*/*value)/(2.0 + /*(alpha/beta)**/hy);
#else
                p05[M][i] = p05[M-1][i] + hy/**(gamma/beta)*/*value;
#endif
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                const double alpha = condition.alpha();
                const double beta  = condition.beta();
                //const double gamma = condition.gamma();
                const double gamma = 1.0;

#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
                p05[M][i] = (3.5*p05[M-1][i] - 2.0*p05[M-2][i] + 0.5*p05[M-3][i] + hy*(gamma/beta)*value)/(2.0 + (alpha/beta)*hy);
#else
                p05[M][i] = (p05[M-1][i] + hy*(gamma/beta)*value)/(1.0 + hy*(alpha/beta));
#endif
            }
        }

        /**************************************************** x direction apprx ***************************************************/

        /**************************************************** y direction apprx ***************************************************/

        for (n=xmin+1, sn.i=n, sn.x=n*hx, i=1; n<=xmax-1; ++n, sn.i=n, sn.x=n*hx, ++i)
        {
            for (m=ymin+1, sn.j=m, sn.y=m*hy, j=1; m<=ymax-1; ++m, sn.j=m, sn.y=m*hy, ++j)
            {
                double fx = f(sn, tn00);/*double fx = f(sn, tn10);*/
                //double fx = 0.5 * ( f(sn, tn00)+f(sn, tn10) );
                dy[j] = k24*p05[j][i-1] + k25*p05[j][i] + k26*p05[j][i+1] + ht05*fx;
            }

            sn.j = ymin; sn.y = ymin*hy;
            value = boundary(sn, tn10, condition);

            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                s = 1;
                p10[0][i] = /*(gamma/alpha)**/value;
                dy[1] -= k21 * p10[0][i];
                ay[1] = ay[0] = by[0] = cy[0] = dy[0] = ry[0] = 0.0;
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
                s = 0;
#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
                double fx = f(sn, tn00);/*double fx = f(sn, tn10);*/
                //double fx = 0.5 * ( f(sn, tn00)+f(sn, tn10) );

                ay[s]  = 0.0;
                by[s]  = /*beta  **/ b211;
                cy[s]  = /*beta  **/ b213;
                dy[s]  = /*beta  **/ b214 * p05[s][i];
                dy[s] += /*beta  **/ b215 * ((p05[s][i+1]-2.0*p05[s][i]+p05[s][i-1])/(hx*hx));
                dy[s] += /*beta  **/ b216 * ((p05[s][i+1]-p05[s][i-1])/(2.0*hx));
                dy[s] += /*beta  **/ ht05 * fx;
                dy[s] += /*gamma **/ b217 * value;
#else
                // O(h1)
                ay[s]  = 0.0;
                by[s]  = +/*beta*/1.0 + /*alpha **/ hx;
                cy[s]  = -/*beta*/1.0;
                dy[s]  = +/*gamma **/ hx * value;
#endif
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                const double alpha = condition.alpha();
                const double beta  = condition.beta();
                //const double gamma = condition.gamma();
                const double gamma = 1.0;

                s = 0;
#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
                double fx = f(sn, tn00);/*double fx = f(sn, tn10);*/
                //double fx = 0.5 * ( f(sn, tn00)+f(sn, tn10) );

                ay[s]  = 0.0;
                by[s]  = beta  * b211 + alpha * b212;
                cy[s]  = beta  * b213;
                dy[s]  = beta  * b214 * p05[s][i];
                dy[s] += beta  * b215 * ((p05[s][i+1]-2.0*p05[s][i]+p05[s][i-1])/(hx*hx));
                dy[s] += beta  * b216 * ((p05[s][i+1]-p05[s][i-1])/(2.0*hx));
                dy[s] += beta  * ht05 * fx;
                dy[s] += gamma * b217 * value;
#else
                ay[s]  = 0.0;
                by[s]  = +beta + alpha * hx;
                cy[s]  = -beta;
                dy[s]  = +gamma * hx * value;
#endif
            }

            sn.j = ymax; sn.y = ymax*hy;
            value = boundary(sn, tn10, condition);

            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                e = M-1;
                p10[M][i] = /*(gamma/alpha)**/value;
                dy[M-1] -= k23 * p10[M][i];
                cy[M-1] = ay[M] = by[M] = cy[M] = dy[M] = ry[M] = 0.0;
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
                e = M;
#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
                double fx = f(sn, tn00);/*double fx = f(sn, tn10);*/
                //double fx = 0.5 * ( f(sn, tn00)+f(sn, tn10) );

                ay[e]  = /*beta  **/ b221;
                by[e]  = /*beta  **/ b223;
                cy[e]  = 0.0;
                dy[e]  = /*beta  **/ b224 * p05[e][i];
                dy[e] += /*beta  **/ b225 * ((p05[e][i+1]-2.0*p05[e][i]+p05[e][i-1])/(hx*hx));
                dy[e] += /*beta  **/ b226 * ((p05[e][i+1]-p05[e][i-1])/(2.0*hx));
                dy[e] += /*beta  **/ ht05 * fx;
                dy[e] += /*gamma **/ b227 * value;
#else
                // O(h1)
                ay[e]  = -/*beta*/1.0;
                by[e]  = +/*beta*/1.0 + /*alpha **/ hx;
                cy[e]  = 0.0;
                dy[e]  = +/*gamma **/ hx * value;
#endif
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                const double alpha = condition.alpha();
                const double beta  = condition.beta();
                //const double gamma = condition.gamma();
                const double gamma = 1.0;

                e = M;
#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
                double fx = f(sn, tn00);/*double fx = f(sn, tn10);*/
                //double fx = 0.5 * ( f(sn, tn00)+f(sn, tn10) );

                ay[e]  = beta  * b221;
                by[e]  = beta  * b223 + alpha * b222;
                cy[e]  = 0.0;
                dy[e]  = beta  * b224 * p05[e][i];
                dy[e] += beta  * b225 * ((p05[e][i+1]-2.0*p05[e][i]+p05[e][i-1])/(hx*hx));
                dy[e] += beta  * b226 * ((p05[e][i+1]-p05[e][i-1])/(2.0*hx));
                dy[e] += beta  * ht05 * fx;
                dy[e] += gamma * b227 * value;
#else
                // O(h1)
                ay[e]  = -beta;
                by[e]  = +beta + alpha * hx;
                cy[e]  = 0.0;
                dy[e]  = +gamma * hx * value;
#endif
            }

            tomasAlgorithmLeft2Right(ay+s, by+s, cy+s, dy+s, ry+s, e-s+1);
            for (unsigned int j=s; j<=e; j++) p10[j][i] = ry[j];
            //memcpy(p10[i]+s, ry+s, sizeof(double*)*(e-s+1));
        }

        for (m=ymin, sn.j=m, sn.y=m*hy, j=0; m<=ymax; ++m, sn.j=m, sn.y=m*hy, ++j)
        {
            sn.i = xmin; sn.x = xmin*hx;
            value = boundary(sn, tn10, condition);

            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                p10[j][0] = /*(gamma/alpha)**/value;
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
                p10[j][0] = (3.5*p10[j][1] - 2.0*p10[j][2] + 0.5*p10[j][3] + hx/**(gamma/beta)*/*value)/(2.0);
#else
                p10[j][0] = p10[j][1] + hx/**(gamma/beta)*/*value;
#endif
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                const double alpha = condition.alpha();
                const double beta  = condition.beta();
                //const double gamma = condition.gamma();
                const double gamma = 1.0;

#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
                p10[j][0] = (3.5*p10[j][1] - 2.0*p10[j][2] + 0.5*p10[j][3] + hx*(gamma/beta)*value)/(2.0 + (alpha/beta)*hx);
#else
                p10[j][0] = (p10[j][1] + hx*(gamma/beta)*value)/(1.0 + hx*(alpha/beta));
#endif
            }

            sn.i = xmax; sn.x = xmax*hx;
            value = boundary(sn, tn10, condition);

            if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
            {
                p10[j][N] = /*(gamma/alpha)**/value;
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
            {
#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
                p10[j][N] = (3.5*p10[j][N-1] - 2.0*p10[j][N-2] + 0.5*p10[j][N-3] + hx/**(gamma/beta)*/*value)/(2.0);
#else
                p10[j][N] = p10[j][N-1] + hx/**(gamma/beta)*/*value;
#endif
            }
            else if (condition.boundaryCondition() == BoundaryCondition::Robin)
            {
                const double alpha = condition.alpha();
                const double beta  = condition.beta();
                //const double gamma = condition.gamma();
                const double gamma = 1.0;

#ifdef PARABOLIC_IBVP_H_D2V1_BORDER_O2
                p10[j][N] = (3.5*p10[j][N-1] - 2.0*p10[j][N-2] + 0.5*p10[j][N-3] + hx*(gamma/beta)*value)/(2.0 + (alpha/beta)*hx);
#else
                p10[j][N] = (p10[j][N-1] + hx*(gamma/beta)*value)/(1.0 + hx*(alpha/beta));
#endif
            }
        }

        layerInfo(DoubleMatrix(p10, M+1, N+1), tn10);

        /**************************************************** y direction apprx ***************************************************/
        double **_tmp = p00; p00 = p10; p10 = _tmp;
    }

    for (unsigned int i=0; i<=M; i++) { free(p00[i]); free(p05[i]); free(p10[i]); }
    free(p00); free(p05); free(p10);

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

void ILoadedHeatEquationFBVP::setLoadedPoints(const std::vector<LoadedSpacePoint> &loadedPoints)
{
    this->_loadedPoints = loadedPoints;
}

const std::vector<LoadedSpacePoint> ILoadedHeatEquationFBVP::loadedPoints() const
{
    return _loadedPoints;
}

//------------------------------------------------------------------------------------------------------------------------------------------------------//
