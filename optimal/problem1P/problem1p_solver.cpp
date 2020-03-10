#include "problem1p_solver.h"

void tomasAlgorithmRight2LeftBackward(const double *a, const double *b, const double *c, const double *d, double *x, unsigned int N, const double *e);
void tomasAlgorithmRight2LeftForward(const double *a, const double *b, const double *c, const double *d, double *x, unsigned int N, double *e, double f);


using namespace p1p;

void ProblemSolver::Main(int argc, char* argv[])
{
    C_UNUSED(argc);
    C_UNUSED(argv);

    unsigned int length = 101;
    ProblemSolver solver(Dimension(0.01, 0, 100), Dimension(0.01, 0, 100));

    const size_t pointCount = 4;
    solver.mPointCount = pointCount;
    solver.mPoints = new SpacePoint[pointCount];
    solver.mPoints[0].x = 0.2023;
    solver.mPoints[1].x = 0.5543;
    solver.mPoints[2].x = 0.7983;
    solver.mPoints[3].x = 0.91235;

    solver.k = new double [pointCount];
    solver.z = new double [pointCount];

    DoubleVector x(3*pointCount, 0.0);
    for (unsigned int i=0; i<pointCount; i++)
    {
        x[i+pointCount*0] = sin(0.5+i);
        x[i+pointCount*1] = (2.0+i)+sin(i+1.0);
        x[i+pointCount*2] = solver.mPoints[i].x;
    }

    solver.u_vl = new double* [pointCount];
    solver.u_dx = new double* [pointCount];
    for (unsigned int i=0; i<pointCount; i++)
    {
        solver.u_vl[i] = new double[length];
        solver.u_dx[i] = new double[length];
    }

    IPrinter::printSeperatorLine("x");
    IPrinter::printVector(x, "xx: ");

    DoubleVector ng(x.length(), 0.0);
    DoubleVector ag(x.length(), 0.0);

    solver.gradient(x, ag);
    IGradient::Gradient(&solver, 0.01, x, ng);

    //ag[0]   = ng[0]   = 0.0;
    //ag[100] = ng[100] = 0.0;

    //        ag[4] = ag[5] = ag[6] = ag[7] = ag[8] = ag[9] = ag[10] = ag[11] = 0.0;
    //        ng[4] = ng[5] = ng[6] = ng[7] = ng[8] = ng[9] = ng[10] = ng[11] = 0.0;
    //    ag[0] = ag[1] = ag[2] = ag[3] = ag[8] = ag[9] = ag[10] = ag[11] = 0.0;
    //    ng[0] = ng[1] = ng[2] = ng[3] = ng[8] = ng[9] = ng[10] = ng[11] = 0.0;
    //    ag[0] = ag[1] = ag[2] = ag[3] = ag[4] = ag[5] = ag[6] = ag[7] = 0.0;
    //    ng[0] = ng[1] = ng[2] = ng[3] = ng[4] = ng[5] = ng[6] = ng[7] = 0.0;
    //    ag[8] = ag[9] = ag[10] = ag[11] = 0.0;
    //    ng[8] = ng[9] = ng[10] = ng[11] = 0.0;

    IPrinter::printSeperatorLine("gradients", '=');
    IPrinter::printVector(ag, "ag: ");
    IPrinter::printVector(ng, "ng: ");
    IPrinter::printSeperatorLine("normolized", '=');
    DoubleVector ka = ag.mid(0,3); ka.EuclideanNormalize();
    DoubleVector kn = ng.mid(0,3); kn.EuclideanNormalize();
    IPrinter::printVector(ka, "ag: ", ka.length());
    IPrinter::printVector(kn, "ng: ", kn.length());
    IPrinter::printSeperatorLine(nullptr, '=');
    DoubleVector za = ag.mid(4,7); za.EuclideanNormalize();
    DoubleVector zn = ng.mid(4,7); zn.EuclideanNormalize();
    IPrinter::printVector(za, "ag: ", za.length());
    IPrinter::printVector(zn, "ng: ", zn.length());
    IPrinter::printSeperatorLine(nullptr, '=');
    DoubleVector xa = ag.mid(8,11); xa.EuclideanNormalize();
    DoubleVector xn = ng.mid(8,11); xn.EuclideanNormalize();
    IPrinter::printVector(xa, "ag: ", xa.length());
    IPrinter::printVector(xn, "ng: ", xn.length());
    IPrinter::printSeperatorLine(nullptr, '=');

    //ConjugateGradient g;
    //g.setNormalize(false);
    SteepestDescentGradient g;
    g.setNormalize(true);
    g.setFunction(&solver);
    g.setGradient(&solver);
    g.setPrinter(&solver);
    g.setProjection(&solver);
    g.setOptimalityTolerance(-0.0);
    g.setFunctionTolerance(-0.0);
    g.setStepTolerance(-0.0);
    g.setR1MinimizeEpsilon(0.1, 0.0001);
    g.setMaxIterationCount(100);
    g.showExitMessage(true);
    g.calculate(x);
}

ProblemSolver::ProblemSolver(const Dimension &timeDimension, const Dimension &spaceDimensionX)
{
    forward.solver = backward.solver = this;
    forward._userHalfValues = backward._userHalfValues = true;

    setTimeDimension(timeDimension);
    setSpaceDimensionX(spaceDimensionX);

    initialTemperature = 0.0;
    environmentTemperature = 0.5;

    forward.setThermalDiffusivity(thermalDiffusivity);
    forward.setThermalConvection(-thermalConvection);
    forward.setThermalConductivity(0.0);

    backward.setThermalDiffusivity(-thermalDiffusivity);
    backward.setThermalConvection(thermalConvection);
    backward.setThermalConductivity(0.0);

    const_this = const_cast<ProblemSolver*>(this);
}

void ProblemSolver::setTimeDimension(const Dimension &timeDimension)
{
    _timeDimension = timeDimension;
    p0.resize(timeDimension.size());
}

void ProblemSolver::setSpaceDimensionX(const Dimension &spaceDimensionX)
{
    _spaceDimensionX = spaceDimensionX;

    unsigned int size = spaceDimensionX.size();
    U.resize(size, 0.0);
    V.resize(size, 5.0);
}

void ProblemSolver::gradient(const DoubleVector &x, DoubleVector &g) const
{
    unsigned int M = timeDimension().size()-1;
    double ht = timeDimension().step();

    const size_t pc = mPointCount;
    for (unsigned int i=0; i<pc; i++)
    {
        k[i] = x[i+mPointCount*0];
        z[i] = x[i+mPointCount*1];
        mPoints[i].x = x[i+mPointCount*2];
    }
    g.resize(x.length());

    forward.implicit_calculate_D1V1();
    backward.implicit_calculate_D1V1();

    for (unsigned int i=0; i<mPointCount; i++)
    {
        double integral1 = 0.0;
        double integral2 = 0.0;
        double integral3 = 0.0;

        integral1 += 0.5*ht*p0[0]*(u_vl[i][0]-z[i]);
        integral2 += 0.5*ht*p0[0];
        integral3 += 0.5*ht*p0[0]*u_dx[i][0];
        for (unsigned int j=1; j<=M-1; j++)
        {
            integral1 += ht*p0[j]*(u_vl[i][j]-z[i]);
            integral2 += ht*p0[j];
            integral3 += ht*p0[j]*u_dx[i][j];
        }

        integral1 += 0.5*ht*p0[M]*(u_vl[i][M]-z[i]);
        integral2 += 0.5*ht*p0[M];
        integral3 += 0.5*ht*p0[M]*u_dx[i][M];

        g[i+mPointCount*0] = +thermalConductivity1*thermalDiffusivity*integral1;
        g[i+mPointCount*1] = -thermalConductivity1*thermalDiffusivity*integral2*k[i];
        g[i+mPointCount*2] = 0.0;//+thermalConductivity1*thermalDiffusivity*integral3*k[i];
    }
}

double ProblemSolver::fx(const DoubleVector &x) const
{
    double sum = integral(x);
    return sum;
}

double ProblemSolver::integral(const DoubleVector &x) const
{
    const double hx = _spaceDimensionX.step();
    const unsigned int N = _spaceDimensionX.size() - 1;

    const uint32_t pointCount = mPointCount;
    for (unsigned int i=0; i<pointCount; i++)
    {
        k[i] = x[i+pointCount*0];
        z[i] = x[i+pointCount*1];
        mPoints[i].x = x[i+pointCount*2];
    }

    forward.implicit_calculate_D1V1();

    double integ_sum = 0.0;
    integ_sum += 0.5*(U[0]-V[0])*(U[0]-V[0]);
    for (unsigned int n=1; n<=N-1; n++)
    {
        integ_sum += (U[n]-V[n])*(U[n]-V[n]);
    }
    integ_sum += 0.5*(U[N]-V[N])*(U[N]-V[N]);
    integ_sum *= hx;

    return integ_sum;
}

void ProblemSolver::project(DoubleVector &x, unsigned int)
{
    unsigned int s = mPointCount*2;
    unsigned int f = mPointCount*3;
//    for (unsigned int i=0; i<mPointCount; i++)
//    {
//        if (x[i] >= +5.00) x[i] = +5.00;
//        if (x[i] <= -5.00) x[i] = -5.00;
//    }

//    for (unsigned int i=mPointCount; i<2*mPointCount; i++)
//    {
//        if (x[i] >= +5.00) x[i] = +5.00;
//        if (x[i] <= -5.00) x[i] = -5.00;
//    }

    for (unsigned int i=s; i<f; i++)
    {
        if (x[i] <= 0.05) x[i] = 0.05;
        if (x[i] >= 0.95) x[i] = 0.95;
    }
}

void ProblemSolver::print(unsigned int iteration, const DoubleVector &x, const DoubleVector &g, double f, double alpha, GradientMethod::MethodResult result) const
{
    printf("%6d %10.6f ", iteration, f);
    IPrinter::printVector(x, "xx: ");
}

/*********************************************************************************************************************************/

double ProblemSolver::frw_initial(const SpaceNodePDE &, InitialCondition) const { return initialTemperature; }

double ProblemSolver::frw_boundary(const SpaceNodePDE &sn, const TimeNodePDE &, BoundaryConditionPDE &c) const
{
    if (sn.i == spaceDimensionX().min())
    {
        c = BoundaryConditionPDE(BoundaryCondition::Robin, thermalConductivity1, -1.0, thermalConductivity1);
        return 0.0;
    }
    if (sn.i == spaceDimensionX().max())
    {
        c = BoundaryConditionPDE(BoundaryCondition::Robin, thermalConductivity2, +1.0, thermalConductivity2);
        return environmentTemperature;
    }

    throw std::runtime_error("frw_boundary error");
}

double ProblemSolver::frw_f(const SpaceNodePDE &, const TimeNodePDE &) const { return thermalConvection * environmentTemperature; }

void ProblemSolver::frw_layerInfo(const DoubleVector &u, const TimeNodePDE &tn) const
{
    if (static_cast<int>(tn.i) == timeDimension().max())
    {
        const_cast<ProblemSolver*>(this)->U = u;
    }

    auto const_solver = const_cast<ProblemSolver *>(this);
    auto hx = spaceDimensionX().step();
    auto N = spaceDimensionX().size()-1;
    auto xmin = spaceDimensionX().min();
    auto xmax = spaceDimensionX().max();
    const double sigma = hx;

    double w1 = WEIGHT;
    double w2 = 1.0-WEIGHT;

    if (static_cast<int>(tn.i) == timeDimension().min())
    {
        unsigned int j=0;
        for (int n=xmin; n<=xmax; n++, j++)
        {
            double x = n*hx;
            for (unsigned int i=0; i<mPointCount; i++)
            {
                unsigned int idx = static_cast<unsigned>(round(mPoints[i].x*N));
                const_solver->u_vl[i][tn.i] = w2*exp(-((x-mPoints[i].x)*(x-mPoints[i].x))/(2.0*sigma*sigma))/(sigma*sqrt(2.0*M_PI))*hx*u[j];
                const_solver->u_dx[i][tn.i] = w2*((u[idx+1]-u[idx-1])/(2.0*hx));
            }
        }
    }
    else
    {
        unsigned int j=0;
        for (int n=xmin; n<=xmax; n++, j++)
        {
            double x = n*hx;
            for (unsigned int i=0; i<mPointCount; i++)
            {
                const_solver->u_vl[i][tn.i-1] += w2*exp(-((x-mPoints[i].x)*(x-mPoints[i].x))/(2.0*sigma*sigma))/(sigma*sqrt(2.0*M_PI))*hx*u[j];
                const_solver->u_vl[i][tn.i-0]  = w1*exp(-((x-mPoints[i].x)*(x-mPoints[i].x))/(2.0*sigma*sigma))/(sigma*sqrt(2.0*M_PI))*hx*u[j];

                unsigned int idx = static_cast<unsigned>(round(mPoints[i].x*N));
                const_solver->u_dx[i][tn.i-1] += w2*((u[idx+1]-u[idx-1])/(2.0*hx));
                const_solver->u_dx[i][tn.i]    = w1*((u[idx+1]-u[idx-1])/(2.0*hx));
            }
        }
    }
}

/*********************************************************************************************************************************/

double ProblemSolver::bcw_final(const SpaceNodePDE &sn, FinalCondition) const
{
    unsigned int i = static_cast<unsigned int>(sn.i);
    return -2.0*(U[i] - V[i]);
}

double ProblemSolver::bcw_boundary(const SpaceNodePDE &sn, const TimeNodePDE &, BoundaryConditionPDE &c) const
{
    if (sn.i == spaceDimensionX().min())
    {
        c = BoundaryConditionPDE(BoundaryCondition::Robin, thermalConductivity1, -1.0, 0.0);
        return 0.0;
    }

    if (sn.i == spaceDimensionX().max())
    {
        c = BoundaryConditionPDE(BoundaryCondition::Robin, thermalConductivity2, +1.0, 0.0);
        return 0.0;
    }

    throw std::runtime_error("bcw_boundary error:");
}

double ProblemSolver::bcw_f(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }

void ProblemSolver::bcw_layerInfo(const DoubleVector &p, const TimeNodePDE &tn) const
{
    double w1 = WEIGHT;
    double w2 = 1.0-WEIGHT;

    auto const_solver = const_cast<ProblemSolver *>(this);

    if (static_cast<int>(tn.i) == timeDimension().max())
    {
        const_solver->p0[tn.i] = w2*p[0];
    }
    else
    {
        const_solver->p0[tn.i+1] += w1*p[0];
        const_solver->p0[tn.i+0]  = w2*p[0];
    }
    //const_solver->p0[tn.i] = p[0];
}

/*********************************************************************************************************************************/

double HeatEquationIBVP::initial(const SpaceNodePDE &sn, InitialCondition c) const { return solver->frw_initial(sn, c); }

double HeatEquationIBVP::boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &c) const { return solver->frw_boundary(sn, tn, c); }

double HeatEquationIBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const { return solver->frw_f(sn, tn); }

void HeatEquationIBVP::layerInfo(const DoubleVector &u, const TimeNodePDE &tn) const { return solver->frw_layerInfo(u, tn); }

Dimension HeatEquationIBVP::timeDimension() const { return solver->timeDimension(); }

Dimension HeatEquationIBVP::spaceDimensionX() const { return solver->spaceDimensionX(); }

Dimension HeatEquationIBVP::spaceDimensionY() const { return solver->spaceDimensionX(); }

Dimension HeatEquationIBVP::spaceDimensionZ() const { return solver->spaceDimensionX(); }

void HeatEquationIBVP::implicit_calculate_D1V1() const
{
    const double *k = solver->k;
    const double *z = solver->z;
    const SpacePoint *mPoints = solver->mPoints;
    const unsigned int mPointsCount = solver->mPointCount;

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
    const double b17 = -((2.0*a*ht)/hx)*w1 + b*ht*w1;
    const double b18 = -((2.0*a*ht)/hx)*w2 + b*ht*w2;

    // right border condition parameters
    const double b21 = -((2.0*a*ht)/(hx*hx))*w1;
    const double b22 = +((2.0*a*ht)/hx)*w1 + b*ht*w1;
    const double b23 = +1.0 + ((2.0*a*ht)/(hx*hx))*w1 - c*ht*w1;
    const double b24 = +((2.0*a*ht)/(hx*hx))*w2;
    const double b25 = -((2.0*a*ht)/hx)*w2 - b*ht*w2;
    const double b26 = +1.0 - ((2.0*a*ht)/(hx*hx))*w2 + c*ht*w2;
    const double b29 = +((2.0*a*ht)/hx) + b*ht;
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
    double *ex = static_cast<double*>(malloc(sizeof(double)*(N+1)));

    for (unsigned int n=0; n<=N; n++)
    {
        ax[n] = k11;
        bx[n] = k12;
        cx[n] = k13;
    }
    ax[0] = 0.0; cx[N] = 0.0;

    DoubleVector u0(N+1);
    DoubleVector u1(N+1);

    /***********************************************************************************************/
    /***********************************************************************************************/

    SpaceNodePDE sn;
    unsigned int i = 0;
    for (int n=xmin; n<=xmax; n++, i++)
    {
        sn.i = n; sn.x = n*hx;
        u0[i] = initial(sn, InitialCondition::InitialValue);
    }
    layerInfo(u0, TimeNodePDE(0, 0.0));

    /***********************************************************************************************/
    /***********************************************************************************************/

    for (unsigned int ln=1; ln<=L; ln++)
    {
        TimeNodePDE tn0; tn0.i = ln-1; tn0.t = tn0.i*ht;
        TimeNodePDE tn1; tn1.i = ln;   tn1.t = tn1.i*ht;

        unsigned int i = 1;
        for (int n=xmin+1; n<=xmax-1; n++, i++)
        {
            sn.i = n; sn.x = n*hx;
            dx[i] = k21*u0[i-1] + k22*u0[i] + k23*u0[i+1];

            if (_userHalfValues)
            {
                dx[i] += ht*f(sn, tn1);
            }
            else
            {
                dx[i] += ht_w1*f(sn, tn1)+ht_w2*f(sn, tn0);
                throw std::runtime_error("error!");
            }
        }

        unsigned int s=0, e=N;
        BoundaryConditionPDE condition; double value, alpha, beta, gamma;
        BoundaryConditionPDE condition0; double value0;

        sn.i = xmin; sn.x = xmin*hx;
        value = boundary(sn, tn1, condition);
        //value0 = boundary(sn, tn0, condition);
        alpha = condition.alpha();
        beta  = condition.beta();
        gamma = condition.gamma();
        if (condition.boundaryCondition() == BoundaryCondition::Dirichlet) { assert(condition.boundaryCondition() != BoundaryCondition::Dirichlet); }
        else if (condition.boundaryCondition() == BoundaryCondition::Neumann) { assert(condition.boundaryCondition() != BoundaryCondition::Neumann); }
        else if (condition.boundaryCondition() == BoundaryCondition::Robin)
        {
            s = 0;
            ax[0]  = 0.0;
            bx[0]  = beta*b11 + alpha*b12;
            cx[0]  = beta*b13;
            dx[0]  = beta*b14*u0[0] + alpha*b15*u0[0] + beta*b16*u0[1];

            value = 0.0; for (unsigned int m=0; m<mPointsCount; m++) value -= k[m]*z[m];

            for (unsigned int j=0; j<=N; j++)
            {
                ex[j] = 0.0;
                for (unsigned int m=0; m<mPointsCount; m++)
                {
                    if (j == static_cast<unsigned int>(mPoints[m].x/hx))
                    {
                        ex[j] -= (gamma*b19)*k[m];
                    }
                }
            }
            ex[0] += bx[0]; bx[0] = 0.0;
            ex[1] += cx[0]; cx[0] = 0.0;

            if (_userHalfValues)
            {
                dx[0] += gamma*b19*value;
                dx[0] += beta*ht*f(sn, tn1);
            }
            else
            {
                dx[0] += gamma*(b17*value+b18*value0);
                dx[0] += beta*(ht_w1*f(sn, tn1)+ht_w2*f(sn, tn0));
                throw std::runtime_error("error!");
            }
        }

        sn.i = xmax; sn.x = xmax*hx;
        value = boundary(sn, tn1, condition);
        value0 = boundary(sn, tn0, condition);
        alpha = condition.alpha();
        beta  = condition.beta();
        gamma = condition.gamma();
        if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
        {
            assert(condition.boundaryCondition() != BoundaryCondition::Dirichlet);

            e = N-1;
            u1[N] = (gamma/alpha)*value;
            dx[N-1] -= k13 * u1[N];
            cx[N-1] = ax[N] = bx[N] = cx[N] = dx[N] = rx[N] = 0.0;
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
        {
            assert(condition.boundaryCondition() != BoundaryCondition::Neumann);

            e = N;
            ax[e]  = beta*b21;
            bx[e]  = beta*b23;
            cx[e]  = 0.0;
            dx[e]  = beta*b24*u0[e-1] + beta*b26*u0[e];

            if (_userHalfValues)
            {
                dx[e] += beta*ht*f(sn, tn1);
                dx[e] += gamma*b29*value;
            }
            else
            {
                dx[e] += beta*(ht_w1*f(sn, tn1)+ht_w2*f(sn, tn0));
                dx[e] += gamma*(b27*value+b28*value0);
            }
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Robin)
        {
            e = N;
            ax[e]  = beta*b21;
            bx[e]  = alpha*b22 + beta*b23;
            cx[e]  = 0.0;
            dx[e]  = beta*b24*u0[e-1] + alpha*b25*u0[e] + beta*b26*u0[e];

            if (_userHalfValues)
            {
                dx[e]  += beta*ht*f(sn, tn1);
                dx[e] += gamma*b29*value;
            }
            else
            {
                dx[e]  += beta*(ht_w1*f(sn, tn1)+ht_w2*f(sn, tn0));
                dx[e] += gamma*(b27*value+b28*value0);
            }
        }

        tomasAlgorithmRight2LeftForward(ax+s, bx+s, cx+s, dx+s, rx+s, e-s+1, ex+s, dx[0]);
        for (unsigned int n=s; n<=e; n++) u1[n] = rx[n];

        layerInfo(u1, tn1);

        u0 = u1;
    }

    u0.clear();
    u1.clear();

    free(rx);
    free(ex);
    free(dx);
    free(cx);
    free(bx);
    free(ax);
}

/*********************************************************************************************************************************/

double HeatEquationFBVP::final(const SpaceNodePDE &sn, FinalCondition c) const { return solver->bcw_final(sn, c); }

double HeatEquationFBVP::boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &c) const { return solver->bcw_boundary(sn, tn, c); }

double HeatEquationFBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const { return solver->bcw_f(sn, tn); }

void HeatEquationFBVP::layerInfo(const DoubleVector &p, const TimeNodePDE &tn) const { return solver->bcw_layerInfo(p, tn); }

Dimension HeatEquationFBVP::timeDimension() const { return solver->timeDimension(); }

Dimension HeatEquationFBVP::spaceDimensionX() const { return solver->spaceDimensionX(); }

Dimension HeatEquationFBVP::spaceDimensionY() const { return solver->spaceDimensionX(); }

Dimension HeatEquationFBVP::spaceDimensionZ() const { return solver->spaceDimensionX(); }

void HeatEquationFBVP::implicit_calculate_D1V1() const
{
    const double *k = solver->k;
    //const double *z = solver->z;
    const SpacePoint *mPoints = solver->mPoints;
    const unsigned int mPointCount = solver->mPointCount;

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
    const double b19 = -((2.0*a*ht)/hx) + b*ht;
    const double b17 = -((2.0*a*ht)/hx)*w1 + b*ht*w1;
    const double b18 = -((2.0*a*ht)/hx)*w2 + b*ht*w2;

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
    const double ht_w1 = +ht*w1;
    const double ht_w2 = +ht*w2;

    double *ax = static_cast<double*>(malloc(sizeof(double)*(N+1)));
    double *bx = static_cast<double*>(malloc(sizeof(double)*(N+1)));
    double *cx = static_cast<double*>(malloc(sizeof(double)*(N+1)));
    double *dx = static_cast<double*>(malloc(sizeof(double)*(N+1)));
    double *rx = static_cast<double*>(malloc(sizeof(double)*(N+1)));
    double *ex = static_cast<double*>(malloc(sizeof(double)*(N+1)));

    for (unsigned int n=0; n<=N; n++)
    {
        ax[n] = k11;
        bx[n] = k12;
        cx[n] = k13;
    }
    ax[0] = 0.0; cx[N] = 0.0;

    DoubleVector p0(N+1);
    DoubleVector p1(N+1);

    /***********************************************************************************************/
    /***********************************************************************************************/

    SpaceNodePDE sn;
    unsigned int i = 0;
    for (int n=xmin; n<=xmax; n++, i++)
    {
        sn.i = n; sn.x = sn.i*hx;
        p0[i] = final(sn, FinalCondition::FinalValue);
    }
    layerInfo(p0, TimeNodePDE(M, M*ht));

    /***********************************************************************************************/
    /***********************************************************************************************/

    for (unsigned int ln=M-1, size_ln = static_cast<unsigned int>(0)-1; ln != size_ln; ln--)
    {
        TimeNodePDE tn0; tn0.i = ln+1; tn0.t = tn0.i*ht;
        TimeNodePDE tn1; tn1.i = ln;   tn1.t = tn1.i*ht;

        unsigned int i = 1;
        for (int n=xmin+1; n<=xmax-1; n++, i++)
        {
            sn.i = n; sn.x = n*hx;
            dx[i] = k21*p0[i-1] + k22*p0[i] + k23*p0[i+1];

            if (_userHalfValues)
            {
                dx[i] += ht*f(sn, tn1);
            }
            else
            {
                dx[i] += ht_w1*f(sn, tn1)+ht_w2*f(sn, tn0);
                throw std::runtime_error("_userHalfValues");
            }
        }

        unsigned int s=0, e=N;
        BoundaryConditionPDE condition; double alpha, beta, gamma, value;
        BoundaryConditionPDE condition0; double value0;

        sn.i = xmin; sn.x = xmin*hx;
        value = boundary(sn, tn1, condition);
        value0 = boundary(sn, tn0, condition);
        alpha = condition.alpha();
        beta  = condition.beta();
        gamma = condition.gamma();
        //        printf(">>> %f %f %f\n", alpha, beta, gamma);

        if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
        {
            s = 1;
            p1[0] = (gamma/alpha)*value;
            dx[1] -= k11 * p1[0];
            ax[1] = ax[0] = bx[0] = cx[0] = dx[0] = rx[0] = 0.0;

            throw std::runtime_error("BoundaryCondition::Dirichlet");
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
        {
            s = 0;
            ax[s]  = 0.0;
            bx[s]  = beta*b11;
            cx[s]  = beta*b13;
            dx[s]  = beta*b14*p0[s] + beta*b16*p0[s+1];

            if (_userHalfValues)
            {
                dx[s] += beta*w1*f(sn, tn1);
                dx[s] += gamma*b19*value;
            }
            else
            {
                dx[s] += beta*(ht_w1*f(sn, tn1)+ht_w2*f(sn, tn0));
                dx[s] += gamma*(b17*value+b18*value0);
                throw std::runtime_error("_userHalfValues");
            }

            throw std::runtime_error("BoundaryCondition::Neumann");
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Robin)
        {
            s = 0;
            ax[0]  = 0.0;
            bx[0]  = beta*b11 + alpha*b12;
            cx[0]  = beta*b13;
            dx[0]  = beta*b14*p0[0] + alpha*b15*p0[0] + beta*b16*p0[1];

            if (_userHalfValues)
            {
                dx[0] += gamma*b19*value;
                dx[0] += beta*ht*f(sn, tn1);
            }
            else
            {
                dx[0] +=  gamma*(b17*value+b18*value0);
                dx[0] += beta*(ht_w1*f(sn, tn1)+ht_w2*f(sn, tn0));
                throw std::runtime_error("error!");
            }
        }


        ///////////////////////////////////////////////////////////////////////////////////////
        i = 0;
        for (int n=xmin; n<=xmax; n++, i++)
        {
            ex[i] = 0.0;
            for (unsigned int j=0; j<mPointCount; j++)
            {
                ex[i] -= a*alpha*k[j]*(exp( -((n*hx-mPoints[j].x)*(n*hx-mPoints[j].x))/(2.0*hx*hx) )/(hx*sqrt(2.0*M_PI)));
            }
        }
        ///////////////////////////////////////////////////////////////////////////////////////

        sn.i = xmax; sn.x = xmax*hx;
        value = boundary(sn, tn1, condition);
        value0 = boundary(sn, tn0, condition);
        alpha = condition.alpha();
        beta  = condition.beta();
        gamma = condition.gamma();
        if (condition.boundaryCondition() == BoundaryCondition::Dirichlet)
        {
            e = N-1;
            p1[N] = (gamma/alpha)*value;
            dx[N-1] -= k13 * p1[N];
            cx[N-1] = ax[N] = bx[N] = cx[N] = dx[N] = rx[N] = 0.0;

            throw std::runtime_error("BoundaryCondition::Dirichlet");
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Neumann)
        {
            e = N;
            ax[e]  = beta*b21;
            bx[e]  = beta*b23;
            cx[e]  = 0.0;
            dx[e]  = beta*b24*p0[e-1] + beta*b26*p0[e];

            if (_userHalfValues)
            {
                dx[e] += beta*w1*f(sn, tn1);
                dx[e] += gamma*b29*value;
            }
            else
            {
                dx[e] += beta*(ht_w1*f(sn, tn1)+ht_w2*f(sn, tn0));
                dx[e] += gamma*(b27*value+b28*value0);
            }
            throw std::runtime_error("BoundaryCondition::Neumann");
        }
        else if (condition.boundaryCondition() == BoundaryCondition::Robin)
        {
            e = N;
            ax[N]  = beta*b21;
            bx[N]  = alpha*b22 + beta*b23;
            cx[N]  = 0.0;
            dx[N]  = beta*b24*p0[N-1] + alpha*b25*p0[N] + beta*b26*p0[N];

            if (_userHalfValues)
            {
                dx[N] += beta*ht*f(sn, tn1);
                dx[N] += gamma*b29*value;
            }
            else
            {
                dx[N] += beta*(ht_w1*f(sn, tn1)+ht_w2*f(sn, tn0));
                dx[N] += gamma*(b27*value+b28*value0);
                throw std::runtime_error("_userHalfValues");
            }
        }

        tomasAlgorithmRight2LeftBackward(ax+s, bx+s, cx+s, dx+s, rx+s, e-s+1, ex);
        for (unsigned int n=s; n<=e; n++) p1[n] = rx[n];

        layerInfo(p1, tn1);

        p0 = p1;
    }

    p0.clear();
    p1.clear();

    free(rx);
    free(ex);
    free(dx);
    free(cx);
    free(bx);
    free(ax);
}

/*********************************************************************************************************************************/

void tomasAlgorithmRight2LeftBackward(const double *a, const double *b, const double *c, const double *d, double *x, unsigned int N, const double *e)
{
    const unsigned int M = N-1;
    double *alpha = static_cast<double*>(malloc(sizeof(double)*N));
    double *betta = static_cast<double*>(malloc(sizeof(double)*N));
    double *gamma = static_cast<double*>(malloc(sizeof(double)*N));
    unsigned int i = 1;
    double m = 0.0;

    /* Прямой ход метода прогонки. Определение прогоночных коэффициентов. */
    if (fabs(b[0]) <= DBL_EPSILON) { fprintf(stderr, "ERROR: %6d b0=%20.8f\n", 0, fabs(b[0])); assert(!(fabs(b[0]) <= 0.0)); }
    if (fabs(b[M]) <= DBL_EPSILON) { fprintf(stderr, "ERROR: %6d bN=%20.8f\n", 0, fabs(b[M])); assert(!(fabs(b[M]) <= 0.0)); }

    if (fabs(c[0]) < DBL_EPSILON) { fprintf(stderr, "ERROR: %6d c0=%20.8f<0.0\n", 0, fabs(c[0])); assert(!(fabs(c[0]) < 0.0)); }
    if (fabs(a[M]) < DBL_EPSILON) { fprintf(stderr, "ERROR: %6d aN=%20.8f<0.0\n", 0, fabs(a[M])); assert(!(fabs(a[M]) < 0.0)); }

    if (fabs(b[0]) < fabs(c[0])) { fprintf(stderr, "ERROR: %6d %20.8f %20.8f\n", 0, fabs(b[0]), fabs(c[0])); assert(!(fabs(b[0]) < fabs(c[0]))); }
    if (fabs(b[M]) < fabs(a[M])) { fprintf(stderr, "ERROR: %6d %20.8f %20.8f\n", 0, fabs(b[M]), fabs(a[M])); assert(!(fabs(b[M]) < fabs(a[M]))); }
    if (fabs(b[M]) < fabs(e[M])) { fprintf(stderr, "ERROR: %6d %20.8f %20.8f\n", 0, fabs(b[M]), fabs(e[M])); assert(!(fabs(b[M]) < fabs(e[M]))); }

    alpha[M-1] = -a[M]/b[M];
    betta[M-1] = +d[M]/b[M];
    gamma[M-1] = -e[M]/b[M];

    //if (fabs(gamma[M-1]) > 1.0) { fprintf(stderr, "ERROR: %6d %20.8f %20.8f %20.8f %20.8f %20.8f\n", M-1,
    //                                      fabs(gamma[M-1]), e[M], c[i], betta[i], 1.0); assert(!(fabs(gamma[M-1]) > 1.0)); }

    //printf(">> %4d %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n", M, a[M], b[M], c[M], d[M], e[M], alpha[M-1], betta[M-1], gamma[M-1]);
    for (i=M-1; i>=1; i--)
    {
        if (fabs(b[i]) < fabs(a[i])+fabs(c[i])) { fprintf(stderr, "ERROR: %6d %20.8f %20.8f %20.8f\n", i, fabs(a[i]), fabs(b[i]), fabs(c[i])); assert(fabs(b[i]) < fabs(a[i])+fabs(c[i])); }
        if (fabs(a[i]) <= DBL_EPSILON) { fprintf(stderr, "ERROR: %6d %20.8f\n", i, fabs(a[i])); assert(fabs(a[i]) <= 0.0); }
        if (fabs(c[i]) <= DBL_EPSILON) { fprintf(stderr, "ERROR: %6d %20.8f\n", i, fabs(c[i])); assert(fabs(c[i]) <= 0.0); }

        m = b[i] + c[i]*alpha[i];
        alpha[i-1] = -a[i]/m;
        betta[i-1] = +(d[i]-c[i]*betta[i])/m;
        gamma[i-1] = -(e[i]+c[i]*betta[i])/m;

        //printf(">> %4d %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n", i, a[i], b[i], c[i], d[i], e[i], alpha[i-1], betta[i-1], gamma[i-1]);
        //if (fabs(gamma[i-1]) > 1.0) { fprintf(stderr, "ERROR: %6d %20.8f %20.8f %20.8f %20.8f %20.8f\n", i+1,
        //                                      fabs(gamma[i]), e[i+1], c[i+1], betta[i], 0.0); assert(!(fabs(gamma[i-1]) > 1.0)); }
    }

    /**********************************************************************/

    /* Обратный ход метода прогонки. Обратный ход метода прогонки начинается
     * с вычисления х0. */
    /* Остальные значения неизвестных находятся рекуррентно. */
    m = b[0] + c[0]*(alpha[0]+gamma[0]);
    x[0] = (d[0]-c[0]*betta[0])/m;
    for (i=0; i<M; i++)
    {
        x[i+1] = alpha[i]*x[i] + betta[i] + gamma[0]*x[0];
    }

    /**********************************************************************/

    free(gamma); gamma=nullptr;
    free(betta); betta=nullptr;
    free(alpha); alpha=nullptr;
}

void tomasAlgorithmRight2LeftForward(const double *a, const double *b, const double *c, const double *d, double *x, unsigned int N, double *e, double f)
{
    const unsigned int M = N-1;
    double *alpha = static_cast<double*>( malloc(sizeof(double)*N) );
    double *betta = static_cast<double*>( malloc(sizeof(double)*N) );
    unsigned int i = 1;
    double m = 0.0;

    /* Прямой ход метода прогонки. Определение прогоночных коэффициентов. */
    //    if (fabs(b[0]) <= DBL_EPSILON) { fprintf(stderr, "ERROR: %6d b0=%20.8f\n", 0, fabs(b[0])); assert(fabs(b[0]) <= 0.0); }
    if (fabs(b[M]) <= DBL_EPSILON) { fprintf(stderr, "ERROR: %6d bN=%20.8f\n", 0, fabs(b[M])); assert(!(fabs(b[M]) <= 0.0)); }

    //    if (fabs(c[0]) < DBL_EPSILON) { fprintf(stderr, "ERROR: %6d c0=%20.8f<0.0\n", 0, fabs(c[0])); assert(fabs(c[0]) < 0.0); }
    if (fabs(a[M]) < DBL_EPSILON) { fprintf(stderr, "ERROR: %6d aN=%20.8f<0.0\n", 0, fabs(a[M])); assert(!(fabs(a[M]) < 0.0)); }

    //    if (fabs(b[0]) < fabs(c[0])) { fprintf(stderr, "ERROR: %6d %20.8f %20.8f\n", 0, fabs(b[0]), fabs(c[0])); assert(fabs(b[0]) < fabs(c[0])); }
    if (fabs(b[M]) < fabs(a[M])) { fprintf(stderr, "ERROR: %6d %20.8f %20.8f\n", 0, fabs(b[M]), fabs(a[M])); assert(!(fabs(b[M]) < fabs(a[M]))); }

    alpha[M-1] = -a[M]/b[M];
    betta[M-1] = +d[M]/b[M];
    for (i=M-1; i>=1; i--)
    {
        if (fabs(b[i]) < fabs(a[i])+fabs(c[i])) { fprintf(stderr, "ERROR: %6d %20.8f %20.8f %20.8f\n", i, fabs(a[i]), fabs(b[i]), fabs(c[i])); assert(!(fabs(b[i]) < fabs(a[i])+fabs(c[i]))); }
        if (fabs(a[i]) <= DBL_EPSILON) { fprintf(stderr, "ERROR: %6d %20.8f\n", i, fabs(a[i])); assert(!(fabs(a[i]) <= 0.0)); }
        if (fabs(c[i]) <= DBL_EPSILON) { fprintf(stderr, "ERROR: %6d %20.8f\n", i, fabs(c[i])); assert(!(fabs(c[i]) <= 0.0)); }

        m = b[i] + c[i]*alpha[i];
        alpha[i-1] = -a[i]/m;
        betta[i-1] = +(d[i]-c[i]*betta[i])/m;

        e[i] += alpha[i]*e[i+1];
        f    -= betta[i]*e[i+1];
    }

    /**********************************************************************/

    /* Обратный ход метода прогонки. Обратный ход метода прогонки начинается
     * с вычисления х0. */
    /* Остальные значения неизвестных находятся рекуррентно. */
    m = e[0] + e[1]*alpha[0];
    x[0] = (f-e[1]*betta[0])/m;
    for (i=0; i<M; i++)
    {
        x[i+1] = alpha[i]*x[i] + betta[i];
    }

    /**********************************************************************/

    free(betta); betta=nullptr;
    free(alpha); alpha=nullptr;
}
