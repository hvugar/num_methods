#include "solver2.h"

using namespace p3p0;

void Solver2::Main(int /*argc*/, char **/*argv*/)
{
    Solver2 s;
    s.forward.setThermalDiffusivity(1.0);
    s.forward.setThermalConvection(0.0);
    s.forward.setThermalConductivity(0.0);

    s.backward.setThermalDiffusivity(-1.0);
    s.backward.setThermalConvection(0.0);
    s.backward.setThermalConductivity(0.0);

    const size_t time_size = s.timeDimension().size();

    /******************************************************************
     *                     Calculating V(x) function
     ******************************************************************/
    {
        DoubleVector x(s.timeDimension().size()*2);
        for (unsigned int i=0; i<time_size; i++)
        {
            x[0*time_size + i] = (i*s.timeDimension().step());
            x[1*time_size + i] = (i*s.timeDimension().step())*(i*s.timeDimension().step());
        }
        IPrinter::printVector(x.mid(0*time_size, 1*time_size-1));
        IPrinter::printVector(x.mid(1*time_size, 2*time_size-1));
        s.convert1(x, 2, s.timeDimension().size());
        s.forward.implicit_calculate_D1V1();
        IPrinter::print(s.U);
        s.V = s.U;
    }

    /******************************************************************
     *                     Calculating V(x) function
     ******************************************************************/
    DoubleVector x(s.timeDimension().size()*2, 2.0);
    x[0*time_size] = 0.0;
    x[1*time_size] = 0.0;
    DoubleVector g(x.length());
    DoubleVector g1(x.length());

    puts("---");
    IPrinter::printVector(x.mid(0*time_size, 1*time_size-1));
    IPrinter::printVector(x.mid(1*time_size, 2*time_size-1));

    puts("---");
    puts("OK");
    s.gradient(x, g);
    puts("OK");
    IPrinter::printVector(g.mid(0*time_size, 1*time_size-1).L2Normalize());
    IPrinter::printVector(g.mid(1*time_size, 2*time_size-1).L2Normalize());

    puts("---");
    IGradient::Gradient(&s, 0.01, x, g1);
    IPrinter::printVector(g1.mid(0*time_size, 1*time_size-1).L2Normalize());
    IPrinter::printVector(g1.mid(1*time_size, 2*time_size-1).L2Normalize());

    //    ConjugateGradient gm;
    SteepestDescentGradient gm;
    gm.setFunction(&s);
    gm.setGradient(&s);
    gm.setPrinter(&s);
    //gm.setProjection(&s);
    //gm.setGradientNormalizer(&prob);
    gm.setOptimalityTolerance(0.0);
    gm.setStepTolerance(0.0);
    gm.setFunctionTolerance(0.0);
    gm.setR1MinimizeEpsilon(0.01, 0.0001);
//    gm.setMaxIterationCount(10);
    gm.setNormalize(true);
    gm.showExitMessage(true);
    gm.calculate(x);

    puts("---");
    IPrinter::print(x.mid(0*time_size, 1*time_size-1), time_size);
    IPrinter::print(x.mid(1*time_size, 2*time_size-1), time_size);
    puts("+++");
}

HeatEquationIBVP::HeatEquationIBVP() {}

HeatEquationIBVP::HeatEquationIBVP(const HeatEquationIBVP &) {}

HeatEquationIBVP & HeatEquationIBVP::operator =(const HeatEquationIBVP &) { return *this; }

HeatEquationIBVP::~HeatEquationIBVP() {}

auto HeatEquationIBVP::q(const TimeNodePDE &tn) const -> DoubleVector
{
    DoubleVector q;
    q << s->q[0][tn.i];
    q << s->q[1][tn.i];
    return q;
}

auto HeatEquationIBVP::v(size_t /*i*/, const PointNodeODE &/*tn*/, SpacePoint &/*vl*/) const -> void { }

auto HeatEquationIBVP::z(const TimeNodePDE &tn) const -> DoubleVector
{
    DoubleVector z;
    z << 0.05 + 0.9*sin(M_PI*tn.t);
    z << 0.95 - 0.9*sin(M_PI*tn.t);
    return z;
}

auto HeatEquationIBVP::initial(const SpaceNodePDE &/*sn*/, InitialCondition /*cn*/) const -> double { return 0.0; }

auto HeatEquationIBVP::boundary(const SpaceNodePDE &sn, const TimeNodePDE &/*tn*/, BoundaryConditionPDE &bc) const -> double
{
    if (sn.i == spaceDimensionX().min()) { bc = BoundaryConditionPDE::Robin(s->lambda1(), -1.0, s->lambda1()); return s->theta(); }
    if (sn.i == spaceDimensionX().max()) { bc = BoundaryConditionPDE::Robin(s->lambda1(), +1.0, s->lambda1()); return s->theta(); }
    throw std::exception();
}

auto HeatEquationIBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double
{
    const auto sigma = 0.01;
    DoubleVector _z = z(tn);
    DoubleVector _q = q(tn);
    auto fx = 0.0;
    fx += _q[0] * DeltaFunction::gaussian(sn.x, _z[0], sigma);
    fx += _q[1] * DeltaFunction::gaussian(sn.x, _z[1], sigma);
    return fx;
}

auto HeatEquationIBVP::layerInfo(const DoubleVector &u, const TimeNodePDE &tn) const -> void
{
    //printf("%4d | >> ", tn.i); IPrinter::printVector(u);
    //DoubleVector _z = z(tn);
    //DoubleVector _q = q(tn);
    if (static_cast<int>(tn.i) == timeDimension().max()) { s->U = u; }
}

auto HeatEquationIBVP::timeDimension() const -> Dimension { return s->timeDimension(); }

auto HeatEquationIBVP::spaceDimensionX() const -> Dimension { return s->spaceDimensionX(); }

auto HeatEquationIBVP::spaceDimensionY() const -> Dimension { return spaceDimensionX(); }

auto HeatEquationIBVP::spaceDimensionZ() const -> Dimension { return spaceDimensionX(); }

/**************************************************************************************************************************/

HeatEquationFBVP::HeatEquationFBVP() {}

HeatEquationFBVP::HeatEquationFBVP(const HeatEquationFBVP &) {}

HeatEquationFBVP & HeatEquationFBVP::operator =(const HeatEquationFBVP &) { return *this; }

HeatEquationFBVP::~HeatEquationFBVP() {}

auto HeatEquationFBVP::final(const SpaceNodePDE &sn, FinalCondition /*ic*/) const -> double
{
    size_t n = static_cast<size_t>(sn.i);
    return -2.0 * s->mu(sn) * ( s->U[n] - s->V[n] );
}

auto HeatEquationFBVP::boundary(const SpaceNodePDE &sn, const TimeNodePDE &/*tn*/, BoundaryConditionPDE &bc) const -> double
{
    if (sn.i == spaceDimensionX().min()) { bc = BoundaryConditionPDE::Robin(s->lambda1(), -1.0, s->lambda1()); return s->theta(); }
    if (sn.i == spaceDimensionX().max()) { bc = BoundaryConditionPDE::Robin(s->lambda1(), +1.0, s->lambda1()); return s->theta(); }
    throw std::exception();
}

auto HeatEquationFBVP::f(const SpaceNodePDE &/*sn*/, const TimeNodePDE &/*tn*/) const -> double
{
    return 0.0;
}

auto HeatEquationFBVP::HeatEquationFBVP::layerInfo(const DoubleVector &p, const TimeNodePDE &tn) const -> void
{
    printf("%4d | << ", tn.i); IPrinter::printVector(p);

    const auto hx = spaceDimensionX().step();
    const auto Nx = spaceDimensionX().size();
    DoubleVector _z = h->z(tn);
    s->p[0][tn.i] = DeltaFunction::lumpedPoint4(p, _z[0], hx, Nx);
    s->p[1][tn.i] = DeltaFunction::lumpedPoint4(p, _z[1], hx, Nx);
}

auto HeatEquationFBVP::timeDimension() const -> Dimension { return s->timeDimension(); }

auto HeatEquationFBVP::spaceDimensionX() const -> Dimension { return s->spaceDimensionX(); }

auto HeatEquationFBVP::spaceDimensionY() const -> Dimension { return spaceDimensionX(); }

auto HeatEquationFBVP::spaceDimensionZ() const -> Dimension { return spaceDimensionX(); }

/**************************************************************************************************************************/

Common::~Common() {}

auto Common::convert1(const DoubleVector &x, size_t size, size_t length) -> void
{
    if (q != nullptr)
    {
        for (size_t i=0; i<size; i++)
        {
            q[i].clear();
        }
        q = nullptr;
    }

    q = new DoubleVector[size];

    for (size_t i=0; i<size; i++)
    {
        q[i].resize(length);
        for (size_t n=0; n<length; n++)
        {
            q[i][n] = x[i*length+n];
        }
    }
}

auto Common::convert2(size_t size, size_t length, DoubleVector &x) const -> void
{
    x.clear();
    x.resize(size*length);

    if (q != nullptr)
    {
        for (size_t i=0; i<size; i++)
        {
            const DoubleVector &qi = q[i];
            for (size_t n=0; n<length; n++)
            {
                x[i*length+n] = qi[n];
            }
        }
    }
}

/**************************************************************************************************************************/

Solver2::Solver2()
{
    const auto time_min = 0;
    const auto time_max = 1000;
    const auto time_len = 1001;
    const auto time_stp = 0.01;

    const auto dimX_min = 0;
    const auto dimX_max = 100;
    const auto dimX_stp = 0.01;
    const auto dimX_len = 101;

    _timeDimension = Dimension(time_stp, time_min, time_max);
    _spaceDimensionX = Dimension(dimX_stp, dimX_min, dimX_max);

    U.resize(dimX_len, 0.0);
    V.resize(dimX_len, 0.0);

    p = new DoubleVector[heatSourceNumber];
    for (size_t i=0; i<heatSourceNumber; i++) { p[i].resize(time_len); }

    forward.s = this;
    backward.s = this;
    backward.h = &forward;
}

auto Solver2::gradient(const DoubleVector &x, DoubleVector &g) const -> void
{
    const auto time_size = timeDimension().size();
    const_cast<Solver2*>(this)->convert1(x, 2, time_size);

    printf("%d %d\n", time_size, x.length());
    forward.implicit_calculate_D1V1();
    puts("OK2");
    backward.implicit_calculate_D1V1();
    puts("OK3");

    g.resize(x.length(), 0.0);

    for (size_t i=0; i<heatSourceNumber; i++)
    {
        const DoubleVector &pi = p[i];
        for (size_t ln=0; ln<time_size; ln++) { g[i*time_size+ln] = -pi[ln]; }
        g[i*time_size] = 0.0;
    }
}

auto Solver2::fx(const DoubleVector &x) const -> double
{
    const auto time_size = timeDimension().size();
    const_cast<Solver2*>(this)->convert1(x, 2, time_size);
    forward.implicit_calculate_D1V1();

    const auto N = spaceDimensionX().size()-1;
    const auto h = spaceDimensionX().step();

    double sum = (U[0]-V[0])*(U[0]-V[0]);
    for (unsigned int i=1; i<N; i++)
    {
        sum += 0.5*(U[i]-V[i])*(U[i]-V[i]);
    }
    sum += (U[N]-V[N])*(U[N]-V[N]);

    return sum*h;
}

auto Solver2::print(unsigned int iteration, const DoubleVector &x, const DoubleVector &/*g*/, double f, double /*alpha*/, GradientBasedMethod::MethodResult /*result*/) const -> void
{
    const size_t time_size = timeDimension().size();

    IPrinter::printSeperatorLine();
    printf_s("I[%4d] %10.8f %10.8f\n", iteration, f, fx(x));
    IPrinter::printVector(x.mid(0*time_size, 1*time_size-1), "q1");
    IPrinter::printVector(x.mid(1*time_size, 2*time_size-1), "q2");
    IPrinter::printVector(U, "U:");
    IPrinter::printVector(V, "V:");
    IPrinter::printSeperatorLine();
}

/**************************************************************************************************************************/



