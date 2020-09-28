#include "solver2.h"

using namespace p3p0;

void Solver2::Main(int /*argc*/, char **/*argv*/)
{
    Solver2 s;
    double a = 0.01;
    s.forward.setThermalDiffusivity(a);
    s.forward.setThermalConvection(0.0);
    s.forward.setThermalConductivity(0.0);

    s.backward.setThermalDiffusivity(-a);
    s.backward.setThermalConvection(0.0);
    s.backward.setThermalConductivity(0.0);

    const auto time_size = s.timeDimension().size();
    const auto time_step = s.timeDimension().step();
    const auto vctr_size = time_size;
    const unsigned int start[] = { 0*vctr_size+0, 1*vctr_size+0 };
    const unsigned int finsh[] = { 1*vctr_size-1, 2*vctr_size-1 };

    /******************************************************************
     *                     Calculating V(x) function
     ******************************************************************/
    {
        DoubleVector x(vctr_size*s.heatSourceNumber);
        for (unsigned int ln=0; ln<vctr_size; ln++)
        {
            double t = ln * time_step;
            DoubleVector _qNorm = s.qNorm1(t);
            x[0*vctr_size + ln] = _qNorm[0];
            x[1*vctr_size + ln] = _qNorm[1];
        }
        IPrinter::printVector(s._w, s._p, x.mid(start[0], finsh[0]), "q1:");
        IPrinter::printVector(s._w, s._p, x.mid(start[1], finsh[1]), "q2:");
        s.convert1(x, s.heatSourceNumber, vctr_size);
        s.forward.implicit_calculate_D1V1();
        IPrinter::printVector(s._w, s._p, s.U, "V: ");
        s.V = s.U;
    }

    /******************************************************************
     *                     Calculating V(x) function
     ******************************************************************/
    DoubleVector x(vctr_size*s.heatSourceNumber, 2.0);
    //x[0*vctr_size] = 0.0;
    //x[1*vctr_size] = 0.0;

    //    for (unsigned int ln=0; ln<vctr_size; ln++)
    //    {
    //        double t = ln * time_step;
    //        DoubleVector _qNorm = s.qNorm1(t);
    //        x[0*vctr_size + ln] = _qNorm[0] * 1.1;//( rand()%2==0 ? 0.8 : 1.2 );
    //        x[1*vctr_size + ln] = _qNorm[1] * 1.1;//( rand()%2==0 ? 0.8 : 1.2 );
    //    }

    DoubleVector g(x.length());
    DoubleVector g1(x.length());

    puts("---");
    IPrinter::printVector(s._w, s._p, x.mid(start[0], finsh[0]), "q1:");
    IPrinter::printVector(s._w, s._p, x.mid(start[1], finsh[1]), "q2:");

    puts("---");
    s.gradient(x, g);
    IPrinter::printVector(s._w, s._p, g.mid(start[0], finsh[0]).L2Normalize(), "g1:");
    IPrinter::printVector(s._w, s._p, g.mid(start[1], finsh[1]).L2Normalize(), "g2:");

    puts("---");
    IGradient::Gradient(&s, 0.01, x, g1);
    IPrinter::printVector(s._w, s._p, g1.mid(start[0], finsh[0]).L2Normalize(), "g1:");
    IPrinter::printVector(s._w, s._p, g1.mid(start[1], finsh[1]).L2Normalize(), "g2:");

    //    puts("Starting optimization...");
    //    //ConjugateGradient gm;
    //    SteepestDescentGradient gm;
    //    gm.setFunction(&s);
    //    gm.setGradient(&s);
    //    gm.setPrinter(&s);
    //    //gm.setProjection(&s);
    //    //gm.setGradientNormalizer(&prob);
    //    gm.setOptimalityTolerance(0.000001);
    //    gm.setStepTolerance(0.000001);
    //    gm.setFunctionTolerance(0.000001);
    //    gm.setR1MinimizeEpsilon(0.01, 0.001);
    //    //gm.setMaxIterationCount(10);
    //    gm.setNormalize(false);
    //    gm.showExitMessage(true);
    //    gm.calculate(x);

    //    puts("Optimization is finished.");
    //    IPrinter::printVector(s._w, s._p, x.mid(0*vctr_size, 1*vctr_size-1));
    //    IPrinter::printVector(s._w, s._p, x.mid(1*vctr_size, 2*vctr_size-1));
    //    puts("+++");
}

HeatEquationIBVP::HeatEquationIBVP() {}

HeatEquationIBVP::HeatEquationIBVP(const HeatEquationIBVP &h) : IHeatEquationIBVP (h) { *this = *this; }

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
    z << 0.05 + 0.9*fabs(sin(M_PI*tn.t));
    z << 0.95 - 0.9*fabs(sin(M_PI*tn.t));
    //z << 0.25;// + 0.9*fabs(sin(M_PI*tn.t));
    //z << 0.75;// - 0.9*fabs(sin(M_PI*tn.t));
    return z;
}

auto HeatEquationIBVP::initial(const SpaceNodePDE &/*sn*/, InitialCondition /*cn*/) const -> double { return 0.0; }

auto HeatEquationIBVP::boundary(const SpaceNodePDE &sn, const TimeNodePDE &/*tn*/, BoundaryConditionPDE &bc) const -> double
{
    if (sn.i == spaceDimensionX().min()) { bc = BoundaryConditionPDE::Robin(s->lambda1(), -1.0, s->lambda1()); return s->theta(); }
    if (sn.i == spaceDimensionX().max()) { bc = BoundaryConditionPDE::Robin(s->lambda1(), +1.0, s->lambda1()); return s->theta(); }
    //if (sn.i == spaceDimensionX().min()) { bc = BoundaryConditionPDE::Neumann(+1.0, 0.0); return 0.0; }
    //if (sn.i == spaceDimensionX().max()) { bc = BoundaryConditionPDE::Neumann(+1.0, 0.0); return 0.0; }
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
    //if (static_cast<int>(tn.i) == 1) { IPrinter::printVector(u); }

    if (static_cast<int>(tn.i) == timeDimension().max()) { s->U = u; }
}

auto HeatEquationIBVP::timeDimension() const -> Dimension { return s->timeDimension(); }

auto HeatEquationIBVP::spaceDimensionX() const -> Dimension { return s->spaceDimensionX(); }

auto HeatEquationIBVP::spaceDimensionY() const -> Dimension { return spaceDimensionX(); }

auto HeatEquationIBVP::spaceDimensionZ() const -> Dimension { return spaceDimensionX(); }

/**************************************************************************************************************************/

HeatEquationFBVP::HeatEquationFBVP() {}

HeatEquationFBVP::HeatEquationFBVP(const HeatEquationFBVP &h) : IHeatEquationFBVP (h) {}

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
    //if (sn.i == spaceDimensionX().min()) { bc = BoundaryConditionPDE::Neumann(+1.0, 0.0); return 0.0; }
    //if (sn.i == spaceDimensionX().max()) { bc = BoundaryConditionPDE::Neumann(+1.0, 0.0); return 0.0; }
    throw std::exception();
}

auto HeatEquationFBVP::f(const SpaceNodePDE &/*sn*/, const TimeNodePDE &/*tn*/) const -> double
{
    return 0.0;
}

auto HeatEquationFBVP::HeatEquationFBVP::layerInfo(const DoubleVector &psi, const TimeNodePDE &tn) const -> void
{
    const auto hx = spaceDimensionX().step();
    const auto Nx = spaceDimensionX().size();
    DoubleVector _z = h->z(tn);
    s->p[0][tn.i] = DeltaFunction::lumpedPoint4(psi, _z[0], hx, Nx);
    s->p[1][tn.i] = DeltaFunction::lumpedPoint4(psi, _z[1], hx, Nx);
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
    const auto time_max = 200;
    const auto time_len = time_max+1;
    const auto time_stp = 0.005;

    const auto dimX_min = 0;
    const auto dimX_max = 100;
    const auto dimX_len = dimX_max+1;
    const auto dimX_stp = 0.01;

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
    const auto time_step = timeDimension().step();
    const_cast<Solver2*>(this)->convert1(x, heatSourceNumber, time_size);

    forward.implicit_calculate_D1V1();
    backward.implicit_calculate_D1V1();

    g.resize(x.length(), 0.0);

    for (size_t ln=0; ln<time_size; ln++)
    {
        double t = ln * time_step;
        DoubleVector _qNorm = qNorm1(t);
        for (size_t i=0; i<heatSourceNumber; i++)
        {
            g[i*time_size+ln] = -p[i][ln] + 2.0*epsilon*(q[i][ln] - _qNorm[i]*no_norm);
        }
    }

    //for (size_t i=0; i<heatSourceNumber; i++) { g[i*time_size] = 0.0; }
}

auto Solver2::fx(const DoubleVector &x) const -> double
{
    return integral(x) + norm(x)*epsilon;
}

auto Solver2::integral(const DoubleVector &x) const -> double
{
    const auto time_size = timeDimension().size();
    const_cast<Solver2*>(this)->convert1(x, heatSourceNumber, time_size);
    forward.implicit_calculate_D1V1();
    const auto N = spaceDimensionX().size()-1;

    double sum = 0.0;

    sum += (U[0]-V[0])*(U[0]-V[0]);
    for (unsigned int i=1; i<N; i++)
    {
        sum += 0.5*(U[i]-V[i])*(U[i]-V[i]);
    }
    sum += (U[N]-V[N])*(U[N]-V[N]);

    return sum*spaceDimensionX().step();
}

auto Solver2::norm(const DoubleVector &x) const -> double
{
    const auto time_size = timeDimension().size();
    const auto time_step = timeDimension().step();
    const_cast<Solver2*>(this)->convert1(x, heatSourceNumber, time_size);

    double _norm = 0.0;
    size_t ln = 0;
    double t = ln * time_step;
    DoubleVector _qNorm = qNorm1(t);
    for (size_t i=0; i<heatSourceNumber; i++) { _norm += (q[i][ln]-_qNorm[i]*no_norm)*(q[i][ln]-_qNorm[i]*no_norm); }
    for (ln=1; ln<time_size-1; ln++)
    {
        t = ln * time_step;
        _qNorm = qNorm1(t);
        for (unsigned int i=0; i<heatSourceNumber; i++) { _norm += 0.5*(q[i][ln]-_qNorm[i]*no_norm)*(q[i][ln]-_qNorm[i]*no_norm); }
    }
    ln = time_size-1;
    t = ln * time_step;
    _qNorm = qNorm1(t);
    for (size_t i=0; i<heatSourceNumber; i++) { _norm += (q[i][ln]-_qNorm[i]*no_norm)*(q[i][ln]-_qNorm[i]*no_norm); }

    return _norm*time_step;
}

auto Solver2::qNorm1(double t) const -> DoubleVector
{
    DoubleVector q;
    q << (t+2.0)*(t+2.0);
    q << (t+2.0)*(t+2.0)*(t+2.0);
    return q;
}

auto Solver2::print(unsigned int it, const DoubleVector &x, const DoubleVector &g, double /*f*/, double _alpha, GradientBasedMethod::MethodResult result) const -> void
{
    const size_t time_size = timeDimension().size();
    const auto vector_size = time_size;
    IPrinter::printSeperatorLine();

    const auto _fx = fx(x);
    const auto _integral = integral(x);
    const auto _norm = norm(x);
    printf_s("I[%4d] fx: %10.8f int: %10.8f nrm: %10.8f %10.8f %d\n", it, _fx, _integral, _norm, _alpha, result);
    IPrinter::printVector(_w, _p, x.mid(0*vector_size, 1*vector_size-1), "q1");
    IPrinter::printVector(_w, _p, x.mid(1*vector_size, 2*vector_size-1), "q2");
    IPrinter::printVector(_w, _p, g.mid(0*vector_size, 1*vector_size-1), "g1");
    IPrinter::printVector(_w, _p, g.mid(1*vector_size, 2*vector_size-1), "g2");
    IPrinter::printVector(_w, _p, U, "U:");
    IPrinter::printVector(_w, _p, V, "V:");
}

/**************************************************************************************************************************/
/**************************************************************************************************************************/

using namespace p3p2;

CommonParameters::~CommonParameters() {}

auto CommonParameters::q(const TimeNodePDE &tn) const -> DoubleVector { return mq[tn.i]; }

auto CommonParameters::z(const TimeNodePDE &tn) const -> DoubleVector
{
    DoubleVector z;
    //z << 0.05 + 0.9*fabs(sin(M_PI*tn.t));
    //z << 0.95 - 0.9*fabs(sin(M_PI*tn.t));
    z << 0.30;// + 0.9*fabs(sin(M_PI*tn.t));
    z << 0.70;// - 0.9*fabs(sin(M_PI*tn.t));
    return z;
    //return mz[tn.i];
}

auto CommonParameters::convert1(const DoubleVector &x, size_t heatSourceNumber, size_t time_size) -> void
{
    for (size_t ln=0; ln<time_size; ln++)
    {
        mq[ln].resize(heatSourceNumber, 0.0);
        for (size_t i=0; i<heatSourceNumber; i++)
        {
            mq[ln][i] = x[i*time_size+ln];
        }
    }
}

auto CommonParameters::convert2(size_t size, size_t length, DoubleVector &x) const -> void
{
    x.clear();
    x.resize(size*length);

    if (mq != nullptr)
    {
        for (size_t ln=0; ln<length; ln++)
        {
            for (size_t i=0; i<heatSourceNumber; i++)
            {
                x[i*length+ln] = mq[ln][i];
            }
        }
    }
}

auto CommonParameters::qNorm1(double t) const -> DoubleVector
{
    DoubleVector q;
    q << (t+2.0)*(t+2.0);
    q << (t+2.0)*(t+2.0)*(t+2.0);
    return q;
}

/**************************************************************************************************************************/

p3p2::HeatEquationIBVP1::~HeatEquationIBVP1() {}

auto p3p2::HeatEquationIBVP1::calculate_forward(const DoubleVector &/*x*/) const -> void
{
    IHeatEquationIBVP::implicit_calculate_D1V1();
}

auto p3p2::HeatEquationIBVP1::initial(const SpaceNodePDE &/*sn*/, InitialCondition /*ic*/) const -> double { return m_b; }

auto p3p2::HeatEquationIBVP1::boundary(const SpaceNodePDE &sn, const TimeNodePDE &/*tn*/, BoundaryConditionPDE &bc) const -> double
{
    if (sn.i == spaceDimensionX().min()) { bc = BoundaryConditionPDE::Robin(lambda1(), -1.0, lambda1()); return theta(); }
    if (sn.i == spaceDimensionX().max()) { bc = BoundaryConditionPDE::Robin(lambda1(), +1.0, lambda1()); return theta(); }
    //if (sn.i == spaceDimensionX().min()) { bc = BoundaryConditionPDE::Neumann(+1.0, 0.0); return 0.0; }
    //if (sn.i == spaceDimensionX().max()) { bc = BoundaryConditionPDE::Neumann(+1.0, 0.0); return 0.0; }
    throw std::exception();
}

auto p3p2::HeatEquationIBVP1::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double
{
    const auto sigma = 0.01;
    DoubleVector _z = z(tn);
    DoubleVector _q = q(tn);
    auto fx = 0.0;
    fx += _q[0] * DeltaFunction::gaussian(sn.x, _z[0], sigma);
    fx += _q[1] * DeltaFunction::gaussian(sn.x, _z[1], sigma);
    return fx;
}

auto p3p2::HeatEquationIBVP1::layerInfo(const DoubleVector &u, const TimeNodePDE &tn) const -> void
{
    IPrinter::printVector(u, "HeatEquationIBVP1");
    if (static_cast<int>(tn.i) == timeDimension().max()) { const_cast<p3p2::HeatEquationIBVP1*>(this)->U = u; }
}

auto p3p2::HeatEquationIBVP1::timeDimension() const -> Dimension { return CommonParameters::timeDimension(); }

auto p3p2::HeatEquationIBVP1::spaceDimensionX() const -> Dimension { return CommonParameters::spaceDimensionX(); }

auto p3p2::HeatEquationIBVP1::spaceDimensionY() const -> Dimension { return spaceDimensionX(); }

auto p3p2::HeatEquationIBVP1::spaceDimensionZ() const -> Dimension { return spaceDimensionX(); }

/**************************************************************************************************************************/

p3p2::HeatEquationFBVP1::~HeatEquationFBVP1() {}

auto p3p2::HeatEquationFBVP1::calculate_backward(const DoubleVector &/*x*/) const -> void
{
    IHeatEquationFBVP::implicit_calculate_D1V1();
}

auto p3p2::HeatEquationFBVP1::final(const SpaceNodePDE &sn, FinalCondition /*fc*/) const -> double
{
    size_t n = static_cast<size_t>(sn.i);
    return -2.0 * mu(sn) * ( U[n] - V[n] );}

auto p3p2::HeatEquationFBVP1::boundary(const SpaceNodePDE &sn, const TimeNodePDE &/*tn*/, BoundaryConditionPDE &bc) const -> double
{
    if (sn.i == spaceDimensionX().min()) { bc = BoundaryConditionPDE::Robin(lambda1(), -1.0, lambda1()); return theta(); }
    if (sn.i == spaceDimensionX().max()) { bc = BoundaryConditionPDE::Robin(lambda1(), +1.0, lambda1()); return theta(); }
    //if (sn.i == spaceDimensionX().min()) { bc = BoundaryConditionPDE::Neumann(+1.0, 0.0); return 0.0; }
    //if (sn.i == spaceDimensionX().max()) { bc = BoundaryConditionPDE::Neumann(+1.0, 0.0); return 0.0; }
    throw std::exception();

}

auto p3p2::HeatEquationFBVP1::f(const SpaceNodePDE &/*sn*/, const TimeNodePDE &/*tn*/) const -> double { return 0.0; }

auto p3p2::HeatEquationFBVP1::layerInfo(const DoubleVector &psi, const TimeNodePDE &tn) const -> void
{
    IPrinter::printVector(psi, "HeatEquationFBVP1");
//    std::cout << "auto p3p2::HeatEquationFBVP1::layerInfo(const DoubleVector &psi, const TimeNodePDE &tn) const -> void" << std::endl;
//    const auto hx = spaceDimensionX().step();
//    const auto Nx = spaceDimensionX().size();
//    DoubleVector _z = z(tn);
//    mp[0][tn.i] = DeltaFunction::lumpedPoint4(psi, _z[0], hx, Nx);
//    mp[1][tn.i] = DeltaFunction::lumpedPoint4(psi, _z[1], hx, Nx);
}

auto p3p2::HeatEquationFBVP1::timeDimension() const -> Dimension { return CommonParameters::timeDimension(); }

auto p3p2::HeatEquationFBVP1::spaceDimensionX() const -> Dimension { return CommonParameters::spaceDimensionX(); }

auto p3p2::HeatEquationFBVP1::spaceDimensionY() const -> Dimension { return spaceDimensionX(); }

auto p3p2::HeatEquationFBVP1::spaceDimensionZ() const -> Dimension { return spaceDimensionX(); }

/**************************************************************************************************************************/

p3p2::Functional::Functional(double thermalDiffusivity, double thermalConvection, double thermalConductivity)
{
    IHeatEquationIBVP::setThermalDiffusivity(thermalDiffusivity);
    IHeatEquationIBVP::setThermalConvection(thermalConvection);
    IHeatEquationIBVP::setThermalConductivity(thermalConductivity);

    IHeatEquationFBVP::setThermalDiffusivity(-thermalDiffusivity);
    IHeatEquationFBVP::setThermalConvection(thermalConvection);
    IHeatEquationFBVP::setThermalConductivity(thermalConductivity);
}

p3p2::Functional::~Functional() {}

auto p3p2::Functional::fx(const DoubleVector &x) const -> double
{
    return integral(x) + norm(x)*epsilon;
}

auto p3p2::Functional::integral(const DoubleVector &x) const -> double
{
    const auto time_size = CommonParameters::timeDimension().size();
    const_cast<p3p2::Functional*>(this)->convert1(x, heatSourceNumber, time_size);
    //HeatEquationIBVP1::implicit_calculate_D1V1();
    const auto N = CommonParameters::spaceDimensionX().size()-1;

    double sum = 0.0;

    sum += (U[0]-V[0])*(U[0]-V[0]);
    for (unsigned int i=1; i<N; i++)
    {
        sum += 0.5*(U[i]-V[i])*(U[i]-V[i]);
    }
    sum += (U[N]-V[N])*(U[N]-V[N]);

    return sum*CommonParameters::spaceDimensionX().step();
}

auto p3p2::Functional::norm(const DoubleVector &x) const -> double
{
    const auto time_size = CommonParameters::timeDimension().size();
    const auto time_step = CommonParameters::timeDimension().step();
    const_cast<p3p2::Functional*>(this)->convert1(x, heatSourceNumber, time_size);

    double _norm = 0.0;
    size_t ln = 0;
    DoubleVector _qNorm = qNorm1(ln * time_step);
    for (size_t i=0; i<heatSourceNumber; i++) { _norm += sqr(mq[ln][i]-_qNorm[i]*no_norm); }
    for (ln=1; ln<time_size-1; ln++)
    {
        _qNorm = qNorm1(ln * time_step);
        for (unsigned int i=0; i<heatSourceNumber; i++) { _norm += 0.5*sqr(mq[ln][i]-_qNorm[i]*no_norm); }
    }
    ln = time_size-1;
    _qNorm = qNorm1(ln * time_step);
    for (size_t i=0; i<heatSourceNumber; i++) { _norm += sqr(mq[ln][i]-_qNorm[i]*no_norm); }

    return _norm*time_step;
}

auto p3p2::Functional::gradient(const DoubleVector &x, DoubleVector &g) const -> void
{
    const auto time_size = CommonParameters::timeDimension().size();
    const auto time_step = CommonParameters::timeDimension().step();
    const_cast<p3p2::Functional*>(this)->convert1(x, heatSourceNumber, time_size);

    //IHeatEquationIBVP::implicit_calculate_D1V1();
    //IHeatEquationFBVP::implicit_calculate_D1V1();

    g.resize(x.length(), 0.0);

    for (size_t ln=0; ln<time_size; ln++)
    {
        DoubleVector _qNorm = qNorm1(ln * time_step);
        for (size_t i=0; i<heatSourceNumber; i++)
        {
            g[i*time_size+ln] = -mp[ln][i] + 2.0*epsilon*(mq[ln][i] - _qNorm[i]*no_norm);
        }
    }

    //for (size_t i=0; i<heatSourceNumber; i++) { g[i*time_size] = 0.0; }
}

auto p3p2::Functional::print(unsigned int it, const DoubleVector &x, const DoubleVector &g, double f, double _alpha, GradientBasedMethod::MethodResult result) const -> void
{
    const size_t time_size = CommonParameters::timeDimension().size();
    const auto vector_size = time_size;
    IPrinter::printSeperatorLine();

    const auto _fx = fx(x);
    const auto _integral = integral(x);
    const auto _norm = norm(x);
    printf_s("I[%4d] fx: %10.8f int: %10.8f nrm: %10.8f %10.8f %d\n", it, _fx, _integral, _norm, _alpha, result);
    IPrinter::printVector(_w, _p, x.mid(0*vector_size, 1*vector_size-1), "q1");
    IPrinter::printVector(_w, _p, x.mid(1*vector_size, 2*vector_size-1), "q2");
    IPrinter::printVector(_w, _p, g.mid(0*vector_size, 1*vector_size-1), "g1");
    IPrinter::printVector(_w, _p, g.mid(1*vector_size, 2*vector_size-1), "g2");
    IPrinter::printVector(_w, _p, U, "U:");
    IPrinter::printVector(_w, _p, V, "V:");
}

void p3p2::Functional::Main(int /*argc*/, char **/*argv*/)
{
    double thermalDiffusivity = 0.01, thermalConductivity = 0.00, thermalConvection = 0.00;

    const auto time_min = 0;
    const auto time_max = 200;
    const auto time_size = time_max+1;
    const auto time_step = 0.005;

    const auto dimX_min = 0;
    const auto dimX_max = 100;
    const auto dimX_size = dimX_max+1;
    const auto dimX_step = 0.01;

    const auto vctr_size = time_size;
    const unsigned int start[] = { 0*vctr_size+0, 1*vctr_size+0 };
    const unsigned int finsh[] = { 1*vctr_size-1, 2*vctr_size-1 };

    Dimension _timeDimension(time_step, time_min, time_max);
    Dimension _spaceDimensionX(dimX_step, dimX_min, dimX_max);

    Functional fx(thermalDiffusivity, thermalConductivity, thermalConvection);
    //HeatEquationIBVP1 fx;
    fx.setTimeDimension(_timeDimension);
    fx.setSpaceDimensionX(_spaceDimensionX);

    /******************************************************************
     *                     Calculating V(x) function
     ******************************************************************/
    {
        DoubleVector x(vctr_size*fx.heatSourceNumber);
        for (unsigned int ln=0; ln<vctr_size; ln++)
        {
            double t = ln * time_step;
            DoubleVector _qNorm = fx.qNorm1(t);
            x[0*vctr_size + ln] = 2.0;//_qNorm[0];
            x[1*vctr_size + ln] = 2.0;//_qNorm[1];
        }
        IPrinter::printVector(fx._w, fx._p, x.mid(start[0], finsh[0]), "q1:");
        IPrinter::printVector(fx._w, fx._p, x.mid(start[1], finsh[1]), "q2:");
        fx.convert1(x, fx.heatSourceNumber, vctr_size);
        fx.calculate_forward(x);
        //fx.HeatEquationIBVP1::implicit_calculate_D1V1();
        IPrinter::printVector(fx._w, fx._p, fx.U, "V: ");
        fx.V = fx.U;
    }
    return;

    /******************************************************************
     *                     Calculating V(x) function
     ******************************************************************/
//    DoubleVector x(vctr_size*fx.heatSourceNumber, 2.0);
//    //x[0*vctr_size] = 0.0;
//    //x[1*vctr_size] = 0.0;

//    //    for (unsigned int ln=0; ln<vctr_size; ln++)
//    //    {
//    //        double t = ln * time_step;
//    //        DoubleVector _qNorm = s.qNorm1(t);
//    //        x[0*vctr_size + ln] = _qNorm[0] * 1.1;//( rand()%2==0 ? 0.8 : 1.2 );
//    //        x[1*vctr_size + ln] = _qNorm[1] * 1.1;//( rand()%2==0 ? 0.8 : 1.2 );
//    //    }

//    DoubleVector g(x.length());
//    DoubleVector g1(x.length());

//    puts("---");
//    IPrinter::printVector(fx._w, fx._p, x.mid(start[0], finsh[0]), "q1:");
//    IPrinter::printVector(fx._w, fx._p, x.mid(start[1], finsh[1]), "q2:");

//    puts("---");
//    fx.gradient(x, g);
//    IPrinter::printVector(fx._w, fx._p, g.mid(start[0], finsh[0]).L2Normalize(), "g1:");
//    IPrinter::printVector(fx._w, fx._p, g.mid(start[1], finsh[1]).L2Normalize(), "g2:");

//    puts("---");
//    IGradient::Gradient(&fx, 0.01, x, g1);
//    IPrinter::printVector(fx._w, fx._p, g1.mid(start[0], finsh[0]).L2Normalize(), "g1:");
//    IPrinter::printVector(fx._w, fx._p, g1.mid(start[1], finsh[1]).L2Normalize(), "g2:");

    //    puts("Starting optimization...");
    //    //ConjugateGradient gm;
    //    SteepestDescentGradient gm;
    //    gm.setFunction(&s);
    //    gm.setGradient(&s);
    //    gm.setPrinter(&s);
    //    //gm.setProjection(&s);
    //    //gm.setGradientNormalizer(&prob);
    //    gm.setOptimalityTolerance(0.000001);
    //    gm.setStepTolerance(0.000001);
    //    gm.setFunctionTolerance(0.000001);
    //    gm.setR1MinimizeEpsilon(0.01, 0.001);
    //    //gm.setMaxIterationCount(10);
    //    gm.setNormalize(false);
    //    gm.showExitMessage(true);
    //    gm.calculate(x);

    //    puts("Optimization is finished.");
    //    IPrinter::printVector(s._w, s._p, x.mid(0*vctr_size, 1*vctr_size-1));
    //    IPrinter::printVector(s._w, s._p, x.mid(1*vctr_size, 2*vctr_size-1));
    //    puts("+++");
}
