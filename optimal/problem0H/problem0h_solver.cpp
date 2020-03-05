#include "problem0h_solver.h"
#include <r1minimize.h>

using namespace h0p;

void ProblemSolver::Main(int argc, char **argv)
{
#ifdef USE_LIB_IMAGING
    QGuiApplication app(argc, argv);
#endif

    //checkingForwardProblem();
    //compareGradients();
    optimization();
}

ProblemSolver::ProblemSolver()
{
    forward = new WaveEquationIBVP(this);
    backward = new WaveEquationFBVP(this);
}

ProblemSolver::ProblemSolver(const Dimension &timeDimension, const Dimension &spaceDimensionX, const Dimension &spaceDimensionY)
{
    this->_timeDimension = timeDimension;
    this->_spaceDimensionX = spaceDimensionX;
    this->_spaceDimensionY = spaceDimensionY;

    forward = new WaveEquationIBVP(this);
    backward = new WaveEquationFBVP(this);
}

ProblemSolver::ProblemSolver(const ProblemSolver &)
{
}

ProblemSolver::~ProblemSolver()
{
    delete forward;
    delete backward;
}

void ProblemSolver::checkingForwardProblem()
{
    int N = 100;
    int M = 100;

    ProblemSolver fw1;

    fw1.externalSource = {SpacePoint(0.25, 0.36), 0.05, 0.05, 0.05, 0.01};
    fw1.source_number = 2;

    fw1.setDimension(Dimension(0.01, 0, 2000), Dimension(0.01, 0, N), Dimension(0.01, 0, M));
    fw1.optimalParameters[0].initialize(Dimension(0.01, 0, 2000), Dimension(0.01, 0, N), Dimension(0.01, 0, M)).distribute(SpacePoint(0.308, 0.608));
    fw1.optimalParameters[1].initialize(Dimension(0.01, 0, 2000), Dimension(0.01, 0, N), Dimension(0.01, 0, M)).distribute(SpacePoint(0.708, 0.208));

    fw1.f_saveToFilePng = true;
    fw1.f_saveToFileTxt = true;
    fw1.forward->setWaveSpeed(1.0);
    fw1.forward->setWaveDissipation(0.0);
    fw1.forward->implicit_calculate_D2V1();

    return;

    //    Problem0HForward fw2;
    //    fw2.setTimeDimension(Dimension(0.005, 0, 200));
    //    fw2.addSpaceDimension(Dimension(0.01, 0, 100));
    //    fw2.addSpaceDimension(Dimension(0.01, 0, 100));
    //    fw2.source_number = 0;
    //    fw2.ksi = SpacePoint(0.50, 0.50);

    //    fw2.p_sigmaX = 0.01;
    //    fw2.p_sigmaY = 0.01;
    //    fw2.p_sigmaT = 0.01;

    //    DoubleMatrix u2;
    //    fw2.implicit_calculate_D2V1(u2, 1.0, 0.0);
    //    IPrinter::printMatrix(u2);
    //    IPrinter::printSeperatorLine();
}

void ProblemSolver::compareGradients()
{
    ProblemSolver functional;

    functional.externalSource = {SpacePoint(0.25, 0.36), 0.05, 0.05, 0.05, 0.01};

    functional.forward->setWaveSpeed(1.0);
    functional.forward->setWaveDissipation(0.00);
    functional.backward->setWaveSpeed(1.0);
    functional.backward->setWaveDissipation(0.00);

    functional.eps1 = 1.0;
    functional.esp2 = 1.0;
    functional.source_number = 2;
    functional.setDimension(Dimension(0.01, 0, 100), Dimension(0.01, 0, 100), Dimension(0.01, 0, 100));
    //functional.optimalParameters[0].distribute(SpacePoint(0.314, 0.647));
    //functional.optimalParameters[1].distribute(SpacePoint(0.729, 0.232));

    functional.optimalParameters[0].distribute(SpacePoint(0.310, 0.640));
    functional.optimalParameters[1].distribute(SpacePoint(0.720, 0.230));

    for (unsigned int n=0; n<201; n++) functional.optimalParameters[0].pwr_vl[n] = 0.1;//0.001*n;
    for (unsigned int n=0; n<201; n++) functional.optimalParameters[1].pwr_vl[n] = 0.2;//0.002*n;

    DoubleVector x;
    functional.parameterToVector(x);

    const unsigned int L = functional.timeDimension().size()-1;
    const unsigned int s1 = 0*(2*L+1); const unsigned int f1 = s1+2*L;
    const unsigned int s2 = 1*(2*L+1); const unsigned int f2 = s2+2*L;
    const unsigned int s3 = 2*(2*L+1); const unsigned int f3 = s3+3;

    // printing parameter values
    IPrinter::printSeperatorLine();
    IPrinter::printVector(x.mid(s1, f1));
    IPrinter::printVector(x.mid(s2, f2));
    IPrinter::print(x.mid(s3, f3), 4);
    IPrinter::printSeperatorLine();

    // printing analitic gradients
    DoubleVector ga;
    functional.gradient(x, ga);
    IPrinter::printVector(ga.mid(s1, f1).EuclideanNormalize());
    IPrinter::printVector(ga.mid(s2, f2).EuclideanNormalize());
    IPrinter::print(ga.mid(s3, f3).EuclideanNormalize(), 4);
    IPrinter::printSeperatorLine();

    // printing numerical gradients
    DoubleVector gn;
    gn.resize(x.length());
    IGradient::Gradient(&functional, 0.01, x, gn, s1, f1);
    IGradient::Gradient(&functional, 0.01, x, gn, s2, f2);
    IGradient::Gradient(&functional, 0.01, x, gn, s3, f3);
    gn[s1] = gn[s2] = 0.0;
    gn[f1] = gn[f2] = 0.0;
    IPrinter::printVector(gn.mid(s1, f1).EuclideanNormalize());
    IPrinter::printVector(gn.mid(s2, f2).EuclideanNormalize());
    IPrinter::print(gn.mid(s3, f3).EuclideanNormalize(), 4);
    IPrinter::printSeperatorLine();

    //IGradient::Gradient(&functional, 0.001, x, gn, s1, f1);
    //IGradient::Gradient(&functional, 0.001, x, gn, s2, f2);
    IGradient::Gradient(&functional, 0.001, x, gn, s3, f3);
    IPrinter::printVector(gn.mid(s1, f1).EuclideanNormalize());
    IPrinter::printVector(gn.mid(s2, f2).EuclideanNormalize());
    IPrinter::print(gn.mid(s3, f3).EuclideanNormalize(), 4);
    IPrinter::printSeperatorLine();
}

void ProblemSolver::optimization()
{
    (new ProblemSolver)->fx(4.0);

//    double a,b,fxa,fxb;
//    bool unimodal;
//    R1FxMinimizer r1min;
//    r1min.setFunction(new ProblemSolver);
//    //r1min.straightLineSearch(0.5, 0.1, a, b, fxa, fxb, unimodal);
//    double t = 1.0;
//    a = 1.0;
//    b = 21.0;
//    r1min.uniformLineSearch(t, a, b, 100);
}

//--------------------------------------------------------------------------------------------------------------//

auto ProblemSolver::setDimension(const Dimension &timeDimension, const Dimension &dimensionX, const Dimension &dimensionY) -> void
{
    setTimeDimension(timeDimension);
    setSpaceDimensionX(dimensionX);
    setSpaceDimensionY(dimensionY);

    //const unsigned int L = static_cast<unsigned int> ( timeDimension.size() );
    //const unsigned int N = static_cast<unsigned int> ( dimensionX.size() );
    //const unsigned int M = static_cast<unsigned int> ( dimensionY.size() );

    u1.resize(dimensionY.size(), dimensionX.size(), 0.0);
    u2.resize(dimensionY.size(), dimensionX.size(), 0.0);

    optimalParameters.resize(source_number);

    for (unsigned int sn=0; sn<source_number; sn++) optimalParameters[sn].initialize(timeDimension, dimensionX, dimensionY);
}

auto ProblemSolver::fx(const DoubleVector &x) const -> double
{
    vectorToParameter(x);

    forward->implicit_calculate_D2V1();
    double sum = eps1 * integral1(u1) + esp2 * integral2(u2);
    //sum += norm();
    //sum += penalty();
    return sum;
}

auto ProblemSolver::fx(double t) const -> double
{
    IPrinter::printSeperatorLine(nullptr, '=');
    std::cout << "Time: " << t << std::endl;

    ProblemSolver fw1;

    fw1.externalSource = {SpacePoint(0.25, 0.36), 0.10, 0.05, 0.05, 0.01};

    fw1.source_number = 2;

    fw1.waveSpeed = 1.0;
    fw1.waveDissapation = 0.0;

    fw1.eps1 = 1.0;
    fw1.esp2 = 1.0;

    fw1.setDimension(Dimension(0.01, 0, static_cast<int>(t/0.01)), Dimension(0.01, 0, 100), Dimension(0.01, 0, 100));
    fw1.optimalParameters[0].distribute(SpacePoint(0.38, 0.68));
    fw1.optimalParameters[1].distribute(SpacePoint(0.78, 0.28));

    DoubleVector x;
    fw1.parameterToVector(x);

    ConjugateGradient g;
    //SteepestDescentGradient g;
    g.setFunction(&fw1);
    g.setGradient(&fw1);
    g.setPrinter(&fw1);
    g.setProjection(&fw1);
    g.setOptimalityTolerance(0.0);
    g.setFunctionTolerance(0.0);
    g.setStepTolerance(0.0);
    g.setR1MinimizeEpsilon(0.1, 0.01);
    g.setMaxIterationCount(100);
    g.setNormalize(false);
    g.showExitMessage(true);
    g.calculate(x);

//    fw1.vectorToParameter(x);
//    fw1.f_saveToFileTxt = true;
//    fw1.forward->implicit_calculate_D2V1();

//    unsigned int L = 100;
//    unsigned int s1 = 0*(2*L+1); unsigned int f1 = 0*(2*L+1)+2*L;
//    unsigned int s2 = 1*(2*L+1); unsigned int f2 = 1*(2*L+1)+2*L;
//    unsigned int s3 = 2*(2*L+1); unsigned int f3 = 2*(2*L+1)+3;

//    IPrinter::printVector(x.mid(s1, f1).EuclideanNormalize());
//    IPrinter::printVector(x.mid(s2, f2).EuclideanNormalize());
//    IPrinter::print(x.mid(s3, f3).EuclideanNormalize(), 4);
//    IPrinter::printSeperatorLine();

    double res = fw1.fx(x);
    std::cout << "Result: " << res << std::endl;

    return res;
}

auto ProblemSolver::integral1(const DoubleMatrix &) const -> double
{
    //const Dimension &time = Problem0HForward::timeDimension();
    const Dimension &dimX = spaceDimensionX();
    const Dimension &dimY = spaceDimensionY();
    //const unsigned int L = static_cast<unsigned int> ( time.size() );
    const unsigned int N = static_cast<unsigned int> ( dimX.size()-1 );
    const unsigned int M = static_cast<unsigned int> ( dimY.size()-1 );
    //const double ht = time.step();
    const double hx = dimX.step();
    const double hy = dimY.step();

    double udiff = 0.0;
    double usum = 0.0;

    udiff = u1[0][0]; usum += 0.25 * udiff * udiff;// * mu(0, 0);
    udiff = u1[0][N]; usum += 0.25 * udiff * udiff;// * mu(N, 0);
    udiff = u1[M][0]; usum += 0.25 * udiff * udiff;// * mu(0, M);
    udiff = u1[M][N]; usum += 0.25 * udiff * udiff;// * mu(N, M);

    for (unsigned int n=1; n<=N-1; n++)
    {
        udiff = u1[0][n]; usum += 0.5 * udiff * udiff;// * mu(n, 0);
        udiff = u1[M][n]; usum += 0.5 * udiff * udiff;// * mu(n, M);
    }

    for (unsigned int m=1; m<=M-1; m++)
    {
        udiff = u1[m][0]; usum += 0.5 * udiff * udiff;// * mu(0, m);
        udiff = u1[m][N]; usum += 0.5 * udiff * udiff;// * mu(N, m);
    }

    for (unsigned int m=1; m<=M-1; m++)
    {
        for (unsigned int n=1; n<=N-1; n++)
        {
            udiff = u1[m][n]; usum += udiff * udiff;// * mu(n, m);
        }
    }

    return usum*(hx*hy);
}

auto ProblemSolver::integral2(const DoubleMatrix &) const -> double
{
    //const Dimension &time = Problem0HForward::timeDimension();
    const Dimension &dimX = spaceDimensionX();
    const Dimension &dimY = spaceDimensionY();
    //const unsigned int L = static_cast<unsigned int> ( time.size() );
    const unsigned int N = static_cast<unsigned int> ( dimX.size()-1 );
    const unsigned int M = static_cast<unsigned int> ( dimY.size()-1 );
    //const double ht = time.step();
    const double hx = dimX.step();
    const double hy = dimY.step();

    double udiff = 0.0;
    double usum = 0.0;

    udiff = u2[0][0]; usum += 0.25 * udiff * udiff;// * mu(0, 0);
    udiff = u2[0][N]; usum += 0.25 * udiff * udiff;// * mu(N, 0);
    udiff = u2[M][0]; usum += 0.25 * udiff * udiff;// * mu(0, M);
    udiff = u2[M][N]; usum += 0.25 * udiff * udiff;// * mu(N, M);

    for (unsigned int n=1; n<=N-1; n++)
    {
        udiff = u2[0][n]; usum += 0.5 * udiff * udiff;// * mu(n, 0);
        udiff = u2[M][n]; usum += 0.5 * udiff * udiff;// * mu(n, M);
    }

    for (unsigned int m=1; m<=M-1; m++)
    {
        udiff = u2[m][0]; usum += 0.5 * udiff * udiff;// * mu(0, m);
        udiff = u2[m][N]; usum += 0.5 * udiff * udiff;// * mu(N, m);
    }

    for (unsigned int m=1; m<=M-1; m++)
    {
        for (unsigned int n=1; n<=N-1; n++)
        {
            udiff = u2[m][n]; usum += udiff * udiff;// * mu(n, m);
        }
    }

    return usum*(hx*hy);
}

auto ProblemSolver::norm() const -> double { return 0.0; }

auto ProblemSolver::penalty() const -> double { return 0.0; }

auto ProblemSolver::gradient(const DoubleVector &x, DoubleVector &g) const -> void
{
    vectorToParameter(x);

    g.clear(); g.resize(x.length());

    const Dimension &time = timeDimension();
    //const Dimension &dimX = Problem0HForward::spaceDimension(Dimension::DimensionX);
    //const Dimension &dimY = Problem0HForward::spaceDimension(Dimension::DimensionY);
    const unsigned int L = static_cast<unsigned int>(time.size()-1);
    //const unsigned int N = static_cast<unsigned int>(dimX.size());
    //const unsigned int M = static_cast<unsigned int>(dimY.size());
    //const double hx = dimX.step();
    //const double hy = dimY.step();
    const double ht = time.step();

    //const_cast<ProblemSolver*>(this)->_timeDimension.setMax(_timeDimension.max()+1);
    forward->implicit_calculate_D2V1();
    //const_cast<ProblemSolver*>(this)->_timeDimension.setMax(_timeDimension.max()-1);
    backward->implicit_calculate_D2V1();

    unsigned int length = 2*L+1;
    for (unsigned int sn=0; sn<source_number; sn++)
    {
        const Problem0HParameter &optimalParameter = optimalParameters[sn];

        unsigned int offset = sn*length;
        for (unsigned int ln=0; ln<=2*L; ln++)
        {
            g[offset+ln] = -optimalParameter.psi_vl[ln];
        }
        g[offset+0] = 0.0;
        g[offset+2*L] = 0.0;

        //---------------------------------------------------------------------------//

        const unsigned int point_offset = source_number*length;
        const unsigned int ix = point_offset+2*sn+0;
        const unsigned int iy = point_offset+2*sn+1;

        g[ix] = 0.0;
        g[iy] = 0.0;
        unsigned int ln = 0;
        g[ix] += 0.5*optimalParameter.psi_dx[ln] * optimalParameter.pwr_vl[ln];
        g[iy] += 0.5*optimalParameter.psi_dy[ln] * optimalParameter.pwr_vl[ln];
        for (ln=2; ln<=2*(L-1); ln+=2)
        {
            g[ix] += optimalParameter.psi_dx[ln] * optimalParameter.pwr_vl[ln];
            g[iy] += optimalParameter.psi_dy[ln] * optimalParameter.pwr_vl[ln];
        }
        ln = 2*L;
        g[ix] += 0.5*optimalParameter.psi_dx[ln] * optimalParameter.pwr_vl[ln];
        g[iy] += 0.5*optimalParameter.psi_dy[ln] * optimalParameter.pwr_vl[ln];

        g[ix] *= -ht;
        g[iy] *= -ht;
    }
}

auto ProblemSolver::project(DoubleVector &x, unsigned int index) -> void
{
    const Dimension &time = timeDimension();
    const unsigned int L = static_cast<unsigned int>(time.size()-1);
    if (index >= (2*L+1)*source_number)
    {
        if (x[index] < 0.05) x[index] = 0.05;
        if (x[index] > 0.95) x[index] = 0.95;
    }
}

auto ProblemSolver::project(DoubleVector &) const  -> void {}

auto ProblemSolver::print(unsigned int i, const DoubleVector &x, const DoubleVector &g, double f,
                          double alpha, GradientMethod::MethodResult result) const -> void
{
    C_UNUSED(i); C_UNUSED(x); C_UNUSED(g); C_UNUSED(f); C_UNUSED(alpha); C_UNUSED(result);
    const char* msg = nullptr; C_UNUSED(msg);
    if (result == GradientMethod::MethodResult::BREAK_FIRST_ITERATION)    msg = "BREAK_FIRST_ITERATION   ";
    if (result == GradientMethod::MethodResult::FIRST_ITERATION)          msg = "FIRST_ITERATION         ";
    if (result == GradientMethod::MethodResult::BREAK_GRADIENT_NORM_LESS) msg = "BREAK_GRADIENT_NORM_LESS";
    if (result == GradientMethod::MethodResult::BREAK_DISTANCE_LESS)      msg = "BREAK_DISTANCE_LESS     ";
    if (result == GradientMethod::MethodResult::NEXT_ITERATION)           msg = "NEXT_ITERATION          ";

    ProblemSolver* fw = const_cast<ProblemSolver*>(this);
    fw->vectorToParameter(x);
    //const_cast<ProblemSolver*>(this)->_timeDimension.setMax(_timeDimension.max()+1);
    forward->implicit_calculate_D2V1();
    //const_cast<ProblemSolver*>(this)->_timeDimension.setMax(_timeDimension.max()-1);

    const Dimension &time = timeDimension();
    const unsigned int L = static_cast<unsigned int>(time.size()-1);
    const unsigned int o = (2*L+1)*source_number;

    printf("I[%3d]: F:%.6f I:%.6f P:%.6f N:%.5f R:%.3f e:%.3f a:%10.6f | ", i, f, fx(x), 0.0, 0.0, 0.0, 0.0, alpha);
    printf("%12.8f %12.8f | %12.8f %12.8f | ", u1.min(), u1.max(), u2.min(), u2.max());
    printf("eta: %8.6f %8.6f %8.6f %8.6f\n", x[o+0], x[o+1], x[o+2], x[o+3]);
    //IPrinter::printSeperatorLine("-");
}

//--------------------------------------------------------------------------------------------------------------//

WaveEquationIBVP::WaveEquationIBVP(ProblemSolver *solver) { this->solver = solver; }

WaveEquationIBVP::~WaveEquationIBVP() {}

double WaveEquationIBVP::initial(const SpaceNodePDE &sn, InitialCondition cn) const { return solver->frw_initial(sn, cn); }

double WaveEquationIBVP::boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &cn) const { return solver->frw_boundary(sn, tn, cn); }

double WaveEquationIBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const { return solver->frw_f(sn, tn); }

void WaveEquationIBVP::layerInfo(const DoubleMatrix &u, const TimeNodePDE &tn) const { return solver->frw_layerInfo(u, tn); }

Dimension WaveEquationIBVP::timeDimension() const
{
//    return Dimension(solver->timeDimension().step(), solver->timeDimension().min(), solver->timeDimension().max());
    return solver->timeDimension();
}

Dimension WaveEquationIBVP::spaceDimensionX() const { return solver->spaceDimensionX(); }

Dimension WaveEquationIBVP::spaceDimensionY() const { return solver->spaceDimensionY(); }

Dimension WaveEquationIBVP::spaceDimensionZ() const { return solver->spaceDimensionZ(); }

double ProblemSolver::frw_initial(const SpaceNodePDE &, InitialCondition) const { return 0.0; }

double ProblemSolver::frw_boundary(const SpaceNodePDE &, const TimeNodePDE &, BoundaryConditionPDE & condition) const
{
    condition = BoundaryConditionPDE(BoundaryCondition::Dirichlet, 1.0, 0.0, 0.0);
    return 0.0;
}

double ProblemSolver::frw_f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    double pv = p(sn, tn);
    //return pv;
    unsigned int ln = static_cast<unsigned int>(tn.i);

    double pulse1 = optimalParameters[0].pwr_vl[ln] * optimalParameters[0].deltaGrid.weight(sn);
    double pulse2 = optimalParameters[1].pwr_vl[ln] * optimalParameters[1].deltaGrid.weight(sn);

    return pv + pulse1 + pulse2;
}

void ProblemSolver::frw_layerInfo(const DoubleMatrix &u, const TimeNodePDE &tn) const
{
    frw_calculateU1U2(u, tn);

    //unsigned int w = 6;
    //unsigned int p = 3;
    //if (abs(tn.t - 0.1) <= 0.0) { IPrinter::printMatrix(w,p,u); IPrinter::printSeperatorLine(); }
    //if (abs(tn.t - 0.2) <= 0.0) { IPrinter::printMatrix(w,p,u); IPrinter::printSeperatorLine(); }
    //if (abs(tn.t - 0.3) <= 0.0) { IPrinter::printMatrix(w,p,u); IPrinter::printSeperatorLine(); }
    //if (abs(tn.t - 0.4) <= 0.0) { IPrinter::printMatrix(w,p,u); IPrinter::printSeperatorLine(); }
    //if (abs(tn.t - 0.5) <= 0.0) { IPrinter::printMatrix(w,p,u); IPrinter::printSeperatorLine(); }
    //if (abs(tn.t - 0.6) <= 0.0) { IPrinter::printMatrix(w,p,u); IPrinter::printSeperatorLine(); }
    //if (abs(tn.t - 0.7) <= 0.0) { IPrinter::printMatrix(w,p,u); IPrinter::printSeperatorLine(); }
    //if (abs(tn.t - 0.8) <= 0.0) { IPrinter::printMatrix(w,p,u); IPrinter::printSeperatorLine(); }
    //if (abs(tn.t - 0.9) <= 0.0) { IPrinter::printMatrix(w,p,u); IPrinter::printSeperatorLine(); }
    //if (abs(tn.t - 1.0) <= 0.0) { IPrinter::printMatrix(w,p,u); IPrinter::printSeperatorLine(); }

#ifdef USE_LIB_XLSX_WRITER
    frw_saveToExcel(u, tn);
#endif

    if (f_saveToFilePng) frw_saveToImage(u, tn);
    if (f_saveToFileTxt) frw_saveToTextF(u, tn);
}

double ProblemSolver::p(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    if (fabs(tn.i - 10.0*externalSource.sigmaT) == 0.0) return 0.0;

    //const double ht = timeDimension().step();
    const double factor1 = 1.0/(2.0*M_PI*externalSource.sigmaX*externalSource.sigmaY);
    const double sigmax2 = 1.0/(2.0*externalSource.sigmaX*externalSource.sigmaX);
    const double sigmay2 = 1.0/(2.0*externalSource.sigmaY*externalSource.sigmaY);
    const double factor2 = 2.0/(sqrt(2.0*M_PI)*externalSource.sigmaT);
    const double sigmat2 = 1.0/(2.0*externalSource.sigmaT*externalSource.sigmaT);

    const double power = externalSource.power;
    const SpacePoint pnt = externalSource.point;

    double a, b;
    a = factor1 * exp(-(sigmax2*(sn.x-pnt.x)*(sn.x-pnt.x)+sigmay2*(sn.y-pnt.y)*(sn.y-pnt.y)));
    b = factor2 * exp(-(sigmat2*(tn.t-0.0)*(tn.t-0.0)));

    return power*a*b;
}

void ProblemSolver::frw_calculateU1U2(const DoubleMatrix &u, const TimeNodePDE &tn) const
{
    const Dimension &time = forward->timeDimension();

    const unsigned int L = static_cast<unsigned int>(time.size()-1);
    const double ht = time.step();

    //printf("%d %d\n", L, time.size());

    const Dimension &dimX = spaceDimensionX();
    const Dimension &dimY = spaceDimensionY();
    const unsigned int N = static_cast<unsigned int>(dimX.size()-1);
    const unsigned int M = static_cast<unsigned int>(dimY.size()-1);

//    if (tn.i == 2*(L-2))
//    {
//        for (unsigned int m=0; m<=M; m++)
//        {
//            for (unsigned int n=0; n<=N; n++)
//            {
//                u2[m][n] = -u[m][n];
//            }
//        }
//    }

//    if (tn.i == 2*L)
//    {
//        for (unsigned int m=0; m<=M; m++)
//        {
//            for (unsigned int n=0; n<=N; n++)
//            {
//                u1[m][n]  = u[m][n];
//                u2[m][n] += u[m][n];
//                u2[m][n] /= (2.0*ht);
//            }
//        }
//    }


    if (tn.i == 2*(L-2))
    {
        //printf("1 %d %f\n", tn.i, tn.t);
        for (unsigned int m=0; m<=M; m++)
        {
            for (unsigned int n=0; n<=N; n++)
            {
                u2[m][n] = u[m][n];
            }
        }
    }

    if (tn.i == 2*(L-1))
    {
        //printf("2 %d %f\n", tn.i, tn.t);
        for (unsigned int m=0; m<=M; m++)
        {
            for (unsigned int n=0; n<=N; n++)
            {
                u2[m][n] += -4.0*u[m][n];
            }
        }
    }

    if (tn.i == 2*(L-0))
    {
        //printf("3 %d %f\n", tn.i, tn.t);
        for (unsigned int m=0; m<=M; m++)
        {
            for (unsigned int n=0; n<=N; n++)
            {
                u1[m][n]  = u[m][n];
                u2[m][n] += 3.0*u[m][n];
                u2[m][n] /= (2.0*ht);
            }
        }
    }

    //    if (tn.i == 2*(L-2))
    //    {
    //        for (unsigned int m=0; m<=M; m++)
    //        {
    //            for (unsigned int n=0; n<=N; n++)
    //            {
    //                u2[m][n] = 0.0;
    //                u2[m][n] -= u[m][n];
    //            }
    //        }
    //    }

    //    if (tn.i == 2*(L-1))
    //    {
    //        for (unsigned int m=0; m<=M; m++)
    //        {
    //            for (unsigned int n=0; n<=N; n++)
    //            {
    //                u1[m][n] = u[m][n];
    //            }
    //        }
    //    }

    //    if (tn.i == 2*(L-0))
    //    {
    //        for (unsigned int m=0; m<=M; m++)
    //        {
    //            for (unsigned int n=0; n<=N; n++)
    //            {
    //                u2[m][n] += u[m][n];
    //                u2[m][n] /= (2.0*ht);
    //            }
    //        }
    //    }
}

void ProblemSolver::frw_saveToExcel(const DoubleMatrix &u UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM) const
{
#ifdef USE_LIB_XLSX_WRITER
    //    const Dimension &time = Problem0HForward::timeDimension();
    //    const Dimension &dimX = Problem0HForward::spaceDimension(Dimension::DimensionX);
    //    const Dimension &dimY = Problem0HForward::spaceDimension(Dimension::DimensionY);
    //    const unsigned int L = static_cast<unsigned int>(time.size());
    //    const unsigned int N = static_cast<unsigned int>(dimX.size());
    //    const unsigned int M = static_cast<unsigned int>(dimY.size());
    //    const double ht = time.step();
    //    const double hx = dimX.step();
    //    const double hy = dimY.step();

    //    static lxw_workbook *workbook = nullptr;
    //    static lxw_format *format = nullptr;
    //    if (ln == 0)
    //    {
    //        workbook = new_workbook("d:\\chart.xlsx");
    //        format = workbook_add_format(workbook);
    //        format_set_num_format(format, "0.0000000000");
    //    }
    //    if (ln%20==0)
    //    {
    //        std::string name = std::string("layer_")+std::to_string(ln/20);
    //        lxw_worksheet *worksheet = workbook_add_worksheet(workbook, name.data());
    //        for (unsigned int row=0; row<u.rows(); row++)
    //        {
    //            for (unsigned int col=0; col<u.cols(); col++)
    //            {
    //                std::string line = "=Sheet1!$C$"+std::to_string(row+1)+":$C$"+std::to_string(row+1);
    //                worksheet_write_number(worksheet, static_cast<lxw_row_t>(row), static_cast<lxw_col_t>(col), u[row][col], format);
    //            }
    //        }
    //    }
    //    if (ln == 2*L) { workbook_close(workbook); }

    //    if (ln == 2*(L-0))
    //    {
    //        lxw_workbook  *workbook  = new_workbook("d:\\chart.xlsx");
    //        lxw_worksheet *worksheet1 = workbook_add_worksheet(workbook, nullptr);
    //        lxw_chartsheet *chartsheet = workbook_add_chartsheet(workbook, nullptr);
    //        lxw_chart *chart = workbook_add_chart(workbook, LXW_CHART_AREA);
    //        lxw_format *format = workbook_add_format(workbook);
    //        format_set_num_format(format, "0.0000000000");
    //        for (unsigned int row=0; row<u.rows(); row++)
    //        {
    //            for (unsigned int col=0; col<u.cols(); col++)
    //            {
    //                // =Sheet1!$A$1:$A$6
    //                // =Sheet1!$A$1:$CW$1
    //                std::string line = "=Sheet1!$A$"+std::to_string(row+1)+":$CW$"+std::to_string(row+1);
    //                worksheet_write_number(worksheet1, static_cast<lxw_row_t>(row), static_cast<lxw_col_t>(col), u[row][col], format);
    //                //chart_add_series(chart, nullptr, line.data());
    //            }
    //        }
    //        chartsheet_set_chart(chartsheet, chart);
    //        workbook_close(workbook);
    //    }
#endif
}

void ProblemSolver::frw_saveToImage(const DoubleMatrix &u UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM) const
{
#ifdef USE_LIB_IMAGING
    //const Dimension &time = Problem0HForward::timeDimension();
    const Dimension &dimX = spaceDimensionX();
    const Dimension &dimY = spaceDimensionY();
    //const unsigned int L = static_cast<unsigned int>(time.size());
    const unsigned int N = static_cast<unsigned int>(dimX.size());
    const unsigned int M = static_cast<unsigned int>(dimY.size());
    //const double ht = time.step();
    //const double hx = dimX.step();
    //const double hy = dimY.step();

    //QDir path("D:/data2");
    //if (!path.exists()) path.mkdir("D:/data2");

    QString filename = QString("data/problem0H/f/png/f_%1.png").arg(tn.i, 4, 10, QChar('0'));
    QPixmap pixmap;
    visualizeMatrixHeat(u, u.min(), u.max(), pixmap, N+1, M+1);
    pixmap.save(filename);
#endif
}

void ProblemSolver::frw_saveToTextF(const DoubleMatrix &u UNUSED_PARAM, const TimeNodePDE &tn) const
{
    static double MIN = +100000.0;
    static double MAX = -100000.0;

    std::string txt_number = std::to_string(tn.i);
    std::string filename = std::string("data/problem0H/f/txt/f_") + std::string(4 - txt_number.length(), '0') + txt_number + std::string(".txt");
    IPrinter::print(u, filename.c_str());
    if (MIN > u.min()) MIN = u.min();
    if (MAX < u.max()) MAX = u.max();
    printf("Forward: %4d %6.3f | %12.8f %12.8f | %12.8f %12.8f | %4d %4d\n", tn.i, tn.t, u.min(), u.max(), MIN, MAX, 0, 0);
}

//--------------------------------------------------------------------------------------------------------------//

WaveEquationFBVP::WaveEquationFBVP(ProblemSolver *solver) { this->solver = solver; }

WaveEquationFBVP::~WaveEquationFBVP() {}

double WaveEquationFBVP::final(const SpaceNodePDE &sn, FinalCondition condition) const { return solver->bcw_final(sn, condition); }

double WaveEquationFBVP::boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const { return solver->bcw_boundary(sn, tn, condition); }

double WaveEquationFBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const { return solver->bcw_f(sn, tn); }

void WaveEquationFBVP::layerInfo(const DoubleMatrix &p, const TimeNodePDE &tn) const { return solver->bcw_layerInfo(p, tn); }

Dimension WaveEquationFBVP::timeDimension() const { return solver->timeDimension(); }

Dimension WaveEquationFBVP::spaceDimensionX() const { return solver->spaceDimensionX(); }

Dimension WaveEquationFBVP::spaceDimensionY() const { return solver->spaceDimensionY(); }

Dimension WaveEquationFBVP::spaceDimensionZ() const { return solver->spaceDimensionZ(); }

double ProblemSolver::bcw_final(const SpaceNodePDE &sn, FinalCondition condition) const
{
    unsigned int n = static_cast<unsigned int>(sn.i);
    unsigned int m = static_cast<unsigned int>(sn.j);
    if (condition == FinalCondition::FinalValue) { return -2.0*esp2*(u2[m][n]/*-U2[m][n]*/); }
    if (condition == FinalCondition::FinalFirstDerivative) { return +2.0*eps1*(u1[m][n]/*-U1[m][n]*/) + waveDissapation*bcw_final(sn, FinalCondition::FinalValue); }
    throw std::exception();
}

double ProblemSolver::bcw_boundary(const SpaceNodePDE&, const TimeNodePDE&, BoundaryConditionPDE &condition) const
{
    condition = BoundaryConditionPDE(BoundaryCondition::Dirichlet, 1.0, 0.0, 0.0);
    return 0.0;
}

double ProblemSolver::bcw_f(const SpaceNodePDE&, const TimeNodePDE&) const { return 0.0; }

void ProblemSolver::bcw_layerInfo(const DoubleMatrix &p, const TimeNodePDE &tn) const
{
    bcw_saveBackwardInformarion(p, tn);
    //saveToImage(p, ln);
}

void ProblemSolver::bcw_saveBackwardInformarion(const DoubleMatrix &p, const TimeNodePDE &tn) const
{
    ProblemSolver* const_this = const_cast<ProblemSolver*>(this);

    for (unsigned int sn=0; sn<source_number; sn++)
    {
        double psi_vl, psi_dx, psi_dy;
        psi_vl = optimalParameters[sn].deltaGrid.lumpPointGauss(p, psi_dx, psi_dy);

        const_this->optimalParameters[sn].psi_vl[tn.i] = psi_vl;
        const_this->optimalParameters[sn].psi_dx[tn.i] = psi_dx;
        const_this->optimalParameters[sn].psi_dy[tn.i] = psi_dy;
    }
}

void ProblemSolver::bcw_saveToImage(const DoubleMatrix &p, const TimeNodePDE &tn) const
{
#ifdef USE_LIB_IMAGING
    //const Dimension &time = Problem0HBckward::timeDimension();
    const Dimension &dimX = ProblemSolver::spaceDimensionX();
    const Dimension &dimY = ProblemSolver::spaceDimensionY();
    //const unsigned int L = static_cast<unsigned int>(time.size());
    const unsigned int N = static_cast<unsigned int>(dimX.size());
    const unsigned int M = static_cast<unsigned int>(dimY.size());
    //const double ht = time.step();
    //const double hx = dimX.step();
    //const double hy = dimY.step();

    //QString filename1 = QString("data/problem0H/b/txt/b_%1.txt").arg(ln, 4, 10, QChar('0'));
    //IPrinter::print(p,filename1.toLatin1().data());
    //IPrinter::printSeperatorLine();
    //printf("Forward: %4d %0.3f %10.8f %10.8f %4d %4d\n", ln, ln*0.005, u.min(), u.max(), 0, 0);

    QString filename2 = QString("data/problem0H/b/png/b_%1.png").arg(tn.i, 4, 10, QChar('0'));
    QPixmap pixmap;
    visualizeMatrixHeat(p, p.min(), p.max(), pixmap, N+1, M+1);
    pixmap.save(filename2);
    //IPrinter::printSeperatorLine();
#endif
}

//--------------------------------------------------------------------------------------------------------------//

//--------------------------------------------------------------------------------------------------------------//

auto ProblemSolver::vectorToParameter(const DoubleVector &x) const -> void
{
    const Dimension &time = timeDimension();
    const Dimension &dimX = spaceDimensionX();
    const Dimension &dimY = spaceDimensionY();
    const unsigned int L = static_cast<unsigned int> ( time.size() );
    //const unsigned int N = static_cast<unsigned int> ( dimX.size() );
    //const unsigned int M = static_cast<unsigned int> ( dimY.size() );

    ProblemSolver* const_this = const_cast<ProblemSolver*>(this);
    const_this->optimalParameters.resize(source_number);

    const unsigned int length = 2*L-1;
    const unsigned int points_offset = source_number*length;
    for (unsigned int sn=0; sn<source_number; sn++)
    {
        Problem0HParameter &parameter = const_this->optimalParameters[sn];
        parameter.initialize(time, dimX, dimY);
        for (unsigned int ln=0; ln<length; ln++) parameter.pwr_vl[ln] = x.at(sn*length+ln);
        parameter.p.x = x.at(points_offset+2*sn+0);
        parameter.p.y = x.at(points_offset+2*sn+1);
        parameter.distribute(parameter.p);
    }
}

auto ProblemSolver::parameterToVector(DoubleVector &x) const -> void
{
    const Dimension &time = timeDimension();
    //const Dimension &dimX = Problem0HForward::spaceDimension(Dimension::DimensionX);
    //const Dimension &dimY = Problem0HForward::spaceDimension(Dimension::DimensionY);
    const unsigned int L = static_cast<unsigned int> ( time.size() );
    //const unsigned int N = static_cast<unsigned int> ( dimX.size() );
    //const unsigned int M = static_cast<unsigned int> ( dimY.size() );
    //const double ht = time.step();
    //const double hx = dimX.step();
    //const double hy = dimY.step();

    const unsigned int length = 2*L-1;
    const unsigned int size = (length+2)*source_number;
    const unsigned int points_offset = source_number*length;

    x.resize(size);
    for (unsigned int sn=0; sn<source_number; sn++)
    {
        Problem0HParameter &parameter = const_cast<ProblemSolver*>(this)->optimalParameters[sn];
        for (unsigned int ln=0; ln<length; ln++) x[sn*length+ln] = parameter.pwr_vl[ln];
        //x[sn*length+0] = 1.0;
        //x[sn*length+1] = 2.0;
        //x[sn*length+(length-1)] = 0.0;
        x[points_offset+2*sn+0] = parameter.p.x;
        x[points_offset+2*sn+1] = parameter.p.y;
    }

    //printf("parameterToVector << %d %d %d %d %d\n", x.length(), source_number, time.size(), length, size);
    //IPrinter::printSeperatorLine();
    //IPrinter::print(x, x.length());
}
