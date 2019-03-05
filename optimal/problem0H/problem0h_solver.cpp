#include "problem0h_solver.h"

void Problem0HFunctional::Main(int argc, char **argv)
{
    QGuiApplication app(argc, argv);

    Problem0HFunctional functional;
    functional.ksi = SpacePoint(0.25, 0.25);
    functional.a = 1.0;
    functional.gamma = 0.0;
    functional.alpha0 = 1.0;
    functional.setDimension(Dimension(0.005, 0, 200),
                            Dimension(0.010, 0, 100),
                            Dimension(0.010, 0, 100));

    DoubleVector x;
    x.resize(806);

    for (unsigned int i=0; i<=400; i++)
    {
        x[i+000] = 1.0;
        x[i+401] = 1.0;
    }
    x[802+0] = 0.35; x[802+1] = 0.65;
    x[802+2] = 0.75; x[802+3] = 0.25;

    functional.vectorToParameter(x);

    DoubleVector ga;
    functional.gradient(x, ga);
    ga.EuclideanNormalize();
    IPrinter::printVector(ga.mid(0, 400).EuclideanNormalize());
    IPrinter::printVector(ga.mid(401, 801).EuclideanNormalize());
    IPrinter::print(ga.mid(802, 805).EuclideanNormalize(), 4);
    IPrinter::printSeperatorLine();

    DoubleVector gn;
    gn.resize(x.length());
    IGradient::Gradient(&functional, 0.001, x, gn, static_cast<unsigned int>(0x0), static_cast<unsigned int>(400));
    //IGradient::Gradient(&functional, 0.01, x, gn, static_cast<unsigned int>(401), static_cast<unsigned int>(801));

    gn[000] = 0.0;
    gn[401] = 0.0;
    gn[802] = 0.0;
    gn[803] = 0.0;
    gn[804] = 0.0;
    gn[805] = 0.0;

    gn.EuclideanNormalize();
    IPrinter::printVector(gn.mid(0, 400).EuclideanNormalize());
    IPrinter::printVector(gn.mid(401, 801).EuclideanNormalize());
    IPrinter::print(gn.mid(802, 805).EuclideanNormalize(), 4);
    IPrinter::printSeperatorLine();


    //    Benchmark bencmark;
    //    bencmark.tick();
    //    DoubleMatrix u;
    //    functional.Problem0HForward::implicit_calculate_D2V1(u,functional.a, functional.gamma);
    //    bencmark.tock();
    //    std::cout << bencmark.CpuDurationSecond() << bencmark.CpuDurationClock() << std::endl;
    //    bencmark.tick();
    //    functional.Problem0HBckward::implicit_calculate_D2V1(u,functional.a, functional.gamma);
    //    bencmark.tock();
    //    std::cout << bencmark.CpuDurationSecond() << bencmark.CpuDurationClock() << std::endl;

    //    DoubleVector x; x.resize(2*200, 1.0);
    //    DoubleVector g; g.resize(2*200, 1.0);
    //    functional.gradient(x, g);
    //    IPrinter::printVector(g.mid(0, 199));
    //    IPrinter::printVector(g.mid(200, 1199));
}

Problem0HCommon::Problem0HCommon()
{
}

Problem0HCommon::~Problem0HCommon() {}

auto Problem0HFunctional::setDimension(const Dimension &timeDimension, const Dimension &dimensionX, const Dimension &dimensionY) -> void
{
    Problem0HForward::setTimeDimension(timeDimension);
    Problem0HForward::addSpaceDimension(dimensionX);
    Problem0HForward::addSpaceDimension(dimensionY);

    Problem0HBckward::setTimeDimension(timeDimension);
    Problem0HBckward::addSpaceDimension(dimensionX);
    Problem0HBckward::addSpaceDimension(dimensionY);

    //const double hx = spaceDimensionX.step();
    //const double hy = spaceDimensionY.step();
    const unsigned int N = static_cast<unsigned int> ( dimensionX.size() );
    const unsigned int M = static_cast<unsigned int> ( dimensionY.size() );
    const unsigned int L = static_cast<unsigned int> ( timeDimension.size() );

    U1.resize(M+1, N+1);
    U2.resize(M+1, N+1);
    u1.resize(M+1, N+1);
    u2.resize(M+1, N+1);

    ksi = SpacePoint(0.25, 0.25);

    const_this = const_cast<Problem0HFunctional*>(this);
    const_this->psi.resize(2);
    const_this->psi[0].v.resize(2*L+1);
    const_this->psi[1].v.resize(2*L+1);
}

auto Problem0HFunctional::fx(const DoubleVector &x) const -> double
{
    vectorToParameter(x);

    DoubleMatrix u;
    Problem0HForward::implicit_calculate_D2V1(u, Problem0HCommon::a, Problem0HCommon::gamma);
    double sum = 0.0;
    sum += integral1(u1);
    sum += alpha0*integral2(u2);
    //sum += norm();
    //sum += penalty();
    return sum;
}

auto Problem0HFunctional::integral1(const DoubleMatrix &u) const -> double
{
    const Dimension &dimensionX = Problem0HForward::spaceDimension(Dimension::DimensionX);
    const Dimension &dimensionY = Problem0HForward::spaceDimension(Dimension::DimensionY);
    const unsigned int N = static_cast<unsigned int> ( dimensionX.size() );
    const unsigned int M = static_cast<unsigned int> ( dimensionY.size() );
    const double hx = dimensionX.step();
    const double hy = dimensionY.step();

    double udiff = 0.0;
    double usum = 0.0;

    udiff = u[0][0]; usum += 0.25 * udiff * udiff * mu(0, 0);
    udiff = u[0][N]; usum += 0.25 * udiff * udiff * mu(N, 0);
    udiff = u[M][0]; usum += 0.25 * udiff * udiff * mu(0, M);
    udiff = u[M][N]; usum += 0.25 * udiff * udiff * mu(N, M);

    for (unsigned int n=1; n<=N-1; n++)
    {
        udiff = u[0][n]; usum += 0.5 * udiff * udiff * mu(n, 0);
        udiff = u[M][n]; usum += 0.5 * udiff * udiff * mu(n, M);
    }

    for (unsigned int m=1; m<=M-1; m++)
    {
        udiff = u[m][0]; usum += 0.5 * udiff * udiff * mu(0, m);
        udiff = u[m][N]; usum += 0.5 * udiff * udiff * mu(N, m);
    }

    for (unsigned int m=1; m<=M-1; m++)
    {
        for (unsigned int n=1; n<=N-1; n++)
        {
            udiff = u[m][n]; usum += udiff * udiff * mu(n, m);
        }
    }

    return usum*(hx*hy);
}

auto Problem0HFunctional::integral2(const DoubleMatrix &u) const -> double
{
    const Dimension &dimensionX = Problem0HForward::spaceDimension(Dimension::DimensionX);
    const Dimension &dimensionY = Problem0HForward::spaceDimension(Dimension::DimensionY);
    const unsigned int N = static_cast<unsigned int> ( dimensionX.size() );
    const unsigned int M = static_cast<unsigned int> ( dimensionY.size() );
    const double hx = dimensionX.step();
    const double hy = dimensionY.step();

    double udiff = 0.0;
    double usum = 0.0;

    udiff = u[0][0]; usum += 0.25 * udiff * udiff * mu(0, 0);
    udiff = u[0][N]; usum += 0.25 * udiff * udiff * mu(N, 0);
    udiff = u[M][0]; usum += 0.25 * udiff * udiff * mu(0, M);
    udiff = u[M][N]; usum += 0.25 * udiff * udiff * mu(N, M);

    for (unsigned int n=1; n<=N-1; n++)
    {
        udiff = u[0][n]; usum += 0.5 * udiff * udiff * mu(n, 0);
        udiff = u[M][n]; usum += 0.5 * udiff * udiff * mu(n, M);
    }

    for (unsigned int m=1; m<=M-1; m++)
    {
        udiff = u[m][0]; usum += 0.5 * udiff * udiff * mu(0, m);
        udiff = u[m][N]; usum += 0.5 * udiff * udiff * mu(N, m);
    }

    for (unsigned int m=1; m<=M-1; m++)
    {
        for (unsigned int n=1; n<=N-1; n++)
        {
            udiff = u[m][n]; usum += udiff * udiff * mu(n, m);
        }
    }

    return usum*(hx*hy);
}

auto Problem0HFunctional::norm() const -> double { return 0.0; }

auto Problem0HFunctional::penalty() const -> double { return 0.0; }

auto Problem0HFunctional::gradient(const DoubleVector &x, DoubleVector &g) const -> void
{
    vectorToParameter(x);

    g.clear(); g.resize(x.length());

    const Dimension &timeDimension = Problem0HForward::timeDimension();
    const unsigned int L = static_cast<unsigned int> ( timeDimension.size() );

    DoubleMatrix u, p;
    Problem0HForward::implicit_calculate_D2V1(u, a, gamma);
    Problem0HBckward::implicit_calculate_D2V1(p, a, gamma);

    unsigned int size = 2*L;
    for (unsigned int ln=0; ln<=size; ln++)
    {
        g[ln+0*L]   = -psi[0].psi_vl[ln];
        g[ln+2*L+1] = -psi[1].psi_vl[ln];
    }
    g[2*(2*L+1)+0] = 0.0;
    g[2*(2*L+1)+1] = 0.0;
    g[2*(2*L+1)+2] = 0.0;
    g[2*(2*L+1)+3] = 0.0;

    g[000] = 0.0;
    g[401] = 0.0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////

auto Problem0HForward::layerInfo(const DoubleMatrix &u, unsigned int ln) const -> void
{
    const unsigned int L = static_cast<unsigned int>(timeDimension().size());
    const double ht = timeDimension().step();

    Problem0HForward *forward = const_cast<Problem0HForward*>(this);
    if (ln == 2*(L-2)) forward->uL2 = u;
    if (ln == 2*(L-1)) forward->uL1 = u;
    if (ln == 2*(L-0)) forward->uL0 = u;
    if (ln == 2*(L-0))
    {
        forward->u1 = u;
        forward->u2 = (1.0/(2.0*ht))*(uL2 - 4.0*uL1 + 3.0*uL0);
    }

    //QString filename1 = QString("data/problem0H/f/txt/f_%1.txt").arg(ln, 4, 10, QChar('0'));
    //IPrinter::print(u,filename1.toLatin1().data());
    //IPrinter::printSeperatorLine();
    //printf("Forward: %4d %0.3f %10.8f %10.8f %4d %4d\n", ln, ln*0.005, u.min(), u.max(), 0, 0);
    //QString filename2 = QString("data/problem0H/f/png/f_%1.png").arg(ln, 4, 10, QChar('0'));
    //QPixmap pixmap;
    //visualizeMatrixHeat(u, u.min(), u.max(), pixmap, 101, 101);
    //pixmap.save(filename2);
    //IPrinter::printSeperatorLine();
}

auto Problem0HForward::initial(const SpaceNodePDE &, InitialCondition) const -> double { return 0.0; }

auto Problem0HForward::boundary(const SpaceNodePDE &, const TimeNodePDE &) const -> double { return 0.0; }

auto Problem0HForward::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double
{
    double pv = p(sn, tn);

    //static const double sigma = 0.01;
    //static const double cff1 = 1.0/(2.0*M_PI*sigma*sigma);
    //static const double cff2 = 1.0/(2.0*sigma*sigma);

    unsigned int ln = static_cast<unsigned int>(tn.i);
    double _v1 = psi[0].v[ln];
    double _v2 = psi[1].v[ln];


    double pulse1 = _v1 * psi[0].deltaGrid.weight(sn);//cff1 * exp(-cff2*((sn.x-psi[0].p.x)*(sn.x-psi[0].p.x)+(sn.y-psi[0].p.y)*(sn.y-psi[0].p.y)));
    double pulse2 = _v2 * psi[1].deltaGrid.weight(sn);//cff1 * exp(-cff2*((sn.x-psi[1].p.x)*(sn.x-psi[1].p.x)+(sn.y-psi[1].p.y)*(sn.y-psi[1].p.y)));

    //double pulse1, pulse2; pulse1 = pulse2 = 0.0;

    return pv + pulse1 + pulse2;
}

auto Problem0HForward::p(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double
{
    static const double sigma = 0.005*4.0;
    static const double alpha1 = 1.0/(2.0*M_PI*sigma*sigma);
    static const double alpha2 = 1.0/(2.0*sigma*sigma);
    static const double alpha3 = 100.0;
    double pv = 5.0 * alpha1 * exp( -alpha2 * ((sn.x-ksi.x)*(sn.x-ksi.x)+(sn.y-ksi.y)*(sn.y-ksi.y)) - alpha3*tn.t);
    return pv;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////

auto Problem0HBckward::layerInfo(const DoubleMatrix &p, unsigned int ln) const -> void
{
    Problem0HBckward* const_this = const_cast<Problem0HBckward*>(this);
    const_this->psi[0].psi_vl[ln] = psi[0].deltaGrid.consentrateInPoint(p);
    const_this->psi[1].psi_vl[ln] = psi[1].deltaGrid.consentrateInPoint(p);
    //printf("%d %f %f %d %d\n", ln, const_this->psi[0].psi_vl[ln], const_this->psi[1].psi_vl[ln], psi[0].deltaGrid.minX(), psi[0].deltaGrid.maxX());

    //IPrinter::printMatrix(p);
    //IPrinter::printSeperatorLine();

    //QString filename1 = QString("data/problem0H/b/txt/b_%1.txt").arg(ln, 4, 10, QChar('0'));
    //IPrinter::print(p,filename1.toLatin1().data());
    //IPrinter::printSeperatorLine();
    //printf("Forward: %4d %0.3f %10.8f %10.8f %4d %4d\n", ln, ln*0.005, u.min(), u.max(), 0, 0);
    QString filename2 = QString("data/problem0H/b/png/b_%1.png").arg(ln, 4, 10, QChar('0'));
    QPixmap pixmap;
    visualizeMatrixHeat(p, p.min(), p.max(), pixmap, 201, 201);
    pixmap.save(filename2);
    //IPrinter::printSeperatorLine();
}

auto Problem0HBckward::initial(const SpaceNodePDE &sn, InitialCondition condition) const -> double
{
    unsigned int n = static_cast<unsigned int>(sn.i);
    unsigned int m = static_cast<unsigned int>(sn.j);
    if (condition == InitialCondition::InitialValue)
    {
        return -2.0*alpha0*(u2[m][n]-U2[m][n]);
    }
    else
    {
        return 2.0*(u1[m][n]-U1[m][n]) + gamma*initial(sn, InitialCondition::InitialValue);
    }
}

auto Problem0HBckward::boundary(const SpaceNodePDE&, const TimeNodePDE&) const -> double { return 0.0; }

auto Problem0HBckward::f(const SpaceNodePDE&, const TimeNodePDE&) const -> double { return 0.0; }

///////////////////////////////////////////////////////////////////////////////////////////////////////////////


auto Problem0HFunctional::vectorToParameter(const DoubleVector &x) const -> void
{    
    const Dimension &time = Problem0HForward::timeDimension();
    const Dimension &dimX = Problem0HForward::spaceDimension(Dimension::DimensionX);
    const Dimension &dimY = Problem0HForward::spaceDimension(Dimension::DimensionY);

    const unsigned int N = static_cast<unsigned int> ( dimX.size() );
    const unsigned int M = static_cast<unsigned int> ( dimY.size() );
    const unsigned int L = static_cast<unsigned int> ( time.size() );
    const double hy = dimX.step();
    const double hx = dimY.step();
    //const double ht = time.size();

    Problem0HFunctional* const_this = const_cast<Problem0HFunctional*>(this);

    const_this->psi[0].v.resize(2*L+1);
    const_this->psi[1].v.resize(2*L+1);
    const_this->psi[0].psi_vl.resize(2*L+1);
    const_this->psi[1].psi_vl.resize(2*L+1);
    const_this->psi[0].psi_dx.resize(2*L+1);
    const_this->psi[1].psi_dx.resize(2*L+1);
    const_this->psi[0].psi_dy.resize(2*L+1);
    const_this->psi[1].psi_dy.resize(2*L+1);

    unsigned int offset = 2*L;
    for (unsigned int i=0; i<=2*L; i++)
    {
        const_this->psi[0].v[i] = x[i];
        const_this->psi[1].v[i] = x[i+offset];
    }
    offset = 4*L+2;
    const_this->psi[0].p.x = x[offset+0];
    const_this->psi[0].p.y = x[offset+1];
    const_this->psi[1].p.x = x[offset+2];
    const_this->psi[1].p.y = x[offset+3];

    const_this->psi[0].deltaGrid.cleanGrid();
    const_this->psi[0].deltaGrid.initGrid(N, hx, M, hy);
    const_this->psi[0].deltaGrid.distributeGauss(psi[0].p);

    const_this->psi[1].deltaGrid.cleanGrid();
    const_this->psi[1].deltaGrid.initGrid(N, hx, M, hy);
    const_this->psi[1].deltaGrid.distributeGauss(psi[1].p);
}

auto Problem0HFunctional::parameterToVector(DoubleVector &x) const -> void
{
    const unsigned int L = static_cast<unsigned int> ( Problem0HForward::timeDimension().size() );
    Problem0HFunctional* const_this = const_cast<Problem0HFunctional*>(this);
    x.clear();
    x.resize(806);

    unsigned int offset = 2*L;
    for (unsigned int i=0; i<=2*L; i++)
    {
        x[i] = const_this->psi[0].v[i];
        x[i+offset] = const_this->psi[1].v[i];
    }
    offset = 4*L+2;
    x[offset+0] = const_this->psi[0].p.x;
    x[offset+1] = const_this->psi[0].p.y;
    x[offset+2] = const_this->psi[1].p.x;
    x[offset+3] = const_this->psi[1].p.y;
}
