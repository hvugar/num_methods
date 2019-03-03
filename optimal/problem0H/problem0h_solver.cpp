﻿#include "problem0h_solver.h"

void Problem0HFunctional::Main(int argc, char **argv)
{
    QGuiApplication app(argc, argv);

    Problem0HFunctional functional;
    functional.ksi = SpacePoint(0.25, 0.25);
    functional.a = 1.0;
    functional.gamma = 0.0;
    functional.alpha0 = 1.0;
    functional.setDimension(Dimension(0.005, 0, 600), Dimension(0.005, 0, 200), Dimension(0.005, 0, 200));

    Benchmark bencmark;
    bencmark.tick();
    DoubleMatrix u;
    functional.Problem0HForward::implicit_calculate_D2V1(u,functional.a, functional.gamma);
    bencmark.tock();
    std::cout << bencmark.CpuDurationSecond() << bencmark.CpuDurationClock() << std::endl;
    bencmark.tick();
    functional.Problem0HBckward::implicit_calculate_D2V1(u,functional.a, functional.gamma);
    bencmark.tock();
    std::cout << bencmark.CpuDurationSecond() << bencmark.CpuDurationClock() << std::endl;

    DoubleVector x; x.resize(2*600, 1.0);
    DoubleVector g; g.resize(2*600, 1.0);
    functional.gradient(x, g);
    IPrinter::printVector(g.mid(0, 599));
    IPrinter::printVector(g.mid(600, 1199));
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

    U1.resize(M+1, N+1);
    U2.resize(M+1, N+1);
    u1.resize(M+1, N+1);
    u2.resize(M+1, N+1);

    ksi = SpacePoint(0.25, 0.25);
}

auto Problem0HFunctional::fx(const DoubleVector &x) const -> double
{
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

    const double hx = dimensionX.step();
    const double hy = dimensionY.step();
    const unsigned int N = static_cast<unsigned int> ( dimensionX.size() );
    const unsigned int M = static_cast<unsigned int> ( dimensionY.size() );

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

    const double hx = dimensionX.step();
    const double hy = dimensionY.step();
    const unsigned int N = static_cast<unsigned int> ( dimensionX.size() );
    const unsigned int M = static_cast<unsigned int> ( dimensionY.size() );

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
    g.clear(); g.resize(x.length());

    const Dimension &timeDimension = Problem0HForward::timeDimension();
    const double ht = timeDimension.step();
    const unsigned int L = static_cast<unsigned int> ( timeDimension.size() );

    DoubleMatrix u, p;
    Problem0HForward::implicit_calculate_D2V1(u, a, gamma);
    Problem0HBckward::implicit_calculate_D2V1(p, a, gamma);

    for (unsigned int i=0; i<L; i++)
    {
        g[i+0*L] = -Problem0HCommon::ps1[i].vl;
        g[i+1*L] = -Problem0HCommon::ps2[i].vl;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////

auto Problem0HForward::layerInfo(const DoubleMatrix &u, unsigned int ln) const -> void
{
    const unsigned int M = static_cast<unsigned int>(timeDimension().max());
    const double ht = timeDimension().step();

    Problem0HForward *forward = const_cast<Problem0HForward*>(this);
    if (ln == M-2) forward->u2 = u;
    if (ln == M-1) forward->u1 = u;
    if (ln == M)
    {
        forward->u1 = u;
        forward->u2 = u2 - 4.0*u1 + 3.0*u;
        forward->u2 *= (1.0/(2.0*ht));
    }


//    std::string filename = std::string("d:/data/files/f/txt/") + std::to_string(ln) + std::string(".txt");
//    IPrinter::print(u,filename.data());
//    IPrinter::printSeperatorLine();
//    printf("Forward: %4d %0.3f %10.8f %10.8f %4d %4d\n", ln, ln*0.005, u.min(), u.max(), 0, 0);

//    std::string filename1 = std::string("d:/data/files/f/png/")+std::to_string(ln) + std::string(".png");
//    QPixmap pixmap;
//    visualizeMatrixHeat(u, u.min(), u.max(), pixmap, 201, 201);
//    pixmap.save(QString(filename1.data()));

//    IPrinter::printSeperatorLine();
}

auto Problem0HForward::initial1(const SpaceNodePDE &) const -> double { return 0.0; }

auto Problem0HForward::initial2(const SpaceNodePDE &) const -> double { return 0.0; }

auto Problem0HForward::boundary(const SpaceNodePDE &, const TimeNodePDE &) const -> double { return 0.0; }

auto Problem0HForward::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double
{
    double pv = p(sn, tn);

    static const double sigma = 0.005;
    static const double alpha1 = 1.0/(2.0*M_PI*sigma*sigma);
    static const double alpha2 = 1.0/(2.0*sigma*sigma);
    static const double alpha3 = 100.0;

    //    double pulse1 = 1.0 * alpha1 * exp(-alpha2*((sn.x-common->p1.x)*(sn.x-common->p1.x)
    //                                                          +(sn.y-common->p1.y)*(sn.y-common->p1.y)));
    //    double pulse2 = 1.0 * alpha1 * exp(-alpha2*((sn.x-common->p2.x)*(sn.x-common->p2.x)
    //                                                          +(sn.y-common->p2.y)*(sn.y-common->p2.y)));

    return pv;// + pulse1 + pulse2;
}

auto Problem0HForward::p(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double
{
    static const double sigma = 0.005*4.0;
    static const double alpha1 = 1.0/(2.0*M_PI*sigma*sigma);
    static const double alpha2 = 1.0/(2.0*sigma*sigma);
    static const double alpha3 = 100.0;

    double pv = 5.0 * alpha1 * exp( -alpha2 * ((sn.x-ksi.x)*(sn.x-ksi.x)+(sn.y-ksi.y)*(sn.y-ksi.y)) - alpha3*tn.t);
    //    if (pv > 0.0000000001)
    //    {
    //        count1++;
    //        //printf("%4d %5.4f %4d %4.3f %4d %4.3f %20.14f %4d\n", tn.i, tn.t, sn.i, sn.x, sn.j, sn.y, pv, count1);
    //    }
    //    else
    //    {
    //        pv = 0.0;
    //        count2++;
    //    }

    return pv;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////

auto Problem0HBckward::layerInfo(const DoubleMatrix &p, unsigned int ln) const -> void
{   
//    std::string filename = std::string("d:/data/files/b/txt/")+std::to_string(ln) + std::string(".txt");
//    IPrinter::print(p,filename.data());
//    //IPrinter::printSeperatorLine();
//    printf("Backward: %4d %0.3f %10.8f %10.8f %4d %4d\n", ln, ln*0.005, p.min(), p.max(), 0, 0);

//    std::string filename1 = std::string("d:/data/files/b/png/")+std::to_string(ln) + std::string(".png");
//    QPixmap pixmap;
//    visualizeMatrixHeat(p, p.min(), p.max(), pixmap, 201, 201);
//    pixmap.save(QString(filename1.data()));

    //IPrinter::printSeperatorLine();
}

auto Problem0HBckward::initial1(const SpaceNodePDE &sn) const -> double
{
    unsigned int n = static_cast<unsigned int>(sn.i);
    unsigned int m = static_cast<unsigned int>(sn.j);
    return -2.0*alpha0*(u2[m][n]-U2[m][n]);
}

auto Problem0HBckward::initial2(const SpaceNodePDE &sn) const -> double
{
    unsigned int n = static_cast<unsigned int>(sn.i);
    unsigned int m = static_cast<unsigned int>(sn.j);
    return 2.0*(u1[m][n]-U1[m][n]) + gamma*initial1(sn);
}

auto Problem0HBckward::boundary(const SpaceNodePDE&, const TimeNodePDE&) const -> double { return 0.0; }

auto Problem0HBckward::f(const SpaceNodePDE&, const TimeNodePDE&) const -> double { return 0.0; }

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
