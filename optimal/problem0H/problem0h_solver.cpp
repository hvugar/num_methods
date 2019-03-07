#include "problem0h_solver.h"

void Problem0HFunctional::Main(int argc, char **argv)
{
    QGuiApplication app(argc, argv);


    unsigned int N = 100; double hx = 0.01;
    unsigned int M = 100; double hy = 0.01;
    unsigned int L = 100; double ht = 0.01;

    Problem0HFunctional functional;
    functional.ksi = SpacePoint(0.25, 0.25);
    functional.a = 1.0;
    functional.gamma = 0.0;
    functional.alpha0 = 1.0;
    functional.source_number = 2;
    functional.setDimension(Dimension(ht, 0, static_cast<int>(L)),
                            Dimension(hx, 0, static_cast<int>(N)),
                            Dimension(hy, 0, static_cast<int>(M)));
    functional.optimalParameters[0].p = SpacePoint(0.35, 0.65);
    functional.optimalParameters[1].p = SpacePoint(0.75, 0.25);

    DoubleVector x;
    functional.parameterToVector(x);

    unsigned int start1 = 0;
    unsigned int finish1 = 2*L;
    unsigned int start2 = 2*L+1;
    unsigned int finish2 = 2*L+1+2*L;
    unsigned int start3 = 2*L+1+2*L+1;
    unsigned int finish3 = 2*L+1+2*L+4;

    DoubleVector ga;
    functional.gradient(x, ga);
    IPrinter::printVector(ga.mid(start1, finish1).EuclideanNormalize());
    IPrinter::printVector(ga.mid(start2, finish2).EuclideanNormalize());
    IPrinter::print(ga.mid(start3, finish3).EuclideanNormalize(), 4);
    IPrinter::printSeperatorLine();

    DoubleVector gn;
    gn.resize(x.length());
    //IGradient::Gradient(&functional, 0.05, x, gn, static_cast<unsigned int>(start1), static_cast<unsigned int>(finish1));
    //IGradient::Gradient(&functional, 0.05, x, gn, static_cast<unsigned int>(start2), static_cast<unsigned int>(finish2));
    IGradient::Gradient(&functional, 0.01, x, gn, static_cast<unsigned int>(start3), static_cast<unsigned int>(finish3));
    //IGradient::Gradient(&functional, 0.05, x, gn);

    unsigned int length = 2*L+1;
    for (unsigned int sn=0; sn<functional.source_number; sn++)
    {
        gn[sn*length+0] = 0.0;
        gn[sn*length+1] = 0.0;
        gn[sn*length+length-1] = 0.0;

        gn[functional.source_number*length+2*sn+0] = 0.0;
        gn[functional.source_number*length+2*sn+1] = 0.0;
    }

    IPrinter::printVector(gn.mid(start1, finish1).EuclideanNormalize());
    IPrinter::printVector(gn.mid(start2, finish2).EuclideanNormalize());
    IPrinter::print(gn.mid(start3, finish3).EuclideanNormalize(), 4);
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

Problem0HCommon::Problem0HCommon() {}

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

    U1.resize(M+1, N+1, 0.0);
    U2.resize(M+1, N+1, 0.0);
    u1.resize(M+1, N+1, 0.0);
    u2.resize(M+1, N+1, 0.0);

    ksi = SpacePoint(0.25, 0.25);

    const_this = const_cast<Problem0HFunctional*>(this);
    const_this->optimalParameters.resize(source_number);

    unsigned int length = 2*L+1;
    for (unsigned int sn=0; sn<source_number; sn++)
    {
        const_this->optimalParameters[sn].pwr_vl.resize(length, 0.10);
        const_this->optimalParameters[sn].pwr_vl[0] = 0.0;
        const_this->optimalParameters[sn].pwr_vl[1] = 0.0;
        const_this->optimalParameters[sn].pwr_vl[2*L] = 0.0;
        const_this->optimalParameters[sn].psi_vl.resize(length);
        const_this->optimalParameters[sn].psi_dx.resize(length);
        const_this->optimalParameters[sn].psi_dy.resize(length);
    }
}

auto Problem0HFunctional::fx(const DoubleVector &x) const -> double
{
    vectorToParameter(x);

    DoubleMatrix u;
    Problem0HForward::implicit_calculate_D2V1(u, a, gamma);
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
    const unsigned int N = static_cast<unsigned int>(dimensionX.size());
    const unsigned int M = static_cast<unsigned int>(dimensionY.size());
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
    const double ht = static_cast<double>(timeDimension.size());

    DoubleMatrix u, p;
    Problem0HForward::implicit_calculate_D2V1(u, a, gamma);
    Problem0HBckward::implicit_calculate_D2V1(p, a, gamma);

    unsigned int length = 2*L+1;
    for (unsigned int sn=0; sn<source_number; sn++)
    {
        for (unsigned int ln=0; ln<length; ln++)
        {
            g[sn*length+ln] = -optimalParameters[sn].psi_vl[ln];
        }
        g[sn*length+0] = 0.0;
        g[sn*length+1] = 0.0;
        g[sn*length+length-1] = 0.0;

        //////////////////////////////////////////////////////////////////

        unsigned ix = source_number*length+2*sn+0;
        unsigned iy = source_number*length+2*sn+1;
        g[ix] = 0.0;
        g[iy] = 0.0;
        g[ix] += 0.5*optimalParameters[sn].psi_dx[2];
        g[iy] += 0.5*optimalParameters[sn].psi_dy[2];
        for (unsigned int ln=3; ln<=2*(L-1); ln++)
        {
            g[ix] += optimalParameters[sn].psi_dx[ln];
            g[iy] += optimalParameters[sn].psi_dy[ln];
        }
        g[ix] += 0.5*optimalParameters[sn].psi_dx[2*(L-1)+1];
        g[iy] += 0.5*optimalParameters[sn].psi_dy[2*(L-1)+1];

        g[ix] *= 0.5*ht;
        g[iy] *= 0.5*ht;
    }
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
        //        IPrinter::printSeperatorLine((std::string("uL2") + std::to_string(ln)).data()); IPrinter::printMatrix(uL2);
        //        IPrinter::printSeperatorLine((std::string("uL1") + std::to_string(ln)).data()); IPrinter::printMatrix(uL1);
        //        IPrinter::printSeperatorLine((std::string("uL0") + std::to_string(ln)).data()); IPrinter::printMatrix(uL0);

        //lxw_workbook  *workbook  = new_workbook("d:\\chart.xlsx");
        //lxw_worksheet *worksheet1 = workbook_add_worksheet(workbook, nullptr);
        //lxw_chartsheet *chartsheet = workbook_add_chartsheet(workbook, nullptr);
        //lxw_chart *chart = workbook_add_chart(workbook, LXW_CHART_AREA);

        //lxw_format *format = workbook_add_format(workbook);
        //format_set_num_format(format, "0.0000000000");

        //for (unsigned int row=0; row<u.rows(); row++)
        //{
        //    for (unsigned int col=0; col<u.cols(); col++)
        //    {
        //        // =Sheet1!$A$1:$A$6
        //        // =Sheet1!$A$1:$CW$1
        //        std::string line = "=Sheet1!$A$"+std::to_string(row+1)+":$CW$"+std::to_string(row+1);
        //        worksheet_write_number(worksheet1, static_cast<lxw_row_t>(row), static_cast<lxw_col_t>(col), u[row][col], format);
        //        //chart_add_series(chart, nullptr, line.data());
        //    }
        //}
        //chartsheet_set_chart(chartsheet, chart);
        //workbook_close(workbook);

        forward->u1 = u;
        forward->u2 = (1.0/(2.0*ht))*(uL2 - 4.0*uL1 + 3.0*uL0);

        //IPrinter::printMatrix(forward->u2);
    }

    //QString filename1 = QString("data/problem0H/f/txt/f_%1.txt").arg(ln, 4, 10, QChar('0'));
    //IPrinter::print(u,filename1.toLatin1().data());
    //IPrinter::printSeperatorLine();
    //printf("Forward: %4d %0.3f %10.8f %10.8f %4d %4d\n", ln, ln*0.005, u.min(), u.max(), 0, 0);
    //if (ln == 2*(L-2) || ln == 2*(L-1) || ln == 2*(L-0))
    //{
    //    QString filename2 = QString("data/problem0H/f/png/f_%1.png").arg(ln, 4, 10, QChar('0'));
    //    QPixmap pixmap;
    //    visualizeMatrixHeat(u, u.min(), u.max(), pixmap, 101, 101);
    //    pixmap.save(filename2);
    //}
    //IPrinter::printSeperatorLine();
}

auto Problem0HForward::initial(const SpaceNodePDE &, InitialCondition) const -> double { return 0.0; }

auto Problem0HForward::boundary(const SpaceNodePDE &, const TimeNodePDE &) const -> double { return 0.0; }

auto Problem0HForward::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double
{
    double pv = p(sn, tn);
    double pulse1, pulse2; pulse1 = pulse2 = 0.0;
    if (tn.i != 0 && tn.i != 1)
    {
        //static const double sigma = 0.01;
        //static const double cff1 = 1.0/(2.0*M_PI*sigma*sigma);
        //static const double cff2 = 1.0/(2.0*sigma*sigma);

        unsigned int ln = static_cast<unsigned int>(tn.i);
        double _v1 = optimalParameters[0].pwr_vl[ln];
        double _v2 = optimalParameters[1].pwr_vl[ln];

        pulse1 = _v1 * optimalParameters[0].deltaGrid.weight(sn);//cff1 * exp(-cff2*((sn.x-psi[0].p.x)*(sn.x-psi[0].p.x)+(sn.y-psi[0].p.y)*(sn.y-psi[0].p.y)));
        pulse2 = _v2 * optimalParameters[1].deltaGrid.weight(sn);//cff1 * exp(-cff2*((sn.x-psi[1].p.x)*(sn.x-psi[1].p.x)+(sn.y-psi[1].p.y)*(sn.y-psi[1].p.y)));
    }
    return pv + pulse1 + pulse2;
}

auto Problem0HForward::p(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double
{
    static const double sigmaX = 8.0*Problem0HForward::spaceDimension(Dimension::DimensionX).step();
    static const double sigmaY = 8.0*Problem0HForward::spaceDimension(Dimension::DimensionY).step();;
    static const double alpha1 = 1.0/(2.0*M_PI*sigmaX*sigmaY);
    static const double alpha2 = 1.0/(2.0*sigmaX*sigmaY);
    static const double alpha3 = 100.0;
    double pv = alpha1 * exp( -alpha2 * ((sn.x-ksi.x)*(sn.x-ksi.x)+(sn.y-ksi.y)*(sn.y-ksi.y)) - alpha3*tn.t);
    return pv;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////

auto Problem0HBckward::layerInfo(const DoubleMatrix &p, unsigned int ln) const -> void
{
    Problem0HBckward* const_this = const_cast<Problem0HBckward*>(this);

    for (unsigned int sn=0; sn<source_number; sn++)
    {
        double psi_dx, psi_dy;
        double psi_vl = optimalParameters[sn].deltaGrid.consentrateInPoint(p, psi_dx, psi_dy);
        const_this->optimalParameters[sn].psi_vl[ln] = psi_vl;
        const_this->optimalParameters[sn].psi_dx[ln] = psi_dx;
        const_this->optimalParameters[sn].psi_dy[ln] = psi_dy;

        //IPrinter::printSeperatorLine(std::to_string(ln).data());
        //IPrinter::printMatrix(p);
    }

    //IPrinter::printMatrix(p);
    //IPrinter::printSeperatorLine();

    //QString filename1 = QString("data/problem0H/b/txt/b_%1.txt").arg(ln, 4, 10, QChar('0'));
    //IPrinter::print(p,filename1.toLatin1().data());
    //IPrinter::printSeperatorLine();
    //printf("Forward: %4d %0.3f %10.8f %10.8f %4d %4d\n", ln, ln*0.005, u.min(), u.max(), 0, 0);
    //QString filename2 = QString("data/problem0H/b/png/b_%1.png").arg(ln, 4, 10, QChar('0'));
    //QPixmap pixmap;
    //visualizeMatrixHeat(p, p.min(), p.max(), pixmap, 201, 201);
    //pixmap.save(filename2);
    //IPrinter::printSeperatorLine();
}

auto Problem0HBckward::initial(const SpaceNodePDE &sn, InitialCondition condition) const -> double
{
    unsigned int n = static_cast<unsigned int>(sn.i);
    unsigned int m = static_cast<unsigned int>(sn.j);
    if (condition == InitialCondition::InitialValue)
    {
        return -2.0*alpha0*(u2[m][n]/*-U2[m][n]*/);
    }
    else
    {
        return +2.0*(u1[m][n]/*-U1[m][n]*/) + gamma*initial(sn, InitialCondition::InitialValue);
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
    const_this->optimalParameters.resize(source_number);

    const unsigned int length = 2*L+1;
    for (unsigned int sn=0; sn<source_number; sn++)
    {
        const_this->optimalParameters[sn].pwr_vl.resize(length);
        const_this->optimalParameters[sn].psi_vl.resize(length);
        const_this->optimalParameters[sn].psi_dx.resize(length);
        const_this->optimalParameters[sn].psi_dy.resize(length);
        for (unsigned int i=0; i<length; i++)
        {
            const_this->optimalParameters[sn].pwr_vl[i] = x[sn*length+i];
        }
        const_this->optimalParameters[sn].pwr_vl[0] = 0.0;
        const_this->optimalParameters[sn].pwr_vl[1] = 0.0;
        const_this->optimalParameters[sn].pwr_vl[length-1] = 0.0;

        const_this->optimalParameters[sn].p.x = x[source_number*length+2*sn+0];
        const_this->optimalParameters[sn].p.y = x[source_number*length+2*sn+1];

        const_this->optimalParameters[sn].deltaGrid.cleanGrid();
        const_this->optimalParameters[sn].deltaGrid.initGrid(N, hx, M, hy);
        const_this->optimalParameters[sn].deltaGrid.distributeGauss(const_this->optimalParameters[sn].p);
    }
}

auto Problem0HFunctional::parameterToVector(DoubleVector &x) const -> void
{
    const unsigned int L = static_cast<unsigned int> ( Problem0HForward::timeDimension().size() );
    Problem0HFunctional* const_this = const_cast<Problem0HFunctional*>(this);

    unsigned int size = (2*L+3)*source_number;
    x.resize(size);

    unsigned int length = 2*L+1;
    for (unsigned int sn=0; sn<source_number; sn++)
    {
        for (unsigned int i=0; i<length; i++)
        {
            x[sn*length+i] = const_this->optimalParameters[sn].pwr_vl[i];
        }
        x[sn*length+0] = 0.0;
        x[sn*length+1] = 0.0;
        x[sn*length+length-1] = 0.0;

        x[source_number*length+2*sn+0] = const_this->optimalParameters[sn].p.x;
        x[source_number*length+2*sn+1] = const_this->optimalParameters[sn].p.y;
    }
}
