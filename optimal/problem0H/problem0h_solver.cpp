#include "problem0h_solver.h"

void Problem0HFunctional::Main(int argc, char **argv)
{
#ifdef USE_LIB_IMAGING
    QGuiApplication app(argc, argv);
#endif

    //checkingForwardProblem();
    compareGradients();
}

void Problem0HFunctional::compareGradients()
{
    Problem0HFunctional functional;
    functional.ksi = SpacePoint(0.25, 0.25);
    functional.p_sigmaX = 0.05;
    functional.p_sigmaY = 0.05;
    functional.p_sigmaT = 0.01;
    functional.a = 1.0;
    functional.gamma = 0.0;
    functional.epsilon1 = 1.0;
    functional.epsilon2 = 1.0;
    functional.source_number = 2;
    functional.setDimension(Dimension(0.005, 0, 200), Dimension(0.01, 0, 100), Dimension(0.01, 0, 100));
    functional.optimalParameters[0].distribute(SpacePoint(0.300, 0.600));
    functional.optimalParameters[1].distribute(SpacePoint(0.700, 0.200));

    DoubleVector x;
    functional.parameterToVector(x);

    unsigned int L = 200;
    unsigned int s1 = 0*(2*L+1);   unsigned int f1 = 0*(2*L+1)+2*L;
    unsigned int s2 = 1*(2*L+1);   unsigned int f2 = 1*(2*L+1)+2*L;
    unsigned int s3 = 2*(2*L+1)+0; unsigned int f3 = 2*(2*L+1)+3;
    //unsigned int start31 = 2*(2*L+1)+0; unsigned int finish31 = 2*(2*L+1)+1;
    //unsigned int start32 = 2*(2*L+1)+2; unsigned int finish32 = 2*(2*L+1)+3;

    DoubleVector ga;
    functional.gradient(x, ga);

    //IPrinter::printVector(ga.mid(start1, finish1).EuclideanNormalize());
    //IPrinter::printVector(ga.mid(start2, finish2).EuclideanNormalize());
    //IPrinter::print(ga.mid(0, 4).EuclideanNormalize(), 5);
    IPrinter::print(ga.mid(s3, f3).EuclideanNormalize(), 4);
    //IPrinter::print(ga.mid(start31, finish31).EuclideanNormalize(), 2);
    //IPrinter::print(ga.mid(start32, finish32).EuclideanNormalize(), 2);
    //IPrinter::printSeperatorLine();

    DoubleVector gn;
    gn.resize(x.length());

    //IGradient::Gradient(&functional, 0.05, x, gn, static_cast<unsigned int>(s1), static_cast<unsigned int>(f1));
    //IGradient::Gradient(&functional, 0.05, x, gn, static_cast<unsigned int>(s2), static_cast<unsigned int>(f2));
    IGradient::Gradient(&functional, 0.01, x, gn, static_cast<unsigned int>(s3), static_cast<unsigned int>(f3));

    //IGradient::Gradient(&functional, 0.01, x, gn, static_cast<unsigned int>(start31), static_cast<unsigned int>(finish31));
    //IGradient::Gradient(&functional, 0.01, x, gn, static_cast<unsigned int>(start32), static_cast<unsigned int>(finish32));
    IGradient::Gradient(&functional, 0.05, x, gn, static_cast<unsigned int>(s3), static_cast<unsigned int>(s3));

    //unsigned int length = 2*L+1;
    //for (unsigned int sn=0; sn<functional.source_number; sn++)
    //{
    //    gn[sn*length+0] = 0.0;
    //    gn[sn*length+1] = 0.0;
    //    gn[sn*length+length-1] = 0.0;
    //}

    //    IPrinter::printVector(gn.mid(start1, finish1).EuclideanNormalize());
    //    IPrinter::printVector(gn.mid(start2, finish2).EuclideanNormalize());
    //IPrinter::print(gn.mid(0, 4).EuclideanNormalize(), 5);
    IPrinter::print(gn.mid(s3, f3).EuclideanNormalize(), 4);
    //IPrinter::print(gn.mid(start31, finish31).EuclideanNormalize(), 2);
    //IPrinter::print(gn.mid(start32, finish32).EuclideanNormalize(), 2);
    IPrinter::printSeperatorLine();
}

void Problem0HFunctional::checkingForwardProblem()
{
    Problem0HFunctional fw1;
    fw1.source_number = 2;
    fw1.ksi = SpacePoint(0.25, 0.36);
    fw1.p_sigmaX = 0.05;
    fw1.p_sigmaY = 0.05;
    fw1.p_sigmaT = 0.01;

    fw1.setDimension(Dimension(0.005, 0, 2000), Dimension(0.01, 0, 100), Dimension(0.01, 0, 100));
    fw1.optimalParameters[0].distribute(SpacePoint(0.308, 0.608));
    fw1.optimalParameters[1].distribute(SpacePoint(0.708, 0.208));

    DoubleMatrix u1;puts("3");
    fw1.forward(u1, 1.0, 0.0);puts("3");
    IPrinter::printMatrix(u1);
    IPrinter::printSeperatorLine();
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

auto Problem0HForward::layerInfo(const DoubleMatrix &u, const TimeNodePDE &tn) const -> void
{
    //return;

    //if (abs(tn.t - 0.1) <= 0.0) { IPrinter::printMatrix(10,6,u); IPrinter::printSeperatorLine(); }
    //if (abs(tn.t - 0.2) <= 0.0) { IPrinter::printMatrix(10,6,u); IPrinter::printSeperatorLine(); }
    //if (abs(tn.t - 0.3) <= 0.0) { IPrinter::printMatrix(10,6,u); IPrinter::printSeperatorLine(); }
    //if (abs(tn.t - 1.0) <= 0.0) { IPrinter::printMatrix(10,6,u); IPrinter::printSeperatorLine(); }

    calculateU1U2(u, tn);
    //saveToExcel(u, tn);
    //saveToImage(u, tn);
    //saveToTextF(u, tn);
    //if (ln%2 == 0) printf("%4d min: %8.6f max: %8.6f\n", ln, u.min(), u.max());
}

Problem0HCommon::Problem0HCommon() {}

Problem0HCommon::~Problem0HCommon() {}

//--------------------------------------------------------------------------------------------------------------//

auto Problem0HFunctional::setDimension(const Dimension &timeDimension, const Dimension &dimensionX, const Dimension &dimensionY) -> void
{
    Problem0HForward::setTimeDimension(timeDimension);
    Problem0HForward::addSpaceDimension(dimensionX);
    Problem0HForward::addSpaceDimension(dimensionY);

    Problem0HBckward::setTimeDimension(timeDimension);
    Problem0HBckward::addSpaceDimension(dimensionX);
    Problem0HBckward::addSpaceDimension(dimensionY);

    const unsigned int L = static_cast<unsigned int> ( timeDimension.size() );
    const unsigned int N = static_cast<unsigned int> ( dimensionX.size() );
    const unsigned int M = static_cast<unsigned int> ( dimensionY.size() );

    p_sigmaX = 0.04;
    p_sigmaY = 0.04;
    p_sigmaT = 0.04;

    c_sigmaX = N/100;
    c_sigmaY = M/100;

    //U1.resize(M+1, N+1, 0.0);
    //U2.resize(M+1, N+1, 0.0);
    u1.resize(M+1, N+1, 0.0);
    u2.resize(M+1, N+1, 0.0);

    ksi = SpacePoint(0.25, 0.25);

    optimalParameters.resize(source_number);

    for (unsigned int sn=0; sn<source_number; sn++)
    {
        optimalParameters[sn].create(timeDimension, dimensionX, dimensionY);
    }
}

auto Problem0HFunctional::fx(const DoubleVector &x) const -> double
{
    vectorToParameter(x);

    DoubleMatrix u;
    forward(u, a, gamma);
    double sum = 0.0;
    sum += epsilon1 * integral1(Problem0HForward::u1);
    //    sum += epsilon2 * integral2(Problem0HForward::u2);
    sum += epsilon2 * integral1(Problem0HForward::u2);
    //sum += norm();
    //sum += penalty();
    return sum;
}

auto Problem0HFunctional::integral1(const DoubleMatrix &u) const -> double
{
    //const Dimension &time = Problem0HForward::timeDimension();
    const Dimension &dimX = Problem0HForward::spaceDimension(Dimension::DimensionX);
    const Dimension &dimY = Problem0HForward::spaceDimension(Dimension::DimensionY);
    //const unsigned int L = static_cast<unsigned int> ( time.size() );
    const unsigned int N = static_cast<unsigned int> ( dimX.size() );
    const unsigned int M = static_cast<unsigned int> ( dimY.size() );
    //const double ht = time.step();
    const double hx = dimX.step();
    const double hy = dimY.step();

    double udiff = 0.0;
    double usum = 0.0;

    udiff = u[0][0]; usum += 0.25 * udiff * udiff;// * mu(0, 0);
    udiff = u[0][N]; usum += 0.25 * udiff * udiff;// * mu(N, 0);
    udiff = u[M][0]; usum += 0.25 * udiff * udiff;// * mu(0, M);
    udiff = u[M][N]; usum += 0.25 * udiff * udiff;// * mu(N, M);

    for (unsigned int n=1; n<=N-1; n++)
    {
        udiff = u[0][n]; usum += 0.5 * udiff * udiff;// * mu(n, 0);
        udiff = u[M][n]; usum += 0.5 * udiff * udiff;// * mu(n, M);
    }

    for (unsigned int m=1; m<=M-1; m++)
    {
        udiff = u[m][0]; usum += 0.5 * udiff * udiff;// * mu(0, m);
        udiff = u[m][N]; usum += 0.5 * udiff * udiff;// * mu(N, m);
    }

    for (unsigned int m=1; m<=M-1; m++)
    {
        for (unsigned int n=1; n<=N-1; n++)
        {
            udiff = u[m][n]; usum += udiff * udiff;// * mu(n, m);
        }
    }

    return usum*(hx*hy);
}

auto Problem0HFunctional::integral2(const DoubleMatrix &u) const -> double
{
    //const Dimension &time = Problem0HForward::timeDimension();
    const Dimension &dimX = Problem0HForward::spaceDimension(Dimension::DimensionX);
    const Dimension &dimY = Problem0HForward::spaceDimension(Dimension::DimensionY);
    //const unsigned int L = static_cast<unsigned int> ( time.size() );
    const unsigned int N = static_cast<unsigned int> ( dimX.size() );
    const unsigned int M = static_cast<unsigned int> ( dimY.size() );
    //const double ht = time.step();
    const double hx = dimX.step();
    const double hy = dimY.step();

    double udiff = 0.0;
    double usum = 0.0;

    udiff = u[0][0]; usum += 0.25 * udiff * udiff;// * mu(0, 0);
    udiff = u[0][N]; usum += 0.25 * udiff * udiff;// * mu(N, 0);
    udiff = u[M][0]; usum += 0.25 * udiff * udiff;// * mu(0, M);
    udiff = u[M][N]; usum += 0.25 * udiff * udiff;// * mu(N, M);

    for (unsigned int n=1; n<=N-1; n++)
    {
        udiff = u[0][n]; usum += 0.5 * udiff * udiff;// * mu(n, 0);
        udiff = u[M][n]; usum += 0.5 * udiff * udiff;// * mu(n, M);
    }

    for (unsigned int m=1; m<=M-1; m++)
    {
        udiff = u[m][0]; usum += 0.5 * udiff * udiff;// * mu(0, m);
        udiff = u[m][N]; usum += 0.5 * udiff * udiff;// * mu(N, m);
    }

    for (unsigned int m=1; m<=M-1; m++)
    {
        for (unsigned int n=1; n<=N-1; n++)
        {
            udiff = u[m][n]; usum += udiff * udiff;// * mu(n, m);
        }
    }

    return usum*(hx*hy);
}

auto Problem0HFunctional::norm() const -> double { return 0.0; }

auto Problem0HFunctional::penalty() const -> double { return 0.0; }

auto Problem0HFunctional::forward(DoubleMatrix &u, double a, double gamma) const -> void
{
    Problem0HForward::implicit_calculate_D2V1(u, a, gamma);
}

auto Problem0HFunctional::backward(DoubleMatrix &p, double a, double gamma)  const -> void
{
    Problem0HBckward::implicit_calculate_D2V1(p, a, gamma);
}

auto Problem0HFunctional::gradient(const DoubleVector &x, DoubleVector &g) const -> void
{
    vectorToParameter(x);

    g.clear(); g.resize(x.length());

    const Dimension &time = Problem0HForward::timeDimension();
    //const Dimension &dimX = Problem0HForward::spaceDimension(Dimension::DimensionX);
    //const Dimension &dimY = Problem0HForward::spaceDimension(Dimension::DimensionY);
    const unsigned int L = static_cast<unsigned int>(time.size());
    //const unsigned int N = static_cast<unsigned int>(dimX.size());
    //const unsigned int M = static_cast<unsigned int>(dimY.size());
    //const double hx = dimX.step();
    //const double hy = dimY.step();
    const double ht = time.step();

    DoubleMatrix u, p;
    forward(u, a, gamma);
    backward(p, a, gamma);

    unsigned int length = 2*L+1;
    for (unsigned int sn=0; sn<source_number; sn++)
    {
        const Problem0HParameter &optimalParameter = optimalParameters[sn];

        unsigned int offset = sn*length;
        g[offset+0] = 0.0;
        for (unsigned int ln=1; ln<=2*L; ln++)
        {
            g[offset+ln] = -optimalParameter.psi_vl[ln];
        }

        //---------------------------------------------------------------------------//

        const unsigned int point_offset = source_number*length;
        const unsigned int ix = point_offset+2*sn+0;
        const unsigned int iy = point_offset+2*sn+1;

        unsigned int ln = 0;
        g[ix] = 0.0;
        g[iy] = 0.0;
        ln = 0;
        g[ix] += 0.5*optimalParameter.psi_dx[ln] * optimalParameter.pwr_vl[ln];
        g[iy] += 0.5*optimalParameter.psi_dy[ln] * optimalParameter.pwr_vl[ln];
        for (unsigned int i=1; i<=L-1; i++)
        {
            ln = 2*i;
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

//--------------------------------------------------------------------------------------------------------------//



auto Problem0HForward::initial(const SpaceNodePDE &, InitialCondition) const -> double { return 0.0; }

auto Problem0HForward::boundary(const SpaceNodePDE &, const TimeNodePDE &) const -> double { return 0.0; }

auto Problem0HForward::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double
{
    double pv = p(sn, tn);

    unsigned int ln = static_cast<unsigned int>(tn.i);
    double pulse1 = optimalParameters[0].pwr_vl[ln] * optimalParameters[0].deltaGrid.weight(sn);
    double pulse2 = optimalParameters[1].pwr_vl[ln] * optimalParameters[1].deltaGrid.weight(sn);

    return pv + pulse1 + pulse2;
}

auto Problem0HForward::p(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double
{
    //static const double sigmaX = 8.0*Problem0HForward::spaceDimension(Dimension::DimensionX).step();
    //static const double sigmaY = 8.0*Problem0HForward::spaceDimension(Dimension::DimensionY).step();
    const double ht = timeDimension().step();

    //static const double alpha1 = 1.0/(2.0*M_PI*p_sigmaX*p_sigmaY);
    //static const double alpha2 = 1.0/(2.0*p_sigmaX*p_sigmaY);
    //static const double alpha3 = 100.0;
    //double pv = 5.0 * alpha1 * exp( -alpha2 * ((sn.x-ksi.x)*(sn.x-ksi.x)+(sn.y-ksi.y)*(sn.y-ksi.y)) - alpha3*tn.t);
    //return pv;

    const double factor1 = 1.0/(2.0*M_PI*p_sigmaX*p_sigmaY);
    const double sigmax2 = 1.0/(2.0*p_sigmaX*p_sigmaX);
    const double sigmay2 = 1.0/(2.0*p_sigmaY*p_sigmaY);
    //const double factor2 = 2.0/(sqrt(2.0*M_PI)*p_sigmaT);
    //const double sigmat2 = 1.0/(2.0*p_sigmaT*p_sigmaT);

    double a, b;

    a = factor1 * exp(-(sigmax2*(sn.x-ksi.x)*(sn.x-ksi.x)+sigmay2*(sn.y-ksi.y)*(sn.y-ksi.y)));
    //if (abs(tn.t - 0.01) >= 0.0) b = factor2 * exp(-(sigmat2*(tn.t-0.01)*(tn.t-0.01))); else b = 0.0;
    if (abs(tn.t - 0.01) == 0.0) b = 2.0/ht; else b = 0.0;
    return a*b;
}


auto Problem0HForward::calculateU1U2(const DoubleMatrix &u, const TimeNodePDE &tn) const -> void
{
    const Dimension &time = Problem0HForward::timeDimension();
    const unsigned int L = static_cast<unsigned int>(time.size());
    const double ht = time.step();

    Problem0HForward *forward = const_cast<Problem0HForward*>(this);
    //if (ln == 2*(L-0)) { forward->u1 = u; }

    //if (ln == 2*(L-2)) { forward->u2 *= +0.0;}
    //if (ln == 2*(L-1)) { forward->u2  = -1.0*u;}
    //if (ln == 2*(L-0)) { forward->u2 += +1.0*u; forward->u2 *= (1.0/ht);}

    const Dimension &dimX = Problem0HForward::spaceDimension(Dimension::DimensionX);
    const Dimension &dimY = Problem0HForward::spaceDimension(Dimension::DimensionY);
    const unsigned int N = static_cast<unsigned int>(dimX.size());
    const unsigned int M = static_cast<unsigned int>(dimY.size());

    if (tn.i == 2*(L-1))
    {
        for (unsigned int m=0; m<=M; m++)
        {
            for (unsigned int n=0; n<=N; n++)
            {
                forward->u2[m][n] = u[m][n];
            }
        }

        //IPrinter::printMatrix(u);
        //IPrinter::printSeperatorLine();

        //printf("%d\n", ln);
    }

    if (tn.i == 2*(L-0))
    {
        //IPrinter::printSeperatorLine();
        //IPrinter::printMatrix(forward->u2);
        //IPrinter::printSeperatorLine();
        for (unsigned int m=0; m<=M; m++)
        {
            for (unsigned int n=0; n<=N; n++)
            {
                forward->u1[m][n] = u[m][n];
                double prev = forward->u2[m][n];
                forward->u2[m][n] = (u[m][n] - prev)/ht;
            }
        }

        //IPrinter::printMatrix(u);
        //IPrinter::printSeperatorLine();
        //IPrinter::printMatrix(forward->u2);
        //IPrinter::printSeperatorLine();

        //printf("%d\n", ln);
        //IPrinter::printMatrix(forward->u1);
        //IPrinter::printSeperatorLine();
        //IPrinter::printMatrix(forward->u2);
        //IPrinter::printSeperatorLine();
    }


    //    if (ln == 2*(L-2)) { forward->u2  = -1.0*u; }
    //    if (ln == 2*(L-0)) { forward->u2 += +1.0*u; forward->u2 *= (0.5/ht); }

    //    if (ln == 2*(L-2)) { forward->u2  = +1.0*u; }
    //    if (ln == 2*(L-1)) { forward->u2 += -4.0*u; }
    //    if (ln == 2*(L-0)) { forward->u2 += +3.0*u; forward->u2 *= 1.0/(2.0*ht); IPrinter::printMatrix(u2); }

    //    if (ln == 2*(L-3)) { forward->u2  = -2.0*u; }
    //    if (ln == 2*(L-2)) { forward->u2 += +9.0*u; }
    //    if (ln == 2*(L-1)) { forward->u2 += -18.0*u; }
    //    if (ln == 2*(L-0)) { forward->u2 += +11.0*u; forward->u2 *= +(1.0/(6.0*ht)); IPrinter::printMatrix(u2); }

    //    if (ln == 2*(L-4)) { forward->u2  = +3.0*u; }
    //    if (ln == 2*(L-3)) { forward->u2 += -16.0*u; }
    //    if (ln == 2*(L-2)) { forward->u2 += +36.0*u; }
    //    if (ln == 2*(L-1)) { forward->u2 += -48.0*u; }
    //    if (ln == 2*(L-0)) { forward->u2 += +25.0*u; forward->u2 *= +(1.0/(12.0*ht)); /*IPrinter::printMatrix(u2);*/ }
}

auto Problem0HForward::saveToExcel(const DoubleMatrix &u UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM) const -> void
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

auto Problem0HForward::saveToImage(const DoubleMatrix &u UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM) const -> void
{
#ifdef USE_LIB_IMAGING
    //const Dimension &time = Problem0HForward::timeDimension();
    const Dimension &dimX = Problem0HForward::spaceDimension(Dimension::DimensionX);
    const Dimension &dimY = Problem0HForward::spaceDimension(Dimension::DimensionY);
    //const unsigned int L = static_cast<unsigned int>(time.size());
    const unsigned int N = static_cast<unsigned int>(dimX.size());
    const unsigned int M = static_cast<unsigned int>(dimY.size());
    //const double ht = time.step();
    //const double hx = dimX.step();
    //const double hy = dimY.step();

    QString filename = QString("data/problem0H/c/png/f_%1.png").arg(tn.i, 4, 10, QChar('0'));
    QPixmap pixmap;
    visualizeMatrixHeat(u, u.min(), u.max(), pixmap, N+1, M+1);
    pixmap.save(filename);
#endif
}

auto Problem0HForward::saveToTextF(const DoubleMatrix &u UNUSED_PARAM, const TimeNodePDE &tn) const -> void
{
    QString filename = QString("data/problem0H/f/txt/f_%1.txt").arg(tn.i, 4, 10, QChar('0'));
    IPrinter::print(u,filename.toLatin1().data());
    printf("Forward: %4d %0.3f %10.8f %10.8f %4d %4d\n", tn.i, tn.t, u.min(), u.max(), 0, 0);
}

//--------------------------------------------------------------------------------------------------------------//

auto Problem0HBckward::layerInfo(const DoubleMatrix &p, const TimeNodePDE &tn) const -> void
{
    saveBackwardInformarion(p, tn);
    //saveToImage(p, ln);
}

auto Problem0HBckward::initial(const SpaceNodePDE &sn, InitialCondition condition) const -> double
{
    unsigned int n = static_cast<unsigned int>(sn.i);
    unsigned int m = static_cast<unsigned int>(sn.j);
    if (condition == InitialCondition::InitialValue) { return -2.0*epsilon2*(u2[m][n]/*-U2[m][n]*/); }
    if (condition == InitialCondition::FirstDerivative) { return +2.0*epsilon1*(u1[m][n]/*-U1[m][n]*/)
                /*+ gamma*initial(sn, InitialCondition::InitialValue)*/; }
    throw std::exception();
    return NAN;
}

auto Problem0HBckward::boundary(const SpaceNodePDE&, const TimeNodePDE&) const -> double { return 0.0; }

auto Problem0HBckward::f(const SpaceNodePDE&, const TimeNodePDE&) const -> double { return 0.0; }

auto Problem0HBckward::saveBackwardInformarion(const DoubleMatrix &p, const TimeNodePDE &tn) const -> void
{
    Problem0HBckward* const_this = const_cast<Problem0HBckward*>(this);

    for (unsigned int sn=0; sn<source_number; sn++)
    {
        double psi_vl, psi_dx, psi_dy;
        psi_vl = optimalParameters[sn].deltaGrid.lumpPointGauss(p);
        optimalParameters[sn].deltaGrid.lumpPointGauss(p, psi_dx, psi_dy);

        const_this->optimalParameters[sn].psi_vl[tn.i] = psi_vl;
        const_this->optimalParameters[sn].psi_dx[tn.i] = psi_dx;
        const_this->optimalParameters[sn].psi_dy[tn.i] = psi_dy;
    }
}

auto Problem0HBckward::saveToImage(const DoubleMatrix &p, const TimeNodePDE &tn) const -> void
{
#ifdef USE_LIB_IMAGING
    //const Dimension &time = Problem0HBckward::timeDimension();
    const Dimension &dimX = Problem0HBckward::spaceDimension(Dimension::DimensionX);
    const Dimension &dimY = Problem0HBckward::spaceDimension(Dimension::DimensionY);
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

auto Problem0HFunctional::vectorToParameter(const DoubleVector &x) const -> void
{    
    const Dimension &time = Problem0HForward::timeDimension();
    const Dimension &dimX = Problem0HForward::spaceDimension(Dimension::DimensionX);
    const Dimension &dimY = Problem0HForward::spaceDimension(Dimension::DimensionY);
    const unsigned int L = static_cast<unsigned int> ( time.size() );
    const unsigned int N = static_cast<unsigned int> ( dimX.size() );
    const unsigned int M = static_cast<unsigned int> ( dimY.size() );

    Problem0HFunctional* const_this = const_cast<Problem0HFunctional*>(this);
    const_this->optimalParameters.resize(source_number);

    const unsigned int length = 2*L+1;
    for (unsigned int sn=0; sn<source_number; sn++)
    {
        Problem0HParameter &parameter = const_this->optimalParameters[sn];
        parameter.create(time, dimX, dimY);
        for (unsigned int ln=0; ln<length; ln++) parameter.pwr_vl[ln] = x[sn*length+ln];
        parameter.p.x = x[source_number*length+2*sn+0];
        parameter.p.y = x[source_number*length+2*sn+1];
        parameter.distribute(parameter.p);
        parameter.deltaGrid.distributeGauss(parameter.p, c_sigmaX, c_sigmaY);
    }
}

auto Problem0HFunctional::parameterToVector(DoubleVector &x) const -> void
{
    const Dimension &time = Problem0HForward::timeDimension();
    //const Dimension &dimX = Problem0HForward::spaceDimension(Dimension::DimensionX);
    //const Dimension &dimY = Problem0HForward::spaceDimension(Dimension::DimensionY);
    const unsigned int L = static_cast<unsigned int> ( time.size() );
    //const unsigned int N = static_cast<unsigned int> ( dimX.size() );
    //const unsigned int M = static_cast<unsigned int> ( dimY.size() );
    //const double ht = time.step();
    //const double hx = dimX.step();
    //const double hy = dimY.step();

    unsigned int size = (2*L+3)*source_number;
    x.resize(size);
    unsigned int length = 2*L+1;

    for (unsigned int sn=0; sn<source_number; sn++)
    {
        Problem0HParameter &parameter = const_cast<Problem0HFunctional*>(this)->optimalParameters[sn];
        for (unsigned int ln=0; ln<length; ln++)
        {
            x[sn*length+ln] = parameter.pwr_vl[ln];
        }
        //x[sn*length+0] = 0.0;
        //x[sn*length+1] = 0.0;
        //x[sn*length+(length-1)] = 0.0;
        x[source_number*length+2*sn+0] = parameter.p.x;
        x[source_number*length+2*sn+1] = parameter.p.y;
    }
}
