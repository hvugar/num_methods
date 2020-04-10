﻿#include "solver.h"

using namespace p3p;

void Solver::frw_calculate() const
{
    std::vector<DoubleVector> rv;
    forward->solveInitialValueProblem(rv, ODESolverMethod::RUNGE_KUTTA_4);
    rv[0].clear();
    rv[1].clear();
    forward->implicit_calculate_D2V1();
    rv.clear();
}
void Solver::bcw_calculate() const
{
    backward->implicit_calculate_D2V1();
    std::vector<DoubleVector> rv;
    backward->solveFinalValueProblem(rv, ODESolverMethod::RUNGE_KUTTA_4);
    rv[0].clear();
    rv[1].clear();
    rv.clear();
}

void Solver::Main(int argc, char **argv)
{
#ifdef USE_LIB_IMAGING
    QGuiApplication app(argc, argv);
#endif
    Solver solver(Dimension(TIME_STEP, 0, TIME_MAX), Dimension(DIMX_STEP, 0, DIMX_MAX), Dimension(DIMY_STEP, 0, DIMX_MAX));
    int time_max = TIME_MAX;
    unsigned int time_size = TIME_MAX+1;
    double ht = TIME_STEP;

    DoubleVector x(2*time_size);
    for (unsigned int i=0; i<=TIME_MAX; i++)
    {
        x[0*time_size+i] = 1.0;
        x[1*time_size+i] = solver.v(PointNodeODE(i*TIME_STEP, static_cast<int>(i)));
    }

    // printing analitic gradients
    DoubleVector ga;
    solver.gradient(x, ga);
    DoubleVector ga1(11);
    DoubleVector ga2(11);
    for (unsigned int i=0*time_size, j=0; i<1*time_size; i++) { if ((i-0)%((TIME_MAX)/10)==0) { ga1[j] = ga[i]; j++; } }
    for (unsigned int i=0*time_size, j=0; i<2*time_size; i++) { if ((i-1)%((TIME_MAX)/10)==0) { ga2[j] = ga[i]; j++; } }
    ga1[00] = ga2[00] = 0.0;
    ga1[10] = ga2[10] = 0.0;
    IPrinter::printVector(14, 10, ga1.EuclideanNormalize());
    IPrinter::printVector(14, 10, ga2.EuclideanNormalize());
    IPrinter::printSeperatorLine();

    // printing numerical gradients
    DoubleVector gn(x.length());
    DoubleVector gn1(11, 0.0);IPrinter::printSeperatorLine();
    DoubleVector gn2(11);
    IPrinter::printSeperatorLine();
    for (unsigned int i=0*time_size, j=0; i<1*time_size; i++)
    {
        IPrinter::printSeperatorLine();
        if ((i-0)%((TIME_MAX)/10)==0)
        {
            IGradient::Gradient(&solver, 0.010, x, gn, i, i); gn1[j] = gn[i]; j++;
        }
        IPrinter::printSeperatorLine();
    }
    for (unsigned int i=1*time_size, j=0; i<2*time_size; i++) { if ((i-0)%((TIME_MAX)/10)==0) { IGradient::Gradient(&solver, 0.010, x, gn, i, i); gn2[j] = gn[i]; j++; } }
    gn1[00] = gn2[00] = 0.0;
    gn1[10] = gn2[10] = 0.0;
    IPrinter::printVector(14, 10, gn1.EuclideanNormalize());
    IPrinter::printVector(14, 10, gn2.EuclideanNormalize());
    IPrinter::printSeperatorLine();
}

void Solver::example1()
{
    //#ifdef USE_LIB_IMAGING
    //    QGuiApplication app(argc, argv);
    //#endif
    //    IPrinter::printDateTime();
    //    //Solver slv(Dimension(0.005, 0, 200), Dimension(0.01, 100, 200), Dimension(0.01, 200, 300));

    //    int time_max = TIME_MAX;
    //    double ht = TIME_STEP;

    //    Solver solver(Dimension(TIME_STEP, 0, TIME_MAX), Dimension(DIMX_STEP, 0, DIMX_MAX), Dimension(DIMY_STEP, 0, DIMX_MAX));

    //    HeatEquationIBVP forward(&solver);
    //    std::vector<DoubleVector> rv1;
    //    std::vector<DoubleVector> rv2;
    //    forward.solveInitialValueProblem(rv1, IFirstOrderLinearODE::ODESolverMethod::EULER);
    //    forward.solveInitialValueProblem(rv2, IFirstOrderLinearODE::ODESolverMethod::RUNGE_KUTTA_4);

    //    DoubleVector x(TIME_MAX+1);

    //    for (unsigned int i=0; i<=TIME_MAX; i++)
    //    {
    //        solver.mq[i] = 1.0;
    //        solver.mz[i] = SpacePoint( solver.z(PointNodeODE(ht*i, i), 1)+0.03, solver.z(PointNodeODE(ht*i, i), 2) );
    //        x[i] = 0.0;//sin(2.0*M_PI*(i*TIME_STEP)*(i*TIME_STEP));
    //        solver.mv[i] = x[i];
    //    }

    //    for (unsigned int i=0; i<=time_max; i++) { if (i%(time_max/10)==0) { printf("%10.6f ", rv1[i][0]); } } puts("");
    //    for (unsigned int i=0; i<=time_max; i++) { if (i%(time_max/10)==0) { printf("%10.6f ", rv2[i][0]); } } puts("");
    //    for (unsigned int i=0; i<=time_max; i++) { if (i%(time_max/10)==0) { printf("%10.6f ", forward.z(PointNodeODE(ht*i, i), 1)); }} puts("");
    //    for (unsigned int i=0; i<=time_max; i++) { if (i%(time_max/10)==0) { printf("%10.6f ", rv1[i][1]); } } puts("");
    //    for (unsigned int i=0; i<=time_max; i++) { if (i%(time_max/10)==0) { printf("%10.6f ", rv2[i][1]); } } puts("");
    //    for (unsigned int i=0; i<=time_max; i++) { if (i%(time_max/10)==0) { printf("%10.6f ", forward.z(PointNodeODE(ht*i, i), 2)); }} puts("");

    //    for (unsigned int i=0; i<=time_max; i++) { if (i%(time_max/10)==0) { printf("%10.6f ", rv1[i][0]); } } puts("");

    //    forward.setThermalDiffusivity(1.0);
    //    forward.setThermalConvection(0.0);
    //    forward.setThermalConductivity(0.0);
    //    std::cout << "Computing forward problem..." << std::endl;
    //    forward.implicit_calculate_D2V1();
    //    std::cout << "Computed forward problem." << std::endl;

    //    HeatEquationFBVP backward(&solver);
    //    backward.setThermalDiffusivity(-1.0);
    //    backward.setThermalConvection(0.0);
    //    backward.setThermalConductivity(0.0);
    //    std::cout << "Computing backward problem..." << std::endl;
    //    backward.implicit_calculate_D2V1();
    //    std::vector<DoubleVector> r;
    //    backward.solveFinalValueProblem(r);
    //    for (unsigned int i=0; i<=time_max; i++) { solver.psi[i].x = r[i][0]; if (i%(time_max/10)==0) { printf("%10.6f ", r[i][0]); } } puts("");
    //    for (unsigned int i=0; i<=time_max; i++) { solver.psi[i].y = r[i][1]; if (i%(time_max/10)==0) { printf("%10.6f ", r[i][1]); } } puts("");
    //    std::cout << "Computed backward problem." << std::endl;

    //    DoubleVector gx(TIME_MAX+1); for (unsigned int i=1; i<time_max; i++) gx[i] = -(backward.B(PointNodeODE(ht*i, i), 1)*solver.psi[i].x
    //                                                                                   +backward.B(PointNodeODE(ht*i, i), 2)*solver.psi[i].y);
    //    gx.L2Normalize();
    //    //gy.L2Normalize();
    //    for (unsigned int i=0; i<=time_max; i++) { if (i%(time_max/10)==0) { printf("%10.6f ", gx[i]); } } puts("");
    //    //for (unsigned int i=0; i<=time_max; i++) { if (i%(time_max/10)==0) { printf("%10.6f ", gy[i]); } } puts("");


    //    // printing numerical gradients
    //    DoubleVector gn(401);
    //    DoubleVector gn1(11);

    //    for (unsigned int i=0, j=0; i<=400; i++) { if ((i-0)%((401-1)/10)==0) { IGradient::Gradient(&solver, 0.010, x, gn, i, i); gn1[j] = gn[i]; j++; } }
    //    gn1[00] = 0.0;
    //    gn1[10] = 0.0;
    //    IPrinter::printVector(10, 6, gn1);
    //    IPrinter::printSeperatorLine();
}

Solver::Solver()
{
}

Solver::Solver(const Dimension &timeDimension, const Dimension &spaceDimensionX, const Dimension &spaceDimensionY)
{
    setDimension(timeDimension, spaceDimensionX, spaceDimensionY);
    forward = new HeatEquationIBVP(this);
    backward = new HeatEquationFBVP(this);

    forward->setThermalDiffusivity(1.0);
    forward->setThermalConvection(0.0);
    forward->setThermalConductivity(0.0);

    backward->setThermalDiffusivity(-1.0);
    backward->setThermalConvection(0.0);
    backward->setThermalConductivity(0.0);
}

void Solver::setTimeDimension(const Dimension &timeDimension)
{
    this->_timeDimension = timeDimension;
    validate();
}

void Solver::setSpaceDimensionX(const Dimension &spaceDimensionX)
{
    this->_spaceDimensionX = spaceDimensionX;
    validate();
}

void Solver::setSpaceDimensionY(const Dimension &spaceDimensionY)
{
    this->_spaceDimensionY = spaceDimensionY;
    validate();
}

void Solver::validate()
{
    mq.clear(); mq.resize(timeDimension().size());
    mz.clear(); mz.resize(timeDimension().size());
    V.clear(); V.resize(spaceDimensionX().size(), spaceDimensionY().size());
    U.clear(); U.resize(spaceDimensionX().size(), spaceDimensionY().size());

    deltaZ.cleanGrid();
    deltaZ.initGrid(spaceDimensionX(), spaceDimensionY());
    deltaZ.resetAll();

    psi.resize(timeDimension().size());
    phi.resize(timeDimension().size());

    mv.resize(timeDimension().size());
}

Solver::Solver(const Solver &)
{}

Solver::~Solver()
{
}

//--------------------------------------------------------------------------------------------------------------//

auto Solver::setDimension(const Dimension &timeDimension, const Dimension &spaceDimensionX, const Dimension &spaceDimensionY) -> void
{
    _timeDimension = timeDimension;
    _spaceDimensionX = spaceDimensionX;
    _spaceDimensionY = spaceDimensionY;
    validate();
}

auto Solver::fx(const DoubleVector &x) const -> double
{
    vectorToParameter(x);

    frw_calculate();

    double sum = epsilon * integral(V);
    //sum += norm();
    //sum += penalty();
    return sum;
}

auto Solver::integral(const DoubleMatrix &) const -> double
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

    udiff = (U[0][0]-V[0][0]); usum += 0.25 * udiff * udiff;// * mu(0, 0);
    udiff = (U[0][N]-V[0][N]); usum += 0.25 * udiff * udiff;// * mu(N, 0);
    udiff = (U[M][0]-V[M][0]); usum += 0.25 * udiff * udiff;// * mu(0, M);
    udiff = (U[M][N]-V[M][N]); usum += 0.25 * udiff * udiff;// * mu(N, M);

    for (unsigned int n=1; n<=N-1; n++)
    {
        udiff = U[0][n]-V[0][n]; usum += 0.5 * udiff * udiff;// * mu(n, 0);
        udiff = U[M][n]-V[M][n]; usum += 0.5 * udiff * udiff;// * mu(n, M);
    }

    for (unsigned int m=1; m<=M-1; m++)
    {
        udiff = U[m][0]-V[m][0]; usum += 0.5 * udiff * udiff;// * mu(0, m);
        udiff = U[m][N]-V[m][N]; usum += 0.5 * udiff * udiff;// * mu(N, m);
    }

    for (unsigned int m=1; m<=M-1; m++)
    {
        for (unsigned int n=1; n<=N-1; n++)
        {
            udiff = U[m][n]-V[m][n]; usum += udiff * udiff;// * mu(n, m);
        }
    }

    return usum*(hx*hy);
}

auto Solver::norm() const -> double { return 0.0; }

auto Solver::penalty() const -> double { return 0.0; }

auto Solver::gradient(const DoubleVector &x, DoubleVector &g) const -> void
{
    vectorToParameter(x);

    IPrinter::printVector(mq);
    IPrinter::printVector(mv);

    g.clear(); g.resize(x.length());

    std::cout << "Calculating forward..." << std::endl;
    frw_calculate();
    std::cout << "Calculated forward." << std::endl;
    std::cout << "Calculating backward..." << std::endl;
    bcw_calculate();
    std::cout << "Calculated backward." << std::endl;

    const Dimension &time = timeDimension();
    const unsigned int time_size = time.size();
    const double ht = time.step();

//    for (unsigned int i=0; i<time_size; i++)
//    {
//        printf("%d | %f %f | %f %f | %f %f %f\n", i, mz[i].x, mz[i].y, phi[i].x, phi[i].y, psi[i].x, psi[i].y, psi[i].z);
//    }

    for (unsigned int i=0; i<time_size; i++)
    {
        PointNodeODE node(i*ht, i);
        g[0*time_size+i] = -mq[i] * ( phi[i].x*forward->B(node, 1) + phi[i].y*forward->B(node, 2) );
        g[1*time_size+i] = -psi[i].z;
    }
}

auto Solver::project(DoubleVector &x, unsigned int index) -> void
{
}

auto Solver::project(DoubleVector &) const  -> void {}

auto Solver::print(unsigned int i, const DoubleVector &x, const DoubleVector &g, double f, double alpha, GradientMethod::MethodResult result) const -> void
{
    //    C_UNUSED(i); C_UNUSED(x); C_UNUSED(g); C_UNUSED(f); C_UNUSED(alpha); C_UNUSED(result);
    //    const char* msg = nullptr; C_UNUSED(msg);
    //    if (result == GradientMethod::MethodResult::BREAK_FIRST_ITERATION)    msg = "BREAK_FIRST_ITERATION   ";
    //    if (result == GradientMethod::MethodResult::FIRST_ITERATION)          msg = "FIRST_ITERATION         ";
    //    if (result == GradientMethod::MethodResult::BREAK_GRADIENT_NORM_LESS) msg = "BREAK_GRADIENT_NORM_LESS";
    //    if (result == GradientMethod::MethodResult::BREAK_DISTANCE_LESS)      msg = "BREAK_DISTANCE_LESS     ";
    //    if (result == GradientMethod::MethodResult::NEXT_ITERATION)           msg = "NEXT_ITERATION          ";

    //    ProblemSolver* fw = const_cast<ProblemSolver*>(this);
    //    fw->vectorToParameter(x);
    //    frw_calculate();

    //    const Dimension &time = timeDimension();
    //    const unsigned int time_size = time.size();
    //    const unsigned int L = static_cast<unsigned int>(time.size()-1);
    //    const unsigned int o = (2*L+1)*source_number;
    //    const unsigned int offset = source_number*time_size;

    //    printf("I[%3d]: F:%.6f I:%.6f P:%.6f N:%.5f R:%.3f e:%.3f a:%10.6f | ", i, f, fx(x), 0.0, 0.0, 0.0, 0.0, alpha);
    //    printf("%12.8f %12.8f | %12.8f %12.8f | ", u1.min(), u1.max(), u2.min(), u2.max());
    //    printf("eta: %8.6f %8.6f %8.6f %8.6f\n", x[offset+0], x[offset+1], x[offset+2], x[offset+3]);

    //    const unsigned int s1 = 0*time_size; const unsigned int f1 = 1*time_size-1;
    //    const unsigned int s2 = 1*time_size; const unsigned int f2 = 2*time_size-1;
    //    const unsigned int s3 = 2*time_size; const unsigned int f3 = 2*time_size+3;
    //    DoubleVector v1(time_size);
    //    DoubleVector v2(time_size);
    //    for (unsigned int i=s1, j=0; i<=f1; i++, j++) { v1[j] = x[i]; }
    //    for (unsigned int i=s2, j=0; i<=f2; i++, j++) { v2[j] = x[i]; }
    //    IPrinter::printVector(v1);
    //    IPrinter::printVector(v2);
    //    IPrinter::printSeperatorLine("-");
}

//--------------------------------------------------------------------------------------------------------------//

double Solver::frw_initial(const SpaceNodePDE &, InitialCondition) const { return 0.0; }

double Solver::frw_boundary(const SpaceNodePDE &, const TimeNodePDE &, BoundaryConditionPDE & cn) const
{
    cn = BoundaryConditionPDE(BoundaryCondition::Neumann, 0.0, 1.0, 1.0);
    return 0.0;
}

auto Solver::frw_f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double
{
    const unsigned int ln = static_cast<unsigned int>(tn.i/2);
    const double factor1 = 1.0/(2.0*M_PI*SIGMA_X*SIGMA_Y);
    const double sigmax2 = 1.0/(2.0*SIGMA_X*SIGMA_X);
    const double sigmay2 = 1.0/(2.0*SIGMA_Y*SIGMA_Y);

    const SpacePoint pnt = mz.at(ln);
    const double a = factor1 * exp(-(sigmax2*(sn.x-pnt.x)*(sn.x-pnt.x)+sigmay2*(sn.y-pnt.y)*(sn.y-pnt.y)));
    return mq.at(ln)*a;
}

void Solver::frw_layerInfo(const DoubleMatrix &u, const TimeNodePDE &tn) const
{
    if (tn.i/2 == timeDimension().max()) { const_cast<Solver*>(this)->U = u; }

    //frw_saveToImage(u, tn);
}

void Solver::frw_saveToExcel(const DoubleMatrix &u UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM) const
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

void Solver::frw_saveToImage(const DoubleMatrix &u UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM) const
{
#ifdef USE_LIB_IMAGING
    static double MIN = +10000.0;
    static double MAX = -10000.0;
    if (u.max()>MAX) MAX = u.max();
    if (u.min()<MIN) MIN = u.min();
    printf(">>> %6d %14.6f %14.6f %14.6f %14.6f\n", tn.i/2, u.min(), u.max(), MIN, MAX);

    if (tn.i%10==0)
    {
        QString filename = QString("data/problem3P/f/png/%1.png").arg(tn.i/10, 8, 10, QChar('0'));
        QPixmap pixmap;
        visualizeMatrixHeat(u, 0.0, 4.70254, pixmap, DIMX_MAX+1, DIMY_MAX+1);
        pixmap.save(filename);
    }
#endif
}

void Solver::frw_saveToTextF(const DoubleMatrix &u UNUSED_PARAM, const TimeNodePDE &tn) const
{
#ifdef USE_LIB_TEXT
    static double MIN = +100000.0;
    static double MAX = -100000.0;

    std::string txt_number = std::to_string(tn.i);
    std::string filename = std::string("data/problem0H/f/txt/") + std::string(8 - txt_number.length(), '0') + txt_number + std::string(".txt");
    IPrinter::print(u, filename.c_str());
    if (MIN > u.min()) MIN = u.min();
    if (MAX < u.max()) MAX = u.max();
    printf("Forward: %4d %6.3f | %12.8f %12.8f | %12.8f %12.8f | %4d %4d\n", tn.i, tn.t, u.min(), u.max(), MIN, MAX, 0, 0);
#endif
}

//--------------------------------------------------------------------------------------------------------------//

double Solver::bcw_final(const SpaceNodePDE &sn, FinalCondition) const
{
    unsigned int n = static_cast<unsigned int>(sn.i);
    unsigned int m = static_cast<unsigned int>(sn.j);
    double aa = -2.0*(U[m][n]-V[m][n]);
    return aa;
}

double Solver::bcw_boundary(const SpaceNodePDE&, const TimeNodePDE&, BoundaryConditionPDE &cn) const
{
    cn = BoundaryConditionPDE(BoundaryCondition::Neumann, 0.0, 1.0, 1.0);
    return 0.0;
}

double Solver::bcw_f(const SpaceNodePDE&, const TimeNodePDE&) const { return 0.0; }

void Solver::bcw_layerInfo(const DoubleMatrix &p, const TimeNodePDE &tn) const
{
    bcw_saveBackwardInformarion(p, tn);
    bcw_saveToImage(p, tn);
}

void Solver::bcw_saveBackwardInformarion(const DoubleMatrix &p, const TimeNodePDE &tn) const
{
    Solver* const_this = const_cast<Solver*>(this);
    unsigned int ln = static_cast<unsigned int>(tn.i) / 2;

    double psi_vl, psi_dx, psi_dy;

    const_this->deltaZ.resetAll();
    const_this->deltaZ.distributeGauss(mz[ln], 1, 1);
    psi_vl = const_this->deltaZ.lumpPointGauss(p, psi_dx, psi_dy);
    const_this->psi[ln] = SpacePoint(psi_dx, psi_dy, psi_vl);
}

void Solver::bcw_saveToImage(const DoubleMatrix &p, const TimeNodePDE &tn) const
{
#ifdef USE_LIB_IMAGING
    static double MIN = +10000.0;
    static double MAX = -10000.0;
    if (p.max()>MAX) MAX = p.max();
    if (p.min()<MIN) MIN = p.min();
    printf(">>> %6d %14.6f %14.6f %14.6f %14.6f\n", tn.i/2, p.min(), p.max(), MIN, MAX);

    if (tn.i%10==0)
    {
        QString filename = QString("data/problem3P/b/png/%1.png").arg(tn.i/10, 8, 10, QChar('0'));
        QPixmap pixmap;
        //        visualizeMatrixHeat(p, -3.0334, 0.8191, pixmap, DIMX_MAX+1, DIMY_MAX+1);
        visualizeMatrixHeat(p, p.min(), p.max(), pixmap, DIMX_MAX+1, DIMY_MAX+1);
        pixmap.save(filename);
    }
#endif
}

//--------------------------------------------------------------------------------------------------------------//

HeatEquationIBVP::HeatEquationIBVP(Solver *solver) { this->solver = solver; }

HeatEquationIBVP::~HeatEquationIBVP() {}

auto HeatEquationIBVP::initial(const SpaceNodePDE &sn, InitialCondition cn) const -> double { return solver->frw_initial(sn, cn); }

auto HeatEquationIBVP::boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &cn) const -> double { return solver->frw_boundary(sn, tn, cn); }

auto HeatEquationIBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const-> double { return solver->frw_f(sn, tn); }

auto HeatEquationIBVP::layerInfo(const DoubleMatrix &u, const TimeNodePDE &tn) const -> void { return solver->frw_layerInfo(u, tn); }

auto HeatEquationIBVP::weight() const -> double { return 0.5; }

auto HeatEquationIBVP::timeDimension() const -> Dimension { return solver->timeDimension(); }

auto HeatEquationIBVP::spaceDimensionX() const -> Dimension { return solver->spaceDimensionX(); }

auto HeatEquationIBVP::spaceDimensionY() const -> Dimension { return solver->spaceDimensionY(); }

auto HeatEquationIBVP::spaceDimensionZ() const -> Dimension { return Dimension(DIMY_STEP, 0, DIMY_MAX); }

//--------------------------------------------------------------------------------------------------------------//

auto HeatEquationIBVP::A(const PointNodeODE &, unsigned int row, unsigned int col) const -> double
{
    const double mx[2][2] = { {-3.0, -2.0}, { -4.0, -5.0 } };
    return mx[row-1][col-1];
}

auto HeatEquationIBVP::B(const PointNodeODE &node, unsigned int row) const -> double
{
    return solver->zt(node, row)
            - (A(node, row, 1)*solver->z(node, 1)+A(node, row, 2)*solver->z(node, 2))
            - C(node, row)*solver->v(node)
            + C(node, row)*solver->mv[node.i];
}

auto HeatEquationIBVP::C(const PointNodeODE &, unsigned int row) const -> double
{
    const double c[2] = { +5.0, +4.0 };
    return c[row-1];
}

auto HeatEquationIBVP::initial(InitialCondition, unsigned int row) const -> double
{
    const double a[2] = { 0.4, 0.3 };
    return a[row-1];
}

auto HeatEquationIBVP::count() const -> unsigned int { return 2; }

auto HeatEquationIBVP::dimension() const -> Dimension { return solver->timeDimension(); }

auto HeatEquationIBVP::iterationInfo(const DoubleVector &v, const PointNodeODE &node) const -> void
{
    const_cast<HeatEquationIBVP*>(this)->solver->mz[static_cast<unsigned int>(node.i)].x = v[0];
    const_cast<HeatEquationIBVP*>(this)->solver->mz[static_cast<unsigned int>(node.i)].y = v[1];
}

//--------------------------------------------------------------------------------------------------------------//

HeatEquationFBVP::HeatEquationFBVP(Solver *solver) { this->solver = solver; }

HeatEquationFBVP::~HeatEquationFBVP() {}

auto HeatEquationFBVP::final(const SpaceNodePDE &sn, FinalCondition condition) const -> double { return solver->bcw_final(sn, condition); }

auto HeatEquationFBVP::boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const -> double { return solver->bcw_boundary(sn, tn, condition); }

auto HeatEquationFBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double { return solver->bcw_f(sn, tn); }

auto HeatEquationFBVP::layerInfo(const DoubleMatrix &p, const TimeNodePDE &tn) const -> void { return solver->bcw_layerInfo(p, tn); }

auto HeatEquationFBVP::weight() const -> double { return 0.5; }

auto HeatEquationFBVP::timeDimension() const -> Dimension { return solver->timeDimension(); }

auto HeatEquationFBVP::spaceDimensionX() const -> Dimension { return solver->spaceDimensionX(); }

auto HeatEquationFBVP::spaceDimensionY() const -> Dimension { return solver->spaceDimensionY(); }

auto HeatEquationFBVP::spaceDimensionZ() const -> Dimension { return Dimension(DIMY_STEP, 0, DIMY_MAX); }

//--------------------------------------------------------------------------------------------------------------//

auto HeatEquationFBVP::A(const PointNodeODE &, unsigned int row, unsigned int col) const -> double
{
    const double mx[2][2] = { {-3.0, -2.0}, { -4.0, -5.0 } };
    return -mx[row-1][col-1];
}

auto HeatEquationFBVP::B(const PointNodeODE &node, unsigned int row) const -> double
{
    return -(row == 1 ? solver->psi[static_cast<unsigned int>(node.i)].x : solver->psi[static_cast<unsigned int>(node.i)].y);
}

auto HeatEquationFBVP::final(FinalCondition, unsigned int) const -> double
{
    return 0.0;
}

auto HeatEquationFBVP::count() const -> unsigned int { return 2; }

auto HeatEquationFBVP::dimension() const -> Dimension { return solver->timeDimension(); }

auto HeatEquationFBVP::iterationInfo(const DoubleVector &v, const PointNodeODE &node) const -> void
{
    const_cast<HeatEquationFBVP*>(this)->solver->phi[static_cast<unsigned int>(node.i)].x = v[0];
    const_cast<HeatEquationFBVP*>(this)->solver->phi[static_cast<unsigned int>(node.i)].y = v[1];
}

//--------------------------------------------------------------------------------------------------------------//

auto Solver::vectorToParameter(const DoubleVector &x) const -> void
{
    const Dimension &time = timeDimension();
    const unsigned int time_size = time.size();
    for (unsigned int i=0; i<time_size; i++)
    {
        const_cast<Solver*>(this)->mq[i] = x[i];
        const_cast<Solver*>(this)->mv[i] = x[time_size+i];
    }
}

auto Solver::parameterToVector(DoubleVector &x) const -> void
{
    const Dimension &time = timeDimension();
    const unsigned int time_size = time.size();
    x.clear();
    x.resize(2*time_size);
    for (unsigned int i=0; i<time_size; i++)
    {
        x[0*time_size+i] = mq[i];
        x[1*time_size+i] = mv[i];
    }
}

auto Solver::v(const PointNodeODE &node) const -> double
{
    return sin(2.0*M_PI*node.x*node.x);
}

auto Solver::zt(const PointNodeODE &node, unsigned int row) const -> double
{
    const double t1 = node.x;
    const double t2 = node.x*node.x;
    const double a[2] = { 0.5*M_PI*sin(2.0*M_PI*t1) - 1.6*M_PI*sin(4.0*M_PI*t2)*t1,
                          2.4*M_PI*sin(8.0*M_PI*t1) - 0.6*M_PI*sin(6.0*M_PI*t1)};
    return a[row-1];
}

auto Solver::z(const PointNodeODE &node, unsigned int row) const -> double
{
    const double t1 = node.x;
    const double t2 = node.x*node.x;
    const double a[2] = { 0.5*sin(1.0*M_PI*t1)*sin(1.0*M_PI*t1) + 0.4*cos(2.0*M_PI*t2)*cos(2.0*M_PI*t2) + 0.05,
                          0.6*sin(4.0*M_PI*t1)*sin(4.0*M_PI*t1) + 0.2*cos(3.0*M_PI*t1)*cos(3.0*M_PI*t1) + 0.1 };
    return a[row-1];
}

//--------------------------------------------------------------------------------------------------------------//

Problem0HParameter& Problem0HParameter::initialize(const Dimension &time, const Dimension &dimX, const Dimension &dimY)
{
#ifdef NEW_FORM
    destroy();
    deltaGrid.initGrid(static_cast<unsigned int>(dimX.size()-1), dimX.step(), static_cast<unsigned int>(dimY.size()-1), dimY.step());
    const unsigned int time_size = time.size();

    pwr_vl.resize(time_size, 0.0);
    psi_vl.resize(time_size, 0.0);
    psi_dx.resize(time_size, 0.0);
    psi_dy.resize(time_size, 0.0);
    psi_x.resize(time_size, 0.0);
    psi_y.resize(time_size, 0.0);
#else
    destroy();

    deltaGrid.initGrid(static_cast<unsigned int>(dimX.size()), dimX.step(), static_cast<unsigned int>(dimY.size()), dimY.step());

    const unsigned int L = static_cast<unsigned int>(time.size());
    const unsigned int length = 2*L-1;
    pwr_vl.resize(length, 0.0);
    psi_vl.resize(length, 0.0);
    psi_dx.resize(length, 0.0);
    psi_dy.resize(length, 0.0);
    psi_x.resize(length, 0.0);
    psi_y.resize(length, 0.0);
#endif

    return *this;
}

void Problem0HParameter::destroy()
{
    deltaGrid.cleanGrid();
    pwr_vl.clear();
    psi_vl.clear();
    psi_dx.clear();
    psi_dy.clear();
    psi_x.clear();
    psi_y.clear();
}

void Problem0HParameter::distribute(const SpacePoint &p)
{
    this->p = p;
    deltaGrid.distributeGauss(p, 1, 1);
}
