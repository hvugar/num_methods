#include "solver.h"

using namespace p3p;

void Solver::frw_calculate() const
{
    //puts("forward->solveInitialValueProblem...");
    //for (unsigned int i=0; i<=TIME_MAX; i++) printf("%4.2f ", externalSource[i].q); puts("");
    //for (unsigned int i=0; i<=TIME_MAX; i++) printf("%4.2f ", externalSource[i].v); puts("");
    forward->solveInitialValueProblem(ODESolverMethod::RUNGE_KUTTA_2);
//    for (unsigned int i=0; i<=TIME_MAX; i++)
//        printf("%u %14.10f %14.10f\n", i, externalSource[i].z.x, externalSource[i].z.y);

    //for (unsigned int i=0; i<=TIME_MAX; i++) printf("%4.2f ", externalSource[i].z.x); puts("");
    //for (unsigned int i=0; i<=TIME_MAX; i++) printf("%4.2f ", externalSource[i].z.y); puts("");
    //puts("forward->solveInitialValueProblem.");
    //puts("forward->implicit_calculate_D2V...");
    forward->implicit_calculate_D2V1();
    //puts("forward->implicit_calculate_D2V.");
}

void Solver::bcw_calculate() const
{
    //puts("backward->implicit_calculate_D2V1...");
    backward->implicit_calculate_D2V1();
    //puts("backward->implicit_calculate_D2V1.");
    //puts("backward->solveFinalValueProblem...");
    backward->solveFinalValueProblem(ODESolverMethod::RUNGE_KUTTA_2);
    //puts("backward->solveFinalValueProblem.");
}

void Solver::Main(int argc, char **argv)
{
    const unsigned int time_size = TIME_MAX+1;
    const double ht = TIME_STEP;

    Solver solver(Dimension(TIME_STEP, 0, TIME_MAX), Dimension(DIMX_STEP, 0, DIMX_MAX), Dimension(DIMY_STEP, 0, DIMX_MAX));
//    for (unsigned int i=0; i<=TIME_MAX; i++)
//    {
//        double t = i*ht;
//        solver.externalSource[i].z.x = 0.8*t*t+0.1;
//        solver.externalSource[i].z.y = 0.8*t*t+0.1;
//        //solver.externalSource[i].q = 1.0;
//    }
//        solver.frw_calculate();
    //    return;

    //solver.e = SpacePoint(0.2, 0.3);
    DoubleVector x(2*time_size);
    for (unsigned int i=0; i<=TIME_MAX; i++)
    {
        x[i] = 100.2;
        x[i+time_size] = 1.6;
    }

    unsigned int s1 = 0, f1 = TIME_MAX, s2 = TIME_MAX+1, f2 = 2*TIME_MAX+1;
    IPrinter::printVector(8, 4, x.mid(s1, f1));
    IPrinter::printVector(8, 4, x.mid(s2, f2));
    IPrinter::printSeperatorLine();

    // printing analitic gradients
    DoubleVector ga;
//    solver.vectorToParameter(x);
//    solver.frw_calculate();
//    return;
    solver.gradient(x, ga);
    ga[s1] = 0.0;
    ga[f1] = 0.0;
    ga[s2] = 0.0;
    ga[f2] = 0.0;
    DoubleVector gn(x.length());
    IPrinter::printVector(10, 4, ga.mid(s1, f1)/*.EuclideanNormalize()*/);
    IGradient::Gradient(&solver, 0.010, x, gn, s1, f1);
    gn[s1] = 0.0;
    gn[f1] = 0.0;
    IPrinter::printVector(10, 4, gn.mid(s1, f1)/*.EuclideanNormalize()*/);

    IPrinter::printVector(10, 4, ga.mid(s2, f2)/*.EuclideanNormalize()*/);
    IGradient::Gradient(&solver, 0.01, x, gn, s2, f2);
    gn[s2] = 0.0;
    gn[f2] = 0.0;
    IPrinter::printVector(10, 4, gn.mid(s2, f2)/*.EuclideanNormalize()*/);

    //    DoubleVector ga1(11);
    //    //DoubleVector ga2(11);
    //    for (unsigned int i=0*time_size, j=0; i<1*time_size; i++) { if ((i-0)%((TIME_MAX)/10)==0) { ga1[j] = ga[i]; j++; } }
    //    //for (unsigned int i=0*time_size, j=0; i<2*time_size; i++) { if ((i-1)%((TIME_MAX)/10)==0) { ga2[j] = ga[i]; j++; } }
    //    ga1[00] = 0.0;
    //    ga1[10] = 0.0;
    //    //ga2[00] = ga2[10] = 0.0;
    //    IPrinter::printVector(14, 10, ga1.EuclideanNormalize());
    //    //IPrinter::printVector(14, 10, ga2.EuclideanNormalize());
    //    IPrinter::printSeperatorLine();

    //    // printing numerical gradients
    //    DoubleVector gn(x.length());
    //    IGradient::Gradient(&solver, 0.010, x, gn);
    //    DoubleVector gn1(11);
    ////    DoubleVector gn2(11);
    //    for (unsigned int i=0*time_size, j=0; i<1*time_size; i++)
    //    {
    //        if ((i-0)%((TIME_MAX)/10)==0)
    //        {
    //            IGradient::Gradient(&solver, 0.010, x, gn, i, i); gn1[j] = gn[i]; j++;
    //        }
    //    }
    ////    for (unsigned int i=1*time_size, j=0; i<2*time_size; i++)
    ////    {
    ////        if ((i-0)%((TIME_MAX)/10)==0)
    ////        {
    ////            IGradient::Gradient(&solver, 0.010, x, gn, i, i); gn2[j] = gn[i]; j++;
    ////        }
    ////    }
    //    gn1[00] = 0.0;
    //    gn1[10] = 0.0;
    //    //gn2[00] = gn2[10] = 0.0;
    //    //IPrinter::printVector(14, 10, gn.EuclideanNormalize());
    //    IPrinter::printVector(14, 10, gn1.EuclideanNormalize());
    //    //IPrinter::printVector(14, 10, gn2.EuclideanNormalize());
    //    IPrinter::printSeperatorLine();
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

Solver::Solver(const Dimension &timeDimension, const Dimension &spaceDimensionX, const Dimension &spaceDimensionY)
{
    setDimension(timeDimension, spaceDimensionX, spaceDimensionY);

    forward = new HeatEquationIBVP(this);
    backward = new HeatEquationFBVP(this);

    forward->setThermalDiffusivity(0.1);
    forward->setThermalConvection(0.0);
    forward->setThermalConductivity(0.0);

    backward->setThermalDiffusivity(-0.1);
    backward->setThermalConvection(0.0);
    backward->setThermalConductivity(0.0);

    measurePoints[0] = {0.25, 0.25, 0.0};
    measurePoints[1] = {0.25, 0.75, 0.0};
    measurePoints[2] = {0.75, 0.75, 0.0};
    measurePoints[3] = {0.75, 0.75, 0.0};
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
    externalSource.resize(timeDimension().size());

    V.clear(); V.resize(spaceDimensionX().size(), spaceDimensionY().size(), 8.5);
    U.clear(); U.resize(spaceDimensionX().size(), spaceDimensionY().size(), 0.0);

    deltaZ.cleanGrid();
    deltaZ.initGrid(spaceDimensionX(), spaceDimensionY());
    deltaZ.resetAll();

    if (measurePoints) delete [] measurePoints;
    measurePoints = new SpacePoint[measurePointNumber];

    for (SpacePoint* route : heatSourceRoutes) {
        if (route) delete route;
    }
    heatSourceRoutes.clear();

    for (size_t i=0; i<heatSourceNumber; i++)
    {
        SpacePoint* route = new SpacePoint[timeDimension().size()];
        heatSourceRoutes.push_back(route);
    }
}

Solver::Solver(const Solver &)
{}

Solver::~Solver()
{
}

//--------------------------------------------------------------------------------------------------------------//

auto Solver::setDimension(const Dimension &timeDimension, const Dimension &spaceDimensionX, const Dimension &spaceDimensionY) -> void
{
    this->_timeDimension = timeDimension;
    this->_spaceDimensionX = spaceDimensionX;
    this->_spaceDimensionY = spaceDimensionY;
    validate();
}

auto Solver::fx(const DoubleVector &x) const -> double
{
    vectorToParameter(x);

    frw_calculate();

    double sum = integral(V);
    //sum += norm();
    //sum += penalty();
    return sum;
}

auto Solver::integral(const DoubleMatrix &) const -> double
{
    const Dimension &dimensionX = spaceDimensionX();
    const Dimension &dimensionY = spaceDimensionY();
    const unsigned int N = static_cast<unsigned int> ( dimensionX.size()-1 );
    const unsigned int M = static_cast<unsigned int> ( dimensionY.size()-1 );
    const double hx = dimensionX.step();
    const double hy = dimensionY.step();

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

    //for (unsigned int i=0; i<=TIME_MAX; i++) printf("%4.2f ", externalSource[i].q); puts("");
    //for (unsigned int i=0; i<=TIME_MAX; i++) printf("%4.2f ", externalSource[i].v); puts("");

    g.clear(); g.resize(x.length());

    frw_calculate();
    //for (unsigned int i=0; i<=TIME_MAX; i++) printf("%4.2f ", externalSource[i].z.x); puts("");
    //for (unsigned int i=0; i<=TIME_MAX; i++) printf("%4.2f ", externalSource[i].z.y); puts("");
    bcw_calculate();

    const Dimension &time = timeDimension();
    const unsigned int time_size = time.size();
    const double ht = time.step();

    for (unsigned int ln=0; ln<time_size; ln++)
    {
        PointNodeODE node(ln*ht, ln);
        const ExternalSource &es = externalSource[ln];
        g[ln] = -es.psi.z;
        g[ln+time_size] = -( es.phi.x*forward->C(node, 1) + es.phi.y*forward->C(node, 2) );
    }
}

auto Solver::project(DoubleVector &x, unsigned int index) -> void
{}

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
    //cn = BoundaryConditionPDE(BoundaryCondition::Dirichlet, 1.0, 0.0, 1.0);
    cn = BoundaryConditionPDE(BoundaryCondition::Neumann, 0.0, 1.0, 1.0);
    return 0.0;
}

auto Solver::frw_f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double
{
    unsigned int ln = static_cast<unsigned int>(tn.i) / 2;
    double q = externalSource[ln].q;
    const SpacePoint &z = externalSource[ln].z;
    return q * deltaZ.gaussWeight(sn, z, SIGMA_X, SIGMA_Y);

    //    const Dimension &dimensionX = spaceDimensionX();
    //    const Dimension &dimensionY = spaceDimensionY();
    //    const double hx = dimensionX.step();
    //    const double hy = dimensionY.step();

    //    if (sn.i==50 && sn.j==50)
    //        return mq[ln]*(1.0/(hx*hy));
    //    else
    //        return 0.0;

}

void Solver::frw_layerInfo(const DoubleMatrix &u, const TimeNodePDE &tn) const
{
    int ln = tn.i / 2;
    if (ln == timeDimension().max())
    {
        const_cast<Solver*>(this)->U = u;
        //IPrinter::printSeperatorLine(std::to_string(ln).data(), '-');
        //IPrinter::printMatrix(9, 5, U);
        //IPrinter::printSeperatorLine();
        //printf("min: %14.8f max: %14.8f\n", U.min(), U.max());
    }

    frw_saveToImage(u, tn);
}

void Solver::frw_saveToImage(const DoubleMatrix &u UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM) const
{
#ifdef USE_LIB_IMAGING
    static double MIN = +10000.0;
    static double MAX = -10000.0;
    if (u.max()>MAX) MAX = u.max();
    if (u.min()<MIN) MIN = u.min();
    printf(">>> %6d %14.6f %14.6f %14.6f %14.6f\n", tn.i/2, u.min(), u.max(), MIN, MAX);

    //if (tn.i%10==0)
    {
        QString filename = QString("data/problem3P/f/png/%1.png").arg(tn.i/2, 8, 10, QChar('0'));
        QPixmap pixmap;
        visualizeMatrixHeat(u, 0.0, 6.052, pixmap, DIMX_MAX+1, DIMY_MAX+1);
        pixmap.save(filename);
    }
#endif
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
    return -2.0*(U[m][n]-V[m][n]);
}

double Solver::bcw_boundary(const SpaceNodePDE&, const TimeNodePDE&, BoundaryConditionPDE &cn) const
{
    //cn = BoundaryConditionPDE(BoundaryCondition::Dirichlet, 1.0, 0.0, 1.0);
    cn = BoundaryConditionPDE(BoundaryCondition::Neumann, 0.0, 1.0, 1.0);
    return 0.0;
}

double Solver::bcw_f(const SpaceNodePDE&, const TimeNodePDE&) const { return 0.0; }

void Solver::bcw_layerInfo(const DoubleMatrix &p, const TimeNodePDE &tn) const
{
    bcw_saveBackwardInformarion(p, tn);
    //bcw_saveToImage(p, tn);
}

void Solver::bcw_saveBackwardInformarion(const DoubleMatrix &p, const TimeNodePDE &tn) const
{
    Solver* const_this = const_cast<Solver*>(this);
    unsigned int ln = static_cast<unsigned int>(tn.i) / 2;

    const SpacePoint &z = externalSource[ln].z;
    double psi_vl =0.0, psi_dx = 0.0, psi_dy = 0.0;

//    unsigned int rx = static_cast<unsigned int>(round(z.x/(spaceDimensionX().step())));
//    unsigned int ry = static_cast<unsigned int>(round(z.y/(spaceDimensionY().step())));
//    printf(">>> %d %d\n", rx, ry);
//    psi_dx = (p[ry][rx+1]-p[ry][rx-1])/(2.0*spaceDimensionX().step());
//    psi_dy = (p[ry+1][rx]-p[ry-1][rx])/(2.0*spaceDimensionY().step());
//    psi_vl = p[ry][rx];

    //const_this->deltaZ.p() = z;

    //psi_vl = const_this->deltaZ.consentrateInPoint(p, psi_dx, psi_dy, 5);

    const_this->deltaZ.distributeGauss(z, NODE_PER_SIGMA_X, NODE_PER_SIGMA_Y);
    psi_vl = const_this->deltaZ.lumpPointGauss(p, psi_dx, psi_dy);
    const_this->externalSource[ln].psi = SpacePoint(psi_dx, psi_dy, psi_vl);
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

auto HeatEquationIBVP::weight() const -> double { throw std::runtime_error(""); }

auto HeatEquationIBVP::timeDimension() const -> Dimension { return solver->timeDimension(); }

auto HeatEquationIBVP::spaceDimensionX() const -> Dimension { return solver->spaceDimensionX(); }

auto HeatEquationIBVP::spaceDimensionY() const -> Dimension { return solver->spaceDimensionY(); }

auto HeatEquationIBVP::spaceDimensionZ() const -> Dimension { throw std::runtime_error(""); }

//--------------------------------------------------------------------------------------------------------------//

auto HeatEquationIBVP::A(const PointNodeODE &, unsigned int row, unsigned int col) const -> double
{
#ifdef VARIANT_1
    const double mx[2][2] = {
        {-3.0, -2.0}, {-4.0, -5.0}
    };
    return mx[row-1][col-1];
#endif
#ifdef VARIANT_2
    const double vl[4][4] = {
        {+0.0, +0.0, +1.0, +0.0},
        {+0.0, +0.0, +0.0, +1.0},
        {+0.0, +0.0, +0.0, +0.0},
        {+0.0, +0.0, +0.0, +0.0},
    };
    return vl[row-1][col-1];
#endif
#ifdef VARIANT_3
    const double vl[2][2] = { {+0.0, +0.0}, {+0.0, +0.0} };
    return vl[row-1][col-1];
#endif
}

auto HeatEquationIBVP::B(const PointNodeODE &node, unsigned int row) const -> double
{
    //return solver->zt(node, row)
    //        - (A(node, row, 1)*solver->z(node, 1)+A(node, row, 2)*solver->z(node, 2))
    //        - C(node, row)*solver->v(node)
    //        + C(node, row)*solver->mv[node.i];

    double t = node.x;
#ifdef VARIANT_1
    double v = solver->v(node);
    return (row == 1) ?
                0.4*M_PI*sin(2.0*M_PI*t) - 1.8*M_PI*t*sin(4.0*M_PI*t*t)
                + 3.0*(0.40*sin(1.0*M_PI*t)*sin(1.0*M_PI*t)+0.45*cos(2.0*M_PI*t*t)*cos(2.0*M_PI*t*t)+0.05)
                + 2.0*(0.30*sin(4.0*M_PI*t)*sin(4.0*M_PI*t)+0.45*cos(3.0*M_PI*t*t)*cos(3.0*M_PI*t*t)+0.05) - 5.0*sin(2.0*M_PI*t*t) + 5.0*v:
                1.2*M_PI*sin(8.0*M_PI*t) - 2.7*M_PI*t*sin(6.0*M_PI*t*t)
                + 4.0*(0.40*sin(1.0*M_PI*t)*sin(1.0*M_PI*t)+0.45*cos(2.0*M_PI*t*t)*cos(2.0*M_PI*t*t)+0.05)
                + 5.0*(0.30*sin(4.0*M_PI*t)*sin(4.0*M_PI*t)+0.45*cos(3.0*M_PI*t*t)*cos(3.0*M_PI*t*t)+0.05) - 4.0*sin(2.0*M_PI*t*t) + 4.0*v;
#endif
#ifdef VARIANT_2
    puts("B1");
    //const double mx[4] = {+0.0, +0.0, +1.0, +1.0};
    unsigned int ln = static_cast<unsigned int>(node.i);
    double aa = C(node, row-1);//*solver->externalSource[ln].v;
    puts("B1");
    return aa;
#endif
#ifdef VARIANT_3
    const double vl[2] = {+0.0, +0.0};
    unsigned int ln = static_cast<unsigned int>(node.i);
    double ret = C(node, row)*solver->externalSource[ln].v + vl[row-1];
    return ret;
#endif
}

auto HeatEquationIBVP::C(const PointNodeODE &node, unsigned int row) const -> double
{
#ifdef VARIANT_1
    const double c[2] = { +5.0, +4.0 };
    return c[row-1];
#endif
#ifdef VARIANT_2
    puts("C1");
    const double vl[4] = { +0.0, +0.0, +1.0, +1.0 };
    double bb = 0.0;//vl[row-1];
    printf("%f row %d\n", bb, row);
    return bb;
#endif
#ifdef VARIANT_3
    const double vl[2] = {+1.0, +1.0};
    return vl[row-1]*node.x;
#endif
}

auto HeatEquationIBVP::initial(InitialCondition, unsigned int row) const -> double
{
#ifdef VARIANT_1
    const double val[2] = { 0.50, 0.50 };
    return val[row-1];
#endif
#ifdef VARIANT_2
    const double val[4] = { 0.10, 0.10, 0.00, 0.00 };
    return val[row-1];
#endif
#ifdef VARIANT_3
    const double val[2] = { 0.10, 0.10 };
    return val[row-1];
#endif
}

auto HeatEquationIBVP::count() const -> unsigned int
{
#ifdef VARIANT_1
    return 2;
#endif
#ifdef VARIANT_2
    return 4;
#endif
#ifdef VARIANT_3
    return 2;
#endif
}

auto HeatEquationIBVP::dimension() const -> Dimension { return solver->timeDimension(); }

auto HeatEquationIBVP::iterationInfo(const DoubleVector &z, const PointNodeODE &node) const -> void
{
    unsigned int ln = static_cast<unsigned int>(node.i);
    ExternalSource &es = const_cast<Solver*>(solver)->externalSource[ln];
    es.z.x = z[0];
    es.z.y = z[1];
}

//--------------------------------------------------------------------------------------------------------------//

HeatEquationFBVP::HeatEquationFBVP(Solver *solver) { this->solver = solver; }

HeatEquationFBVP::~HeatEquationFBVP() {}

auto HeatEquationFBVP::final(const SpaceNodePDE &sn, FinalCondition condition) const -> double { return solver->bcw_final(sn, condition); }

auto HeatEquationFBVP::boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const -> double { return solver->bcw_boundary(sn, tn, condition); }

auto HeatEquationFBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double { return solver->bcw_f(sn, tn); }

auto HeatEquationFBVP::layerInfo(const DoubleMatrix &p, const TimeNodePDE &tn) const -> void { return solver->bcw_layerInfo(p, tn); }

auto HeatEquationFBVP::weight() const -> double { throw std::runtime_error(""); }

auto HeatEquationFBVP::timeDimension() const -> Dimension { return solver->timeDimension(); }

auto HeatEquationFBVP::spaceDimensionX() const -> Dimension { return solver->spaceDimensionX(); }

auto HeatEquationFBVP::spaceDimensionY() const -> Dimension { return solver->spaceDimensionY(); }

auto HeatEquationFBVP::spaceDimensionZ() const -> Dimension { throw std::runtime_error(""); }

//--------------------------------------------------------------------------------------------------------------//

auto HeatEquationFBVP::A(const PointNodeODE &, unsigned int row, unsigned int col) const -> double
{
#ifdef VARIANT_1
    const double mx[2][2] = { {+3.0, +2.0}, { +4.0, +5.0 } };
    return mx[row-1][col-1];
#endif
#ifdef VARIANT_2
    const double mx[4][4] = {
        {+0.0, +0.0, +0.0, +0.0},
        {+0.0, +0.0, +0.0, +0.0},
        {+0.0, +0.0, +0.0, +0.0},
        {+0.0, +0.0, +0.0, +0.0},
    };
    return mx[row-1][col-1];
#endif
#ifdef VARIANT_3
    const double vl[2][2] = { {+0.0, +0.0}, {+0.0, +0.0} };
    return -vl[row-1][col-1];
#endif
}

auto HeatEquationFBVP::B(const PointNodeODE &node, unsigned int row) const -> double
{
#ifdef VARIANT_1
    unsigned int ln = static_cast<unsigned int>(node.i);
    ExternalSource &es = const_cast<Solver*>(solver)->externalSource[ln];
    return -es.q * (row == 1 ? es.psi.x : es.psi.y);
#endif
#ifdef VARIANT_2
    unsigned int ln = static_cast<unsigned int>(node.i);
    ExternalSource &es = const_cast<Solver*>(solver)->externalSource[ln];
    return - es.q * (row == 1 ? es.psi.x : row == 2 ? es.psi.y : 0.0);
#endif
#ifdef VARIANT_3
    unsigned int ln = static_cast<unsigned int>(node.i);
    ExternalSource &es = const_cast<Solver*>(solver)->externalSource[ln];
    if (row == 1) return -es.q * es.psi.x;
    if (row == 2) return -es.q * es.psi.y;
    throw std::runtime_error("OK");
#endif
}

auto HeatEquationFBVP::final(FinalCondition, unsigned int) const -> double
{
    return 0.0;
}

auto HeatEquationFBVP::count() const -> unsigned int
{
#ifdef VARIANT_1
    return 2;
#endif
#ifdef VARIANT_2
    return 4;
#endif
#ifdef VARIANT_3
    return 2;
#endif
}

auto HeatEquationFBVP::dimension() const -> Dimension { return solver->timeDimension(); }

auto HeatEquationFBVP::iterationInfo(const DoubleVector &phi, const PointNodeODE &node) const -> void
{
    unsigned int ln = static_cast<unsigned int>(node.i);
    ExternalSource &es = const_cast<Solver*>(solver)->externalSource[ln];
    es.phi.x = phi[0];
    es.phi.y = phi[1];
    //if (node.i == 100) { es.psi.x = es.psi.y = 0.0; es.phi.x = es.phi.y = 0.0; }
    //printf("%4d %10.6f %10.6f %10.6f %10.6f\n", ln, es.phi.x, es.phi.y, es.psi.x, es.psi.y);
}

//--------------------------------------------------------------------------------------------------------------//

auto Solver::vectorToParameter(const DoubleVector &x) const -> void
{
    const Dimension &time = timeDimension();
    const unsigned int time_size = time.size();
    for (unsigned int ln=0; ln<time_size; ln++)
    {
        const_cast<Solver*>(this)->externalSource[ln].q = x[ln];
        const_cast<Solver*>(this)->externalSource[ln].v = x[ln+time_size];
    }
}

auto Solver::parameterToVector(DoubleVector &x) const -> void
{
    const Dimension &time = timeDimension();
    const unsigned int time_size = time.size();
    x.clear();
    //x.resize(time_size);
    x.resize(2*time_size);
    for (unsigned int ln=0; ln<time_size; ln++)
    {
        x[ln] = externalSource[ln].q;
        x[ln+time_size] = externalSource[ln].v;
    }
}

auto Solver::v(const PointNodeODE &node) const -> double
{
    return sin(2.0*M_PI*node.x*node.x);
}

auto Solver::zt(const PointNodeODE &node, unsigned int row) const -> double
{
    const double t1 = node.x, t2 = node.x*node.x;
    const double res[2] = { 0.40*M_PI*sin(2.0*M_PI*t1) - 1.80*M_PI*sin(4.0*M_PI*t2)*t1,
                            1.20*M_PI*sin(8.0*M_PI*t1) - 2.70*M_PI*sin(6.0*M_PI*t2)*t1 };
    return res[row-1];
}

auto Solver::z(const PointNodeODE &node, unsigned int row) const -> double
{
    const double t1 = node.x, t2 = node.x*node.x;
    const double res[2] = { 0.40*sin(1.0*M_PI*t1)*sin(1.0*M_PI*t1) + 0.45*cos(2.0*M_PI*t2)*cos(2.0*M_PI*t2) + 0.05,
                            0.30*sin(4.0*M_PI*t1)*sin(4.0*M_PI*t1) + 0.45*cos(3.0*M_PI*t2)*cos(3.0*M_PI*t2) + 0.05 };
    return res[row-1];
}

//--------------------------------------------------------------------------------------------------------------//
