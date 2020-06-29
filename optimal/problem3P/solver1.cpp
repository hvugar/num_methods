#include "solver1.h"

using namespace p3p1;

void Solver1::optimize(int argc, char **argv)
{
    QGuiApplication app(argc, argv);

    Solver1 s(Dimension(0.005, 0, 400), Dimension(0.01, 0, 100), Dimension(0.01, 0, 100));
    s.setPointNumber(2, 4);
    s.forward.solver = &s;
    s.backward.solver = &s;

    double a = 0.0001;
    a = 1.0;
    s.forward.setThermalDiffusivity(a);
    s.forward.setThermalConvection(-s._lambda0);
    s.forward.setThermalConductivity(0.0);
    s.drawImage = false;
    s.forward.implicit_calculate_D2V1();
    s.drawImage = false;
    s.V = s.U;
    //return;

    s.backward.setThermalDiffusivity(-a);
    s.backward.setThermalConvection(s._lambda0);
    s.backward.setThermalConductivity(0.0);

    unsigned int size = static_cast<unsigned int>(s.heatSourceNumber*s.measrPointNumber);

    DoubleVector x;
    s.parameterToVector(x);
    IPrinter::print(x.mid(0*size, 1*size-1), x.mid(0*size, 1*size-1).length(), 8, 4);
    IPrinter::print(x.mid(1*size, 2*size-1), x.mid(1*size, 2*size-1).length(), 8, 4);
    IPrinter::print(x.mid(2*size, 3*size-1), x.mid(2*size, 3*size-1).length(), 8, 4);
    IPrinter::print(x.mid(3*size, 4*size-1), x.mid(3*size, 4*size-1).length(), 8, 4);
    IPrinter::print(x.mid(4*size, 5*size-1), x.mid(4*size, 5*size-1).length(), 8, 4);
    IPrinter::print(x.mid(5*size, 6*size-1), x.mid(5*size, 6*size-1).length(), 8, 4);
    IPrinter::print(x.mid(6*size, 7*size-1), x.mid(6*size, 7*size-1).length(), 8, 4);
    IPrinter::print(x.mid(7*size, 8*size-1), x.mid(7*size, 8*size-1).length(), 8, 4);

    for (size_t i=0; i<s.heatSourceNumber; i++)
    {
        for (size_t j=0; j<s.measrPointNumber; j++)
        {
            // alpha
            s.alpha1[i][j] *= +1.1;
            s.alpha2[i][j] *= +1.1;
            s.alpha3[i][j] *= +1.1;

            // betta
            s.betta1[i][j] *= +1.1;
            s.betta2[i][j] *= +1.1;
            s.betta3[i][j] *= +1.1;
//
            s.nomnU1[i][j] *= +1.1;
            s.nomnU2[i][j] *= +1.1;
        }
    }

    s.parameterToVector(x);

    IPrinter::printSeperatorLine();
    IPrinter::print(x.mid(0*size, 1*size-1), x.mid(0*size, 1*size-1).length(), 8, 4);
    IPrinter::print(x.mid(1*size, 2*size-1), x.mid(1*size, 2*size-1).length(), 8, 4);
    IPrinter::print(x.mid(2*size, 3*size-1), x.mid(2*size, 3*size-1).length(), 8, 4);
    IPrinter::print(x.mid(3*size, 4*size-1), x.mid(3*size, 4*size-1).length(), 8, 4);
    IPrinter::print(x.mid(4*size, 5*size-1), x.mid(4*size, 5*size-1).length(), 8, 4);
    IPrinter::print(x.mid(5*size, 6*size-1), x.mid(5*size, 6*size-1).length(), 8, 4);
    IPrinter::print(x.mid(6*size, 7*size-1), x.mid(6*size, 7*size-1).length(), 8, 4);
    IPrinter::print(x.mid(7*size, 8*size-1), x.mid(7*size, 8*size-1).length(), 8, 4);
    IPrinter::printSeperatorLine();

//    ConjugateGradient g;
    SteepestDescentGradient g;
    g.setFunction(&s);
    g.setGradient(&s);
    g.setPrinter(&s);
    g.setProjection(&s);
    //g.setGradientNormalizer(&prob);
    g.setOptimalityTolerance(-0.00001);
    g.setStepTolerance(-0.00001);
    g.setFunctionTolerance(-0.00001);
    g.setR1MinimizeEpsilon(1.0, 0.01);
    g.setNormalize(true);
    g.showExitMessage(true);
    g.setMaxIterationCount(10);
    s.gradMethod = &g;

    g.calculate(x);

    IPrinter::printSeperatorLine(nullptr, '=');

    IPrinter::printSeperatorLine();
    IPrinter::print(x.mid(0*size, 1*size-1), x.mid(0*size, 1*size-1).length(), 8, 4);
    IPrinter::print(x.mid(1*size, 2*size-1), x.mid(1*size, 2*size-1).length(), 8, 4);
    IPrinter::print(x.mid(2*size, 3*size-1), x.mid(2*size, 3*size-1).length(), 8, 4);
    IPrinter::print(x.mid(3*size, 4*size-1), x.mid(3*size, 4*size-1).length(), 8, 4);
    IPrinter::print(x.mid(4*size, 5*size-1), x.mid(4*size, 5*size-1).length(), 8, 4);
    IPrinter::print(x.mid(5*size, 6*size-1), x.mid(5*size, 6*size-1).length(), 8, 4);
    IPrinter::print(x.mid(6*size, 7*size-1), x.mid(6*size, 7*size-1).length(), 8, 4);
    IPrinter::print(x.mid(7*size, 8*size-1), x.mid(7*size, 8*size-1).length(), 8, 4);
    IPrinter::printSeperatorLine();

}

void Solver1::drawImages(Solver1 &solver) const
{
    solver.drawImage = true;
    solver.minU = DBL_MAX;
    solver.maxU = DBL_MIN;

    solver.saveMinMaxU = true;
    solver.forward.implicit_calculate_D2V1();
    solver.saveMinMaxU = false;

    solver.drawImage = true;

    solver.drawImage = false;
}

void Solver1::Main(int argc, char **argv)
{
    optimize(argc, argv); return;

    QGuiApplication app(argc, argv);

    Solver1 s(Dimension(0.005, 0, 200), Dimension(0.01, 0, 100), Dimension(0.01, 0, 100));
    s.setPointNumber(2, 4);
    s.forward.solver = &s;
    s.backward.solver = &s;

    const double thermalDiffusivity = 0.0001;
    const double thermalConductivity = 0.0;
    const double thermalConvection = -s._lambda0;

    s.forward.setThermalDiffusivity(thermalDiffusivity);
    s.forward.setThermalConductivity(thermalConductivity);
    s.forward.setThermalConvection(thermalConvection);
    s.drawImage = false;
    s.forward.implicit_calculate_D2V1();
    s.drawImage = false;
    s.V = s.U;

    s.backward.setThermalDiffusivity(-thermalDiffusivity);
    s.backward.setThermalConductivity(thermalConductivity);
    s.backward.setThermalConvection(-thermalConvection);
    //s.backward.implicit_calculate_D2V1();

    for (size_t i=0; i<s.heatSourceNumber; i++)
    {
        for (size_t j=0; j<s.measrPointNumber; j++)
        {
            s.alpha1[i][j] *= 1.1;
            s.alpha2[i][j] *= 1.1;
            s.alpha3[i][j] *= 1.1;
            s.betta1[i][j] *= 1.1;
            s.betta2[i][j] *= 1.1;
            s.betta3[i][j] *= 1.1;
            s.nomnU1[i][j] *= 1.1;
            s.nomnU2[i][j] *= 1.1;
        }
    }

    DoubleVector x;
    s.parameterToVector(x);
    DoubleVector g;
    s.gradient(x, g);

    //    IPrinter::print(g, g.length());

    size_t size = s.heatSourceNumber*s.measrPointNumber;
    IPrinter::print(g.mid(0*size, 1*size-1).L2Normalize(), g.mid(0*size, 1*size-1).length(), 8, 4);
    IPrinter::print(g.mid(1*size, 2*size-1).L2Normalize(), g.mid(1*size, 2*size-1).length(), 8, 4);
    IPrinter::print(g.mid(2*size, 3*size-1).L2Normalize(), g.mid(2*size, 3*size-1).length(), 8, 4);
    IPrinter::print(g.mid(3*size, 4*size-1).L2Normalize(), g.mid(3*size, 4*size-1).length(), 8, 4);
    IPrinter::print(g.mid(4*size, 5*size-1).L2Normalize(), g.mid(4*size, 5*size-1).length(), 8, 4);
    IPrinter::print(g.mid(5*size, 6*size-1).L2Normalize(), g.mid(5*size, 6*size-1).length(), 8, 4);
    IPrinter::print(g.mid(6*size, 7*size-1).L2Normalize(), g.mid(6*size, 7*size-1).length(), 8, 4);
    IPrinter::print(g.mid(7*size, 8*size-1).L2Normalize(), g.mid(7*size, 8*size-1).length(), 8, 4);

    puts("---------------------------------------------------------------------------------------");
    DoubleVector g1(x.length(), 0.0);
    IGradient::Gradient(&s, 0.01, x, g1, 0*size, 1*size-1);
    IGradient::Gradient(&s, 0.01, x, g1, 1*size, 2*size-1);
    IGradient::Gradient(&s, 0.01, x, g1, 2*size, 3*size-1);
    IGradient::Gradient(&s, 0.01, x, g1, 3*size, 4*size-1);
    IGradient::Gradient(&s, 0.01, x, g1, 4*size, 5*size-1);
    IGradient::Gradient(&s, 0.01, x, g1, 5*size, 6*size-1);
    IGradient::Gradient(&s, 0.01, x, g1, 6*size, 7*size-1);
    IGradient::Gradient(&s, 0.01, x, g1, 7*size, 8*size-1);

    IPrinter::print(g1.mid(0*size, 1*size-1).L2Normalize(), g1.mid(0*size, 1*size-1).length(), 8, 4);
    IPrinter::print(g1.mid(1*size, 2*size-1).L2Normalize(), g1.mid(1*size, 2*size-1).length(), 8, 4);
    IPrinter::print(g1.mid(2*size, 3*size-1).L2Normalize(), g1.mid(2*size, 3*size-1).length(), 8, 4);
    IPrinter::print(g1.mid(3*size, 4*size-1).L2Normalize(), g1.mid(3*size, 4*size-1).length(), 8, 4);
    IPrinter::print(g1.mid(4*size, 5*size-1).L2Normalize(), g1.mid(4*size, 5*size-1).length(), 8, 4);
    IPrinter::print(g1.mid(5*size, 6*size-1).L2Normalize(), g1.mid(5*size, 6*size-1).length(), 8, 4);
    IPrinter::print(g1.mid(6*size, 7*size-1).L2Normalize(), g1.mid(6*size, 7*size-1).length(), 8, 4);
    IPrinter::print(g1.mid(7*size, 8*size-1).L2Normalize(), g1.mid(7*size, 8*size-1).length(), 8, 4);
}

Solver1::Solver1(const Dimension &timeDimension, const Dimension &spaceDimensionX, const Dimension &spaceDimensionY)
{
    this->_timeDimension = timeDimension;
    this->_spaceDimensionX = spaceDimensionX;
    this->_spaceDimensionY = spaceDimensionY;

    V.resize(spaceDimensionX.size(), spaceDimensionY.size());
}

Solver1::~Solver1() {}

void Solver1::setPointNumber(size_t heatSourceNumber, size_t measrPointNumber)
{
    this->heatSourceNumber = heatSourceNumber;
    this->measrPointNumber = measrPointNumber;

    // q
    alpha1.resize(heatSourceNumber, measrPointNumber, 0.200);
    alpha2.resize(heatSourceNumber, measrPointNumber, 0.250);
    //    alpha1.resize(heatSourceNumber, measrPointNumber, 0.000);
    //    alpha2.resize(heatSourceNumber, measrPointNumber, 0.000);
    alpha3.resize(heatSourceNumber, measrPointNumber, 0.275);

    // v
    betta1.resize(heatSourceNumber, measrPointNumber, 0.200);
    betta2.resize(heatSourceNumber, measrPointNumber, 0.250);
    //    betta1.resize(heatSourceNumber, measrPointNumber, 0.000);
    //    betta2.resize(heatSourceNumber, measrPointNumber, 0.000);
    betta3.resize(heatSourceNumber, measrPointNumber, 0.275);

    nomnU1.resize(heatSourceNumber, measrPointNumber, 0.500);
    nomnU2.resize(heatSourceNumber, measrPointNumber, 0.400);

    for (size_t i=0; i<heatSourceNumber; i++)
    {
        for (size_t j=0; j<measrPointNumber; j++)
        {
            // q
            alpha1[i][j] = 0.015;
            alpha2[i][j] = 0.015;
            alpha3[i][j] = 0.015;

            // v
            betta1[i][j] = 0.04;
            betta2[i][j] = 0.03;
            betta3[i][j] = 0.065+(i*0.015);

            nomnU1[i][j] = 0.0;
            nomnU2[i][j] = 0.0;
        }
    }

    this->measurePoints = new SpacePoint[measrPointNumber];
    //this->measurePointValues = new SpacePoint[measrPointNumber];
    this->measurePoints[0] = SpacePoint(0.50, 0.20);
    this->measurePoints[1] = SpacePoint(0.20, 0.50);
    this->measurePoints[2] = SpacePoint(0.50, 0.80);
    this->measurePoints[3] = SpacePoint(0.80, 0.50);

    sourceParams = new ProblemParams[timeDimension().size()];
    for (size_t ln=0; ln<timeDimension().size(); ln++)
    {
        ProblemParams & pp = sourceParams[ln];
        pp.q = new double[heatSourceNumber];
        pp.v = new double[heatSourceNumber];
        //pp.h = new double[heatSourceNumber];
        pp.zv = new SpacePoint[heatSourceNumber];
        pp.zd = new SpacePoint[heatSourceNumber];
        pp.fv = new SpacePoint[heatSourceNumber];
        pp.fd = new SpacePoint[heatSourceNumber];
        pp.ps = new SpacePoint[heatSourceNumber];

        pp.u = new SpacePoint[measrPointNumber];
    }
}

double Solver1::frw_initial(const SpaceNodePDE &, InitialCondition) const { return _initialValue; }

double Solver1::bcw_final(const SpaceNodePDE &sn, FinalCondition) const
{
    size_t i = static_cast<size_t>(sn.i);
    size_t j = static_cast<size_t>(sn.j);
    return -2.0*(U[j][i]-V[j][i]);
}

double Solver1::frw_boundary(const SpaceNodePDE &, const TimeNodePDE &, BoundaryConditionPDE &cn) const
{
    //cn = BoundaryConditionPDE::Dirichlet(1.0, 1.0); return 0.0;
    //cn = BoundaryConditionPDE::Neumann(1.0, 1.0); return 0.0;
    //cn = BoundaryConditionPDE::Robin(lambda, +1.0, lambda); return _environmentTemperature;
    cn = BoundaryConditionPDE::Neumann(1.0, 0.0); return 0.0;
}

double Solver1::bcw_boundary(const SpaceNodePDE &, const TimeNodePDE &, BoundaryConditionPDE &cn) const
{
    //cn = BoundaryConditionPDE::Robin(lambda, -1.0, 0.0); return 0.0;
    cn = BoundaryConditionPDE::Neumann(1.0, 0.0); return 0.0;
}

void Solver1::frw_saveToImage(const DoubleMatrix &u UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM) const
{
#ifdef USE_LIB_IMAGING
    static double MIN = +10000.0;
    static double MAX = -10000.0;
    if (u.max()>MAX) MAX = u.max();
    if (u.min()<MIN) MIN = u.min();
    //printf(">>> %6d %14.6f %14.6f %14.6f %14.6f\n", tn.i, u.min(), u.max(), MIN, MAX);

    //if (tn.i%10==0)
    {
        QString filename = QString("data/problem3P/f/png/%1.png").arg(tn.i, 8, 10, QChar('0'));
        QPixmap pixmap;
        //visualizeMatrixHeat(u, 0.0, 46.052, pixmap, 101, 101);
        visualizeMatrixHeat(u, 1.50, 2.1807, pixmap, spaceDimensionX().size(), spaceDimensionY().size());
//        visualizeMatrixHeat(u, u.min(), u.max(), pixmap, spaceDimensionX().size(), spaceDimensionY().size());
        pixmap.save(filename);
    }
#endif
}

void Solver1::bcw_saveToImage(const DoubleMatrix &p UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM) const
{
#ifdef USE_LIB_IMAGING
    static double MIN = +10000.0;
    static double MAX = -10000.0;
    if (p.max()>MAX) MAX = p.max();
    if (p.min()<MIN) MIN = p.min();
    //printf(">>> %6d %14.6f %14.6f %14.6f %14.6f\n", tn.i, p.min(), p.max(), MIN, MAX);

    //if (tn.i%10==0)
    {
        QString filename = QString("data/problem3P/b/png/%1.png").arg(tn.i, 8, 10, QChar('0'));
        QPixmap pixmap;
        //visualizeMatrixHeat(u, 0.0, 46.052, pixmap, 101, 101);
        visualizeMatrixHeat(p, p.min(), p.max(), pixmap, spaceDimensionX().size(), spaceDimensionY().size());
        pixmap.save(filename);
    }
#endif
}

auto Solver1::bcw_f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double
{
    double fx = 0.0;
    int ln = static_cast<int>(tn.i);
    PointNodeODE node; node.i = ln; node.x = tn.t;
    const ProblemParams &pp = sourceParams[ln];

    for (size_t j=0; j<measrPointNumber; j++)
    {
        const SpacePoint &mp = measurePoints[j];

        double sum = 0.0;
        for (size_t i=0; i<heatSourceNumber; i++)
        {
            const SpacePoint &zi = pp.zv[i];
            if (isPointOnPlate(zi))
            {
                double dist = sqrt((zi.x - mp.x)*(zi.x - mp.x) + (zi.y - mp.y)*(zi.y - mp.y));
                sum += (alpha1[i][j]*dist*dist + alpha2[i][j]*dist + alpha3[i][j]) * ( pp.ps[i].z );
                sum += (betta1[i][j]*dist*dist + betta2[i][j]*dist + betta3[i][j]) * ( A4(node, 1, i+1)*pp.fv[i].x + A4(node, 2, i+1)*pp.fv[i].y );
            }
            else
            {
                //throw std::runtime_error("bcw_f: unknown error...");
            }
        }
        fx -= sum * DeltaFunction::gaussian(sn, mp, SpacePoint(spaceDimensionX().step()*_factor, spaceDimensionY().step()*_factor));
    }

    return fx;
}

bool Solver1::isPointOnPlate(const SpacePoint &z) const
{
    return (0.05 <= z.x && z.x <= 0.95) && (0.05 <= z.y && z.y <= 0.95);
}

void Solver1::bcw_layerInfo(const DoubleMatrix &p, const TimeNodePDE &tn) const
{
    HeatEquationFBVP &const_backward = const_cast<HeatEquationFBVP&>(this->backward);
    int ln = static_cast<int>(tn.i);
    double ht = timeDimension().step();
    DoubleVector f0(4), f1(4);
    PointNodeODE n0, n1;

    if (ln == timeDimension().max())
    {
        n0.i = static_cast<int>(ln); n0.x = n0.i*ht;
        ProblemParams &pp0 = sourceParams[ln];

        for (size_t i=0; i<heatSourceNumber; i++)
        {
            const_backward.i = i+1;

            /**********************************************************************************************************/

            const_backward.start( f0, n0 );
            pp0.fv[i] = SpacePoint(f0[0], f0[1]);
            pp0.fd[i] = SpacePoint(f0[2], f0[3]);

            /**********************************************************************************************************/

            const SpacePoint &zi = pp0.zv[i];
            if (isPointOnPlate(zi))
            {
                SpacePoint d;
                pp0.ps[i].z = DeltaFunction::lumpedPoint4(p, zi, spaceDimensionX(), spaceDimensionY(), d);
                pp0.ps[i].x = d.x; pp0.ps[i].y = d.y;
            }
            else
            {
                pp0.ps[i].x =pp0.ps[i].y = pp0.ps[i].z = 0.0;
                //throw std::runtime_error("bcw_layerInfo: unknown error...");
            }

            /**********************************************************************************************************/
        }
    }

    if (ln != timeDimension().min())
    {
        n0.i = static_cast<int>(ln+0); n0.x = n0.i*timeDimension().step();
        n1.i = static_cast<int>(ln-1); n1.x = n1.i*timeDimension().step();
        ProblemParams &pp0 = sourceParams[ln];
        ProblemParams &pp1 = sourceParams[ln-1];

        for (size_t i=0; i<heatSourceNumber; i++)
        {
            const_backward.i = i+1;

            f0[0] = pp0.fv[i].x; f0[1] = pp0.fv[i].y; f0[2] = pp0.fd[i].x; f0[3] = pp0.fd[i].y;
            const_backward.next( f0, n0, f1, n1, method );
            pp1.fv[i] = SpacePoint(f1[0], f1[1]);
            pp1.fd[i] = SpacePoint(f1[2], f1[3]);

            /**********************************************************************************************************/

            const SpacePoint &zi = pp1.zv[i];
            if (isPointOnPlate(zi))
            {
                SpacePoint d;
                pp1.ps[i].z = DeltaFunction::lumpedPoint4(p, zi, spaceDimensionX(), spaceDimensionY(), d);
                pp1.ps[i].x = d.x; pp1.ps[i].y = d.y;
            }
            else
            {
                pp1.ps[i].x =pp1.ps[i].y = pp1.ps[i].z = 0.0;
               // throw std::runtime_error("bcw_layerInfo: unknown error...");
            }

            /**********************************************************************************************************/
        }
    }

    //unsigned int ln = static_cast<unsigned int>(tn.i);
    //ProblemParams &pp = sourceParams[ln];
    //std::string msg = "[ " + std::to_string(ln+1) + " ]";
    //msg += " [ " + std::to_string(pp.f[0].x) + ", " + std::to_string(pp.f[0].y) + " ]";
    //msg += " [ " + std::to_string(pp.f[1].x) + ", " + std::to_string(pp.f[1].y) + " ]";
    //msg += " [ " + std::to_string(p.min()) + ", " + std::to_string(p.max()) + " ]";
    //IPrinter::printSeperatorLine(msg.data());
    //IPrinter::printMatrix(u);
    //bcw_saveToImage(p, tn);
}

double Solver1::A1(const PointNodeODE &, size_t r, size_t c, size_t i) const
{
    if (i==1) { double data[2][2] = { { +0.24, +0.00 }, { +0.00, +0.35 } }; return data[r-1][c-1]; }
    if (i==2) { double data[2][2] = { { -0.12, +0.00 }, { +0.00, -0.19 } }; return data[r-1][c-1]; }
    return 0.0;
}

double Solver1::A2(const PointNodeODE &, size_t r, size_t c, size_t i) const
{
    if (i==1) { double data[2][2] = { { +0.00, +0.00 }, { +0.00, +0.00 } }; return data[r-1][c-1]; }
    if (i==2) { double data[2][2] = { { +0.00, +0.00 }, { +0.00, +0.00 } }; return data[r-1][c-1]; }
    return 0.0;
}

double Solver1::A3(const PointNodeODE &, size_t r, size_t i) const
{
    if (i==1) { double data[2] = { +0.00, +0.00 }; return data[r-1]; }
    if (i==2) { double data[2] = { +0.00, +0.00 }; return data[r-1]; }
    return 0.0;
}

double Solver1::A4(const PointNodeODE &, size_t r, size_t i) const
{
    if (i==1) { double data[2] = { +0.79, +0.60 }; return data[r-1]; }
    if (i==2) { double data[2] = { -0.84, +0.83 }; return data[r-1]; }
    return 0.0;
}

auto Solver1::frw_f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double
{
    size_t ln = static_cast<size_t>(tn.i);

    const ProblemParams &pp = sourceParams[ln];

    double fx = 0.0;//_lambda0*_environmentTemperature;
    //if (ln >= 400) return 0.0;

    //if (ln>=10) return 0.0;

    for (size_t i=0; i<heatSourceNumber; i++)
    {
        const SpacePoint &z = pp.zv[i];
        if ( isPointOnPlate(z) )
        {
            fx += pp.q[i] * DeltaFunction::gaussian(sn, z, SpacePoint(spaceDimensionX().step()*_factor, spaceDimensionY().step()*_factor));
        }
    }

    return fx;
}

void Solver1::frw_layerInfo(const DoubleMatrix &u, const TimeNodePDE &tn) const
{
//    if (saveMinMaxU)
//    {
//        if (u.min() < minU) const_cast<Solver1*>(this)->minU = u.min();
//        if (u.max() > maxU) const_cast<Solver1*>(this)->maxU = u.max();
//    }

    HeatEquationIBVP &const_forward = const_cast<HeatEquationIBVP&>(this->forward);
    int ln = static_cast<int>(tn.i);
    double ht = timeDimension().step();
    DoubleVector z0(4), z1(4);
    PointNodeODE n0, n1;

    if (ln == timeDimension().min())
    {
        n0.i = static_cast<int>(ln); n0.x = n0.i*ht;
        ProblemParams &pp0 = sourceParams[ln];

        for (size_t i=0; i<heatSourceNumber; i++)
        {
            const_forward.i = i+1;

            /**********************************************************************************************************/

            const_forward.start( z0, n0 );
            pp0.zv[i] = SpacePoint(z0[0], z0[1]);
            pp0.zd[i] = SpacePoint(z0[2], z0[3]);

            /**********************************************************************************************************/

            const SpacePoint &z = pp0.zv[i];
            if (isPointOnPlate(z))
            {
                double qi = 0.0;
                double vi = 0.0;
                for (size_t j=0; j<measrPointNumber; j++)
                {
                    double nominal1 = nomnU1[i][j];
                    double nominal2 = nomnU2[i][j];
                    SpacePoint d;
                    pp0.u[j].z = DeltaFunction::lumpedPoint4(u, measurePoints[j], spaceDimensionX(), spaceDimensionY(), d);
                    pp0.u[j].x = d.x; pp0.u[j].y = d.y;

                    const SpacePoint &mp = measurePoints[j];
                    double dist = sqrt((z.x - mp.x)*(z.x - mp.x) + (z.y - mp.y)*(z.y - mp.y));
                    qi += (alpha1[i][j]*dist*dist + alpha2[i][j]*dist + alpha3[i][j]) * (pp0.u[j].z - nominal1);
                    vi += (betta1[i][j]*dist*dist + betta2[i][j]*dist + betta3[i][j]) * (pp0.u[j].z - nominal2);
                }
                pp0.q[i] = qi;
                pp0.v[i] = vi;
            }
            else
            {
                pp0.q[i] = 0.0;
                pp0.v[i] = 0.0;
            }

            /**********************************************************************************************************/

            //printf("ln: %d i: %u z0: %10.6f %10.6f %10.6f %10.6f q:%10.6f v:%10.6f u:%10.6f %10.6f\n", ln, i, pp0.zv[i].x, pp0.zv[i].y, pp0.zd[i].x, pp0.zd[i].y, pp0.q[i], pp0.v[i], u.min(), u.max());
        }
    }

    if (ln != timeDimension().max())
    {
        n0.i = static_cast<int>(ln+0); n0.x = n0.i*ht;
        n1.i = static_cast<int>(ln+1); n1.x = n1.i*ht;
        ProblemParams &pp0 = sourceParams[ln];
        ProblemParams &pp1 = sourceParams[ln+1];

        for (size_t i=0; i<heatSourceNumber; i++)
        {
            const_forward.i = i+1;

            z0[0] = pp0.zv[i].x; z0[1] = pp0.zv[i].y; z0[2] = pp0.zd[i].x; z0[3] = pp0.zd[i].y;
            const_forward.next( z0, n0, z1, n1, method );
            pp1.zv[i] = SpacePoint(z1[0], z1[1]);
            pp1.zd[i] = SpacePoint(z1[2], z1[3]);

            /**********************************************************************************************************/

            const SpacePoint &zi = pp1.zv[i];
            if (isPointOnPlate(zi))
            {
                double qi = 0.0;
                double vi = 0.0;
                for (size_t j=0; j<measrPointNumber; j++)
                {
                    double nominal1 = nomnU1[i][j];
                    double nominal2 = nomnU2[i][j];

                    const SpacePoint &mp = measurePoints[j];
                    double dist = sqrt((zi.x - mp.x)*(zi.x - mp.x) + (zi.y - mp.y)*(zi.y - mp.y));
                    qi += (alpha1[i][j]*dist*dist + alpha2[i][j]*dist + alpha3[i][j]) * (pp0.u[j].z - nominal1);
                    vi += (betta1[i][j]*dist*dist + betta2[i][j]*dist + betta3[i][j]) * (pp0.u[j].z - nominal2);

                    SpacePoint d;
                    pp1.u[j].z = DeltaFunction::lumpedPoint4(u, measurePoints[j], spaceDimensionX(), spaceDimensionY(), d);
                    pp1.u[j].x = d.x; pp0.u[j].y = d.y;
                }
                pp1.q[i] = qi;
                pp1.v[i] = vi;
            }
            else
            {
                pp1.q[i] = 0.0;
                pp1.v[i] = 0.0;
            }

            //printf("ln: %d i: %u z0: %10.6f %10.6f %10.6f %10.6f q:%10.6f v:%10.6f u:%10.6f %10.6f\n", ln, i, pp0.zv[i].x, pp0.zv[i].y, pp0.zd[i].x, pp0.zd[i].y, pp0.q[i], pp0.v[i], u.min(), u.max());

            /**********************************************************************************************************/
        }
    }

    if (ln == timeDimension().max())
    {
        const_cast<Solver1*>(this)->U = u;
    }

    //unsigned int ln = static_cast<unsigned int>(tn.i);
    //ProblemParams &pp = sourceParams[ln];
    //std::string msg = "[ " + std::to_string(ln+1) + " ]";
    //msg += " [ " + std::to_string(pp.z[0].x) + ", " + std::to_string(pp.z[0].y) + " ]";
    //msg += " [ " + std::to_string(pp.z[1].x) + ", " + std::to_string(pp.z[1].y) + " ]";
    //msg += " [ " + std::to_string(pp.v[0]) + ", " + std::to_string(pp.v[1]) + " ]";
    //msg += " [ " + std::to_string(pp.q[0]) + ", " + std::to_string(pp.q[1]) + " ]";
    //msg += " [ " + std::to_string(u.min()) + ", " + std::to_string(u.max()) + " ]";
    //IPrinter::printSeperatorLine(msg.data());
    //IPrinter::printMatrix(u);
    if (drawImage) frw_saveToImage(u, tn);
}

void Solver1::gradient(const DoubleVector &x, DoubleVector &g) const
{
    const_cast<Solver1*>(this)->vectorToParameter(x);
    forward.implicit_calculate_D2V1();
    backward.implicit_calculate_D2V1();

    int min = timeDimension().min();
    int max = timeDimension().max();
    double ht = timeDimension().step();

    size_t size = heatSourceNumber*measrPointNumber;
    g.clear();
    g.resize(x.length(), 0.0);

    for (size_t i=0; i<heatSourceNumber; i++)
    {
        for (size_t j=0; j<measrPointNumber; j++)
        {
            const SpacePoint &mp = measurePoints[j];

            g[0*size + i*measrPointNumber  + j] = 0.0;
            g[1*size + i*measrPointNumber  + j] = 0.0;
            g[2*size + i*measrPointNumber  + j] = 0.0;
            g[3*size + i*measrPointNumber  + j] = 0.0;
            g[4*size + i*measrPointNumber  + j] = 0.0;
            g[5*size + i*measrPointNumber  + j] = 0.0;
            g[6*size + i*measrPointNumber  + j] = 0.0;
            g[7*size + i*measrPointNumber  + j] = 0.0;

            int ln = min;
            PointNodeODE node; node.i = ln; node.x = node.i*ht;

            const SpacePoint &zi0 = sourceParams[ln].zv[i];
            if (isPointOnPlate(zi0))
            {
                double dist = sqrt((zi0.x - mp.x)*(zi0.x - mp.x) + (zi0.y - mp.y)*(zi0.y - mp.y));
                double dif1 = sourceParams[ln].u[j].z-nomnU1[i][j];
                double dif2 = sourceParams[ln].u[j].z-nomnU2[i][j];
                double val0 = sourceParams[ln].ps[i].z;
                double val1 = A4(node, 1, i+1)*sourceParams[ln].fv[i].x + A4(node, 2, i+1)*sourceParams[ln].fv[i].y;

                g[0*size + i*measrPointNumber  + j] += 0.5*val0*dif1*dist*dist;
                g[1*size + i*measrPointNumber  + j] += 0.5*val0*dif1*dist;
                g[2*size + i*measrPointNumber  + j] += 0.5*val0*dif1;
                g[3*size + i*measrPointNumber  + j] += 0.5*val1*dif2*dist*dist;
                g[4*size + i*measrPointNumber  + j] += 0.5*val1*dif2*dist;
                g[5*size + i*measrPointNumber  + j] += 0.5*val1*dif2;
                g[6*size + i*measrPointNumber  + j] += 0.5*val0*(alpha1[i][j]*(dist*dist)+alpha2[i][j]*dist+alpha3[i][j]);
                g[7*size + i*measrPointNumber  + j] += 0.5*val1*(betta1[i][j]*(dist*dist)+betta2[i][j]*dist+betta3[i][j]);
            }
            else
            {
            }

            for (int ln=min+1; ln<max; ln++)
            {
                node.i = ln; node.x = node.i*ht;

                const SpacePoint &zi = sourceParams[ln].zv[i];
                if (isPointOnPlate(zi))
                {
                    double dist = sqrt((zi.x - mp.x)*(zi.x - mp.x) + (zi.y - mp.y)*(zi.y - mp.y));
                    double dif1 = sourceParams[ln].u[j].z-nomnU1[i][j];
                    double dif2 = sourceParams[ln].u[j].z-nomnU2[i][j];
                    double val0 = sourceParams[ln].ps[i].z;
                    double val1 = A4(node, 1, i+1)*sourceParams[ln].fv[i].x + A4(node, 2, i+1)*sourceParams[ln].fv[i].y;

                    g[0*size + i*measrPointNumber  + j] += val0*dif1*dist*dist;
                    g[1*size + i*measrPointNumber  + j] += val0*dif1*dist;
                    g[2*size + i*measrPointNumber  + j] += val0*dif1;
                    g[3*size + i*measrPointNumber  + j] += val1*dif2*dist*dist;
                    g[4*size + i*measrPointNumber  + j] += val1*dif2*dist;
                    g[5*size + i*measrPointNumber  + j] += val1*dif2;
                    g[6*size + i*measrPointNumber  + j] += val0*(alpha1[i][j]*(dist*dist)+alpha2[i][j]*dist+alpha3[i][j]);
                    g[7*size + i*measrPointNumber  + j] += val1*(betta1[i][j]*(dist*dist)+betta2[i][j]*dist+betta3[i][j]);
                }
            }

            ln = max;
            node.i = ln; node.x = node.i*ht;

            const SpacePoint &zi1 = sourceParams[ln].zv[i];
            if (isPointOnPlate(zi1))
            {
                double dist = sqrt((zi1.x - mp.x)*(zi1.x - mp.x) + (zi1.y - mp.y)*(zi1.y - mp.y));
                double dif1 = sourceParams[ln].u[j].z-nomnU1[i][j];
                double dif2 = sourceParams[ln].u[j].z-nomnU2[i][j];
                double val0 = sourceParams[ln].ps[i].z;
                double val1 = A4(node, 1, i+1)*sourceParams[ln].fv[i].x + A4(node, 2, i+1)*sourceParams[ln].fv[i].y;

                g[0*size + i*measrPointNumber  + j] += 0.5*val0*dif1*dist*dist;
                g[1*size + i*measrPointNumber  + j] += 0.5*val0*dif1*dist;
                g[2*size + i*measrPointNumber  + j] += 0.5*val0*dif1;
                g[3*size + i*measrPointNumber  + j] += 0.5*val1*dif2*dist*dist;
                g[4*size + i*measrPointNumber  + j] += 0.5*val1*dif2*dist;
                g[5*size + i*measrPointNumber  + j] += 0.5*val1*dif2;
                g[6*size + i*measrPointNumber  + j] += 0.5*val0*(alpha1[i][j]*(dist*dist)+alpha2[i][j]*dist+alpha3[i][j]);
                g[7*size + i*measrPointNumber  + j] += 0.5*val1*(betta1[i][j]*(dist*dist)+betta2[i][j]*dist+betta3[i][j]);
            }

            g[0*size + i*measrPointNumber  + j] *= -ht;
            g[1*size + i*measrPointNumber  + j] *= -ht;
            g[2*size + i*measrPointNumber  + j] *= -ht;

            g[3*size + i*measrPointNumber  + j] *= -ht;
            g[4*size + i*measrPointNumber  + j] *= -ht;
            g[5*size + i*measrPointNumber  + j] *= -ht;

            g[6*size + i*measrPointNumber  + j] *= +ht;
            g[7*size + i*measrPointNumber  + j] *= +ht;
        }
    }
}

double Solver1::fx(const DoubleVector &x) const
{
    const_cast<Solver1*>(this)->vectorToParameter(x);
    forward.implicit_calculate_D2V1();
    return integral(U);
}

auto Solver1::integral(const DoubleMatrix &) const -> double
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

auto Solver1::print(unsigned int i, const DoubleVector &x, const DoubleVector &g,
                             double f, double alpha, GradientMethod::MethodResult result) const -> void
{
    printf_s("%4d %16.8f %16.8f\n", i, f, alpha);
    if (alpha > 0.0)
    const_cast<Solver1*>(this)->gradMethod->setR1MinimizeEpsilon(alpha, alpha/100.0);
    //if (fabs(alpha) < 0.0001)
    //const_cast<Solver1*>(this)->gradMethod->setR1MinimizeEpsilon(0.1, 0.001);
}

auto Solver1::project(DoubleVector &x) const -> void
{

}

auto Solver1::project(DoubleVector &x, unsigned int index) -> void
{}


//--------------------------------------------------------------------------------------------------------------//

HeatEquationIBVP::HeatEquationIBVP(Solver1 *solver) { this->solver = solver; }

HeatEquationIBVP::HeatEquationIBVP(const HeatEquationIBVP &) {}

HeatEquationIBVP & HeatEquationIBVP::operator =(const HeatEquationIBVP &) { return *this; }

HeatEquationIBVP::~HeatEquationIBVP() {}

auto HeatEquationIBVP::initial(const SpaceNodePDE &sn, InitialCondition cn) const -> double { return solver->frw_initial(sn, cn); }

auto HeatEquationIBVP::boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &cn) const -> double { return solver->frw_boundary(sn, tn, cn); }

auto HeatEquationIBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const-> double { return solver->frw_f(sn, tn); }

auto HeatEquationIBVP::layerInfo(const DoubleMatrix &u, const TimeNodePDE &tn) const -> void { solver->frw_layerInfo(u, tn); }

auto HeatEquationIBVP::timeDimension() const -> Dimension { return solver->timeDimension(); }

auto HeatEquationIBVP::spaceDimensionX() const -> Dimension { return solver->spaceDimensionX(); }

auto HeatEquationIBVP::spaceDimensionY() const -> Dimension { return solver->spaceDimensionY(); }

auto HeatEquationIBVP::spaceDimensionZ() const -> Dimension { throw std::runtime_error("Not supportted..."); }

//--------------------------------------------------------------------------------------------------------------//

auto HeatEquationIBVP::A(const PointNodeODE &node, size_t r, size_t c) const -> double { return solver->A1(node, r, c, i); }

auto HeatEquationIBVP::B(const PointNodeODE &node, size_t r, size_t c) const -> double { return solver->A2(node, r, c, i); }

auto HeatEquationIBVP::C(const PointNodeODE &node, size_t r) const -> double
{
    size_t ln = static_cast<size_t>(node.i);
    return solver->A3(node, r, i) + solver->A4(node, r, i) * solver->sourceParams[ln].v[i-1];
}

auto HeatEquationIBVP::initial(InitialCondition c, size_t r) const -> double
{
    if (c == InitialCondition::InitialValue)
    {
        //if (i==1) { const double data[2] = { 0.04, 0.04 }; return data[r-1]; }
        //if (i==2) { const double data[2] = { 0.96, 0.04 }; return data[r-1]; }
        if (i==1) { const double data[2] = { 0.10, 0.10 }; return data[r-1]; }
        if (i==2) { const double data[2] = { 0.90, 0.10 }; return data[r-1]; }
    }
    if (c == InitialCondition::InitialFirstDerivative)
    {
//        if (i==1) { const double data[2] = { +0.35, +0.56 }; return data[r-1]; }
//        if (i==2) { const double data[2] = { -0.54, +0.42 }; return data[r-1]; }
        if (i==1) { const double data[2] = { +0.00, +0.00 }; return data[r-1]; }
        if (i==2) { const double data[2] = { -0.00, +0.00 }; return data[r-1]; }
    }
    throw std::runtime_error(std::string("auto HeatEquationIBVP::initial(InitialCondition, size_t row) const -> double"));
}

auto HeatEquationIBVP::count() const -> size_t { return 2; }

auto HeatEquationIBVP::dimension() const -> Dimension { return solver->timeDimension(); }

auto HeatEquationIBVP::iterationInfo(const DoubleVector &, const PointNodeODE &) const -> void {}

//--------------------------------------------------------------------------------------------------------------//

//--------------------------------------------------------------------------------------------------------------//

HeatEquationFBVP::HeatEquationFBVP(Solver1 *solver) { this->solver = solver; }

HeatEquationFBVP::HeatEquationFBVP(const HeatEquationFBVP &) {}

HeatEquationFBVP & HeatEquationFBVP::operator =(const HeatEquationFBVP &) { return *this; }

HeatEquationFBVP::~HeatEquationFBVP() {}

auto HeatEquationFBVP::final(const SpaceNodePDE &sn, FinalCondition cn) const -> double { return solver->bcw_final(sn, cn); }

auto HeatEquationFBVP::boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &cn) const -> double { return solver->bcw_boundary(sn, tn, cn); }

auto HeatEquationFBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const-> double { return solver->bcw_f(sn, tn); }

auto HeatEquationFBVP::layerInfo(const DoubleMatrix &p, const TimeNodePDE &tn) const -> void { solver->bcw_layerInfo(p, tn); }

auto HeatEquationFBVP::timeDimension() const -> Dimension { return solver->timeDimension(); }

auto HeatEquationFBVP::spaceDimensionX() const -> Dimension { return solver->spaceDimensionX(); }

auto HeatEquationFBVP::spaceDimensionY() const -> Dimension { return solver->spaceDimensionY(); }

auto HeatEquationFBVP::spaceDimensionZ() const -> Dimension { throw std::runtime_error(""); }

//--------------------------------------------------------------------------------------------------------------//

auto HeatEquationFBVP::A(const PointNodeODE &node, size_t r, size_t c) const -> double { return -solver->A1(node, r, c, i); }

auto HeatEquationFBVP::B(const PointNodeODE &node, size_t r, size_t c) const -> double { return +solver->A2(node, r, c, i); }

auto HeatEquationFBVP::C(const PointNodeODE &node, size_t r) const -> double
{
    size_t ln = static_cast<size_t>(node.i);
    const ProblemParams &pp = solver->sourceParams[ln];
    double sum = 0.0;

    const SpacePoint &zi = pp.zv[i-1];
    if (solver->isPointOnPlate(zi))
    {
        double val0 = pp.ps[i-1].z;
        double val1 = solver->A4(node, 1, i)*pp.fv[i-1].x + solver->A4(node, 2, i)*pp.fv[i-1].y;

        for (size_t j=0; j<solver->measrPointNumber; j++)
        {
            const SpacePoint &mp = solver->measurePoints[j];
            double dist = sqrt((zi.x - mp.x)*(zi.x - mp.x) + (zi.y - mp.y)*(zi.y - mp.y));

            double diff1 = pp.u[j].z-solver->nomnU1[i-1][j];
            double diff2 = pp.u[j].z-solver->nomnU2[i-1][j];

            sum += val0 * (2.0*solver->alpha1[i-1][j]*dist + solver->alpha2[i-1][j]) * diff1;
            sum += val1 * (2.0*solver->betta1[i-1][j]*dist + solver->betta2[i-1][j]) * diff2;
        }

        sum += ((r==1) ? pp.ps[i-1].x * pp.q[i-1] : pp.ps[i-1].y * pp.q[i-1]);
    }
    else
    {

    }

    return sum;
}

auto HeatEquationFBVP::final(FinalCondition c, size_t r) const -> double
{
    if (c == FinalCondition::FinalValue)
    {
        if (i==1) { const double data[2] = { 0.00, 0.00 }; return data[r-1]; }
        if (i==2) { const double data[2] = { 0.00, 0.00 }; return data[r-1]; }
    }
    if (c == FinalCondition::FinalFirstDerivative)
    {
        if (i==1) { return 0.0; }
        if (i==2) { return 0.0; }
    }
    return 0.0;
}

auto HeatEquationFBVP::count() const -> size_t { return 2; }

auto HeatEquationFBVP::dimension() const -> Dimension { return solver->timeDimension(); }

//auto HeatEquationFBVP::iterationInfo(const DoubleVector &z, const PointNodeODE &node) const -> void
//{
//}

void Solver1::vectorToParameter(const DoubleVector &x)
{
    size_t size = heatSourceNumber*measrPointNumber;
    for (size_t i=0; i<heatSourceNumber; i++)
    {
        for (size_t j=0; j<measrPointNumber; j++)
        {
            alpha1[i][j] = x[0*size + i*measrPointNumber + j];
            alpha2[i][j] = x[1*size + i*measrPointNumber + j];
            alpha3[i][j] = x[2*size + i*measrPointNumber + j];
            betta1[i][j] = x[3*size + i*measrPointNumber + j];
            betta2[i][j] = x[4*size + i*measrPointNumber + j];
            betta3[i][j] = x[5*size + i*measrPointNumber + j];
            nomnU1[i][j] = x[6*size + i*measrPointNumber + j];
            nomnU2[i][j] = x[7*size + i*measrPointNumber + j];
        }
    }
}

void Solver1::parameterToVector(DoubleVector &x)
{
    size_t size = heatSourceNumber*measrPointNumber;
    x.clear();
    x.resize(size*8);

    for (size_t i=0; i<heatSourceNumber; i++)
    {
        for (size_t j=0; j<measrPointNumber; j++)
        {
            x[0*size + i*measrPointNumber + j] = alpha1[i][j];
            x[1*size + i*measrPointNumber + j] = alpha2[i][j];
            x[2*size + i*measrPointNumber + j] = alpha3[i][j];
            x[3*size + i*measrPointNumber + j] = betta1[i][j];
            x[4*size + i*measrPointNumber + j] = betta2[i][j];
            x[5*size + i*measrPointNumber + j] = betta3[i][j];
            x[6*size + i*measrPointNumber + j] = nomnU1[i][j];
            x[7*size + i*measrPointNumber + j] = nomnU2[i][j];
        }
    }
}


