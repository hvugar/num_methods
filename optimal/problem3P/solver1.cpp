#include "solver1.h"

using namespace p3p1;

void Solver1::Main(int argc, char **argv)
{
    QGuiApplication app(argc, argv);

    Solver1 s(Dimension(0.002, 0, 500), Dimension(0.005, 0, 200), Dimension(0.005, 0, 200));
    s.setPointNumber(2, 4);
    s.forward.solver = &s;
    s.backward.solver = &s;
    s.forward.setThermalDiffusivity(0.01);
    s.forward.setThermalConductivity(0.0);
    s.forward.setThermalConvection(-s.lambda0);
    s.forward.implicit_calculate_D2V1();

    s.backward.setThermalDiffusivity(-0.01);
    s.backward.setThermalConductivity(0.0);
    s.backward.setThermalConvection(s.lambda0);
    s.backward.implicit_calculate_D2V1();
}

Solver1::Solver1(const Dimension &timeDimension, const Dimension &spaceDimensionX, const Dimension &spaceDimensionY)
{
    this->_timeDimension = timeDimension;
    this->_spaceDimensionX = spaceDimensionX;
    this->_spaceDimensionY = spaceDimensionY;
}

Solver1::~Solver1() {}

void Solver1::setPointNumber(size_t heatSourceNumber, size_t measrPointNumber)
{
    this->heatSourceNumber = heatSourceNumber;
    this->measrPointNumber = measrPointNumber;

    alpha1.resize(heatSourceNumber, measrPointNumber, 0.200);
    alpha2.resize(heatSourceNumber, measrPointNumber, 0.250);
    alpha3.resize(heatSourceNumber, measrPointNumber, 0.275);
    betta1.resize(heatSourceNumber, measrPointNumber, 0.200);
    betta2.resize(heatSourceNumber, measrPointNumber, 0.250);
    betta3.resize(heatSourceNumber, measrPointNumber, 0.275);

    this->measurePoints = new SpacePoint[measrPointNumber];
    //this->measurePointValues = new SpacePoint[measrPointNumber];
    this->measurePoints[0] = SpacePoint(0.50, 0.20);
    this->measurePoints[1] = SpacePoint(0.20, 0.50);
    this->measurePoints[2] = SpacePoint(0.50, 0.80);
    this->measurePoints[3] = SpacePoint(0.80, 0.50);

    nU.resize(heatSourceNumber, measrPointNumber, 0.0);

    sourceParams = new ProblemParams[timeDimension().size()];
    for (size_t ln=0; ln<timeDimension().size(); ln++)
    {
        ProblemParams & pp = sourceParams[ln];
        pp.q = new double[heatSourceNumber];
        pp.v = new double[heatSourceNumber];
        pp.z = new SpacePoint[heatSourceNumber];
        pp.zt = new SpacePoint[heatSourceNumber];
        pp.f = new SpacePoint[heatSourceNumber];
        //pp.ft = new SpacePoint[heatSourceNumber];

        pp.u = new SpacePoint[measrPointNumber];
    }
}

double Solver1::frw_initial(const SpaceNodePDE &, InitialCondition) const { return frw_initialValue; }

double Solver1::frw_boundary(const SpaceNodePDE &, const TimeNodePDE &, BoundaryConditionPDE &cn) const
{
    //cn = BoundaryConditionPDE::Robin(lambda, -1.0, lambda);
    cn = BoundaryConditionPDE::Neumann(1.0, 0.0);
    return environmentTemperature;
}

void Solver1::frw_saveToImage(const DoubleMatrix &u UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM) const
{
#ifdef USE_LIB_IMAGING
    static double MIN = +10000.0;
    static double MAX = -10000.0;
    if (u.max()>MAX) MAX = u.max();
    if (u.min()<MIN) MIN = u.min();
//    printf(">>> %6d %14.6f %14.6f %14.6f %14.6f\n", tn.i, u.min(), u.max(), MIN, MAX);

    //if (tn.i%10==0)
    {
        QString filename = QString("data/problem3P/f/png/%1.png").arg(tn.i, 8, 10, QChar('0'));
        QPixmap pixmap;
        //visualizeMatrixHeat(u, 0.0, 46.052, pixmap, 101, 101);
        visualizeMatrixHeat(u, 1.0, 78.0, pixmap, spaceDimensionX().size(), spaceDimensionY().size());
        pixmap.save(filename);
    }
#endif
}

double Solver1::bcw_final(const SpaceNodePDE &sn, FinalCondition condition) const
{
    size_t i = static_cast<size_t>(sn.i);
    size_t j = static_cast<size_t>(sn.j);
    return -2.0*(U[j][i]-V[j][i]);
}

double Solver1::bcw_boundary(const SpaceNodePDE &, const TimeNodePDE &, BoundaryConditionPDE &cn) const
{
    cn = BoundaryConditionPDE(BoundaryCondition::Neumann, 1.0, 0.0, 0.0);
    return 0.0;
}

auto Solver1::bcw_f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double
{
    return 0.0;
}

void Solver1::bcw_layerInfo(const DoubleMatrix &u, const TimeNodePDE &tn) const
{
    //    HeatEquationFBVP *const_this = const_cast<HeatEquationFBVP*>(this);
    //    int ln = static_cast<int>(tn.i);

    //    for (size_t i=0; i<heatSourceNumber; i++)
    //    {
    //        solver->measurePointValues[j].z = DeltaFunction::lumpedPoint4(u, solver->measurePoints[j], spaceDimensionX(), spaceDimensionY());
    //    }
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


//--------------------------------------------------------------------------------------------------------------//

HeatEquationIBVP::HeatEquationIBVP(Solver1 *solver) { this->solver = solver; }

HeatEquationIBVP::HeatEquationIBVP(const HeatEquationIBVP &) {}

HeatEquationIBVP & HeatEquationIBVP::operator =(const HeatEquationIBVP &) { return *this; }

HeatEquationIBVP::~HeatEquationIBVP() {}

auto HeatEquationIBVP::initial(const SpaceNodePDE &sn, InitialCondition cn) const -> double { return solver->frw_initial(sn, cn); }

auto HeatEquationIBVP::boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &cn) const -> double { return solver->frw_boundary(sn, tn, cn); }

auto HeatEquationIBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const-> double { return solver->frw_f(sn, tn); }

auto HeatEquationIBVP::layerInfo(const DoubleMatrix &u, const TimeNodePDE &tn) const -> void
{
    HeatEquationIBVP *const_this = const_cast<HeatEquationIBVP*>(this);
    int ln = static_cast<int>(tn.i);
    DoubleVector z0(4), z1(4);
    PointNodeODE n0, n1;

    n0.i = static_cast<int>(ln+0); n0.x = n0.i*timeDimension().step();
    n1.i = static_cast<int>(ln+1); n1.x = n1.i*timeDimension().step();

    if (ln == timeDimension().min())
    {
        for (size_t i=0; i<solver->heatSourceNumber; i++)
        {
            const_this->i = i+1;
            double qi = 0.0;
            double vi = 0.0;
            ProblemParams &param = solver->sourceParams[ln];

            const_this->start( z0, n0 );
            param.z[i] = SpacePoint(z0[0], z0[1]);
            param.zt[i] = SpacePoint(z0[2], z0[3]);
            const SpacePoint &zi = param.z[i];

            for (size_t j=0; j<solver->measrPointNumber; j++)
            {
                double nominal = solver->nU.at(i, j);
                param.u[j].z = DeltaFunction::lumpedPoint4(u, solver->measurePoints[j], spaceDimensionX(), spaceDimensionY());
                param.u[j].x = param.u[j].y = 0.0;

                const SpacePoint &mp = solver->measurePoints[j];
                double dist = sqrt((zi.x - mp.x)*(zi.x - mp.x) + (zi.y - mp.y)*(zi.y - mp.y));
                qi += (solver->alpha1[i][j]*dist*dist + solver->alpha2[i][j]*dist + solver->alpha3[i][j]) * (param.u[j].z - nominal);
                vi += (solver->betta1[i][j]*dist*dist + solver->betta2[i][j]*dist + solver->betta3[i][j]) * (param.u[j].z - nominal);
            }

            param.q[i] = qi;
            param.v[i] = vi;
        }
    }

    if (ln != timeDimension().max())
    {
        for (size_t i=0; i<solver->heatSourceNumber; i++)
        {
            const_this->i = i+1;
            double qi = 0.0;
            double vi = 0.0;
            ProblemParams &pp0 = solver->sourceParams[ln];
            ProblemParams &pp1 = solver->sourceParams[ln+1];

            z0[0] = pp0.z[i].x; z0[1] = pp0.z[i].y; z0[2] = pp0.zt[i].x; z0[3] = pp0.zt[i].y;
            const_this->next( z0, n0, z1, n1, ODESolverMethod::EULER );
            pp1.z[i]  = SpacePoint(z1[0], z1[1]);
            pp1.zt[i] = SpacePoint(z1[2], z1[3]);
            const SpacePoint &zi = pp1.z[i];

            for (size_t j=0; j<solver->measrPointNumber; j++)
            {
                double nominal = solver->nU.at(i, j);
                pp1.u[j].z = DeltaFunction::lumpedPoint3(u, solver->measurePoints[j], spaceDimensionX(), spaceDimensionY());
                pp1.u[j].x = pp0.u[j].y = 0.0;

                const SpacePoint &mp = solver->measurePoints[j];
                double dist = sqrt((zi.x - mp.x)*(zi.x - mp.x) + (zi.y - mp.y)*(zi.y - mp.y));
                qi += (solver->alpha1[i][j]*dist*dist + solver->alpha2[i][j]*dist + solver->alpha3[i][j]) * (pp1.u[j].z - nominal);
                vi += (solver->betta1[i][j]*dist*dist + solver->betta2[i][j]*dist + solver->betta3[i][j]) * (pp1.u[j].z - nominal);
            }
            pp1.q[i] = qi;
            pp1.v[i] = vi;
        }
    }

    solver->frw_layerInfo(u, tn);
}

auto Solver1::frw_f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double
{
    double fx = lambda0*environmentTemperature;
    size_t ln = static_cast<size_t>(tn.i);

    for (size_t i=0; i<heatSourceNumber; i++)
    {
        const SpacePoint &z = sourceParams[ln].z[i];
        if (0.0 <= z.x <= 1.0 && 0.0 <= z.y <= 1.0)
        {
            fx += sourceParams[ln].q[i] * DeltaFunction::gaussian(sn, z, SpacePoint(spaceDimensionX().step(), spaceDimensionY().step()));
        }
    }
    return fx ;
}

void Solver1::frw_layerInfo(const DoubleMatrix &u, const TimeNodePDE &tn) const
{
    unsigned int ln = static_cast<unsigned int>(tn.i);
    ProblemParams &pp = sourceParams[ln];
    std::string msg = "[ " + std::to_string(ln+1) + " ]";
    msg += " [ " + std::to_string(pp.z[0].x) + ", " + std::to_string(pp.z[0].y) + " ]";
    msg += " [ " + std::to_string(pp.z[1].x) + ", " + std::to_string(pp.z[1].y) + " ]";
    msg += " [ " + std::to_string(pp.v[0]) + ", " + std::to_string(pp.v[1]) + " ]";
    msg += " [ " + std::to_string(pp.q[0]) + ", " + std::to_string(pp.q[1]) + " ]";
    msg += " [ " + std::to_string(u.min()) + ", " + std::to_string(u.max()) + " ]";

    IPrinter::printSeperatorLine(msg.data());
    //    IPrinter::printMatrix(u);
    frw_saveToImage(u, tn);
}

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
    return solver->A3(node, r, i)
            + solver->A4(node, r, i) * solver->sourceParams[ln].v[i-1];
}

auto HeatEquationIBVP::initial(InitialCondition c, size_t r) const -> double
{
    if (c == InitialCondition::InitialValue)
    {
        if (i==1) { const double data[2] = { 0.10, 0.10 }; return data[r-1]; }
        if (i==2) { const double data[2] = { 0.90, 0.10 }; return data[r-1]; }
    }
    if (c == InitialCondition::InitialFirstDerivative)
    {
        if (i==1) { const double data[2] = { 0.00, 0.00 }; return data[r-1]; }
        if (i==2) { const double data[2] = { 0.00, 0.00 }; return data[r-1]; }
    }
    throw std::runtime_error(std::string("auto HeatEquationIBVP::initial(InitialCondition, size_t row) const -> double"));
}

auto HeatEquationIBVP::count() const -> size_t { return 2; }

auto HeatEquationIBVP::dimension() const -> Dimension { return solver->timeDimension(); }

auto HeatEquationIBVP::iterationInfo(const DoubleVector &z, const PointNodeODE &node) const -> void
{}

//--------------------------------------------------------------------------------------------------------------//

//--------------------------------------------------------------------------------------------------------------//

HeatEquationFBVP::HeatEquationFBVP(Solver1 *solver) { this->solver = solver; }

HeatEquationFBVP::HeatEquationFBVP(const HeatEquationFBVP &) {}

HeatEquationFBVP & HeatEquationFBVP::operator =(const HeatEquationFBVP &) { return *this; }

HeatEquationFBVP::~HeatEquationFBVP() {}

auto HeatEquationFBVP::final(const SpaceNodePDE &sn, FinalCondition cn) const -> double { return solver->bcw_final(sn, cn); }

auto HeatEquationFBVP::boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &cn) const -> double { return solver->bcw_boundary(sn, tn, cn); }

auto HeatEquationFBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const-> double { return solver->bcw_f(sn, tn); }

auto HeatEquationFBVP::layerInfo(const DoubleMatrix &p, const TimeNodePDE &tn) const -> void
{
    HeatEquationFBVP *const_this = const_cast<HeatEquationFBVP*>(this);
    int ln = static_cast<int>(tn.i);
    DoubleVector f0(4), f1(4);
    PointNodeODE n0, n1;

    n0.i = static_cast<int>(ln+0); n0.x = n0.i*timeDimension().step();
    n1.i = static_cast<int>(ln-1); n1.x = n1.i*timeDimension().step();

    if (ln == timeDimension().max())
    {
        for (size_t j=0; j<solver->measrPointNumber; j++)
        {

        }
    }

    if (ln != timeDimension().min())
    {}


    solver->bcw_layerInfo(p, tn);
}

auto HeatEquationFBVP::timeDimension() const -> Dimension { return solver->timeDimension(); }

auto HeatEquationFBVP::spaceDimensionX() const -> Dimension { return solver->spaceDimensionX(); }

auto HeatEquationFBVP::spaceDimensionY() const -> Dimension { return solver->spaceDimensionY(); }

auto HeatEquationFBVP::spaceDimensionZ() const -> Dimension { throw std::runtime_error(""); }

//--------------------------------------------------------------------------------------------------------------//

auto HeatEquationFBVP::A(const PointNodeODE &node, size_t r, size_t c) const -> double { return -solver->A1(node, r, c, i); }

auto HeatEquationFBVP::B(const PointNodeODE &node, size_t r, size_t c) const -> double { return +solver->A2(node, r, c, i); }

auto HeatEquationFBVP::C(const PointNodeODE &node, size_t r) const -> double
{
    double res = 0.0;
    size_t ln = static_cast<size_t>(node.i);
    const SpacePoint &z = solver->sourceParams[ln].z[i-1];


    for (size_t j=1; j<solver->measrPointNumber; j++)
    {
    }

    if (i==1)
    {

    }

    if (i==2)
    {

    }


    return 0.0;
}

auto HeatEquationFBVP::final(FinalCondition c, size_t r) const -> double
{
    if (c == FinalCondition::FinalValue)
    {
        if (i==1) { const double data[2] = { 0.10, 0.90 }; return data[r-1]; }
        if (i==2) { const double data[2] = { 0.10, 0.10 }; return data[r-1]; }
    }
    if (c == FinalCondition::FinalFirstDerivative) return 0.0;
    throw std::runtime_error(std::string("auto HeatEquationIBVP::initial(InitialCondition, size_t row) const -> double"));
}

auto HeatEquationFBVP::count() const -> size_t { return 2; }

auto HeatEquationFBVP::dimension() const -> Dimension { return solver->timeDimension(); }

auto HeatEquationFBVP::iterationInfo(const DoubleVector &z, const PointNodeODE &node) const -> void
{
}


