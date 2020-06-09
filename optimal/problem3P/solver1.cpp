#include "solver1.h"

using namespace p3p1;

void Solver1::Main(int argc, char **argv)
{
    QGuiApplication app(argc, argv);
    Solver1 s(Dimension(0.002, 0, 500), Dimension(0.005, 0, 200), Dimension(0.005, 0, 200));
    s.setPointNumber(2, 4);

    HeatEquationIBVP forward(&s);
    forward.setThermalDiffusivity(0.01);
    forward.setThermalConvection(0.0);
    forward.setThermalConductivity(0.0);
    forward.implicit_calculate_D2V1();
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
    this->measurePointValues = new SpacePoint[measrPointNumber];
    this->measurePoints[0] = SpacePoint(0.50, 0.20);
    this->measurePoints[1] = SpacePoint(0.20, 0.50);
    this->measurePoints[2] = SpacePoint(0.50, 0.80);
    this->measurePoints[3] = SpacePoint(0.80, 0.50);

    A = new DoubleMatrix[heatSourceNumber];
    A[0].resize(2, 2); A[0][0][0] = +0.24; A[0][0][1] = 0.0; A[0][1][0] = 0.0; A[0][1][1] = +0.35;
    A[1].resize(2, 2); A[1][0][0] = -0.12; A[1][0][1] = 0.0; A[1][1][0] = 0.0; A[1][1][1] = -0.19;

    B = new DoubleMatrix[heatSourceNumber];
    B[0].resize(2, 2); B[0][0][0] = +0.00; B[0][0][1] = 0.0; B[0][1][0] = 0.0; B[0][1][1] = +0.00;
    B[1].resize(2, 2); B[1][0][0] = +0.00; B[1][0][1] = 0.0; B[1][1][0] = 0.0; B[1][1][1] = +0.00;

    C = new DoubleVector[heatSourceNumber];
    C[0].resize(2); C[0][0] = +0.79; C[0][1] = +0.60;
    C[1].resize(2); C[1][0] = -0.84; C[1][1] = +0.83;

    DoubleMatrix nominU(heatSourceNumber, measrPointNumber, 0.0);

    sourceParams = new HeatSourceParam*[heatSourceNumber];
    sourceParams[0] = new HeatSourceParam[timeDimension().size()];
    sourceParams[1] = new HeatSourceParam[timeDimension().size()];
}

double Solver1::frw_initial(const SpaceNodePDE &, InitialCondition) const { return frw_initialValue; }

double Solver1::frw_boundary(const SpaceNodePDE &, const TimeNodePDE &, BoundaryConditionPDE &cn) const
{
    //cn = BoundaryConditionPDE::Robin(lambda, -1.0, lambda);
    //return environmentTemperature;
    cn = BoundaryConditionPDE::Neumann(1.0, 0.0);
    return 0.0;
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
    unsigned int ln = static_cast<unsigned int>(tn.i);
    DoubleVector z0(4), z1(4);
    PointNodeODE n0, n1;

    if (ln == timeDimension().min())
    {
        n0.i = static_cast<int>(ln+0); n0.x = n0.i*timeDimension().step();
        n1.i = static_cast<int>(ln+1); n1.x = n1.i*timeDimension().step();

        for (size_t i=0; i<solver->heatSourceNumber; i++)
        {
            const_this->i = i+1;
            double qi = 0.0;
            double vi = 0.0;

            const_this->start( z0, n0 );
            solver->sourceParams[i][ln].z = SpacePoint(z0[0], z0[1]);
            solver->sourceParams[i][ln].zt = SpacePoint(z0[2], z0[3]);
            const SpacePoint &zi = solver->sourceParams[i][ln].z;

            for (size_t j=0; j<solver->measrPointNumber; j++)
            {
                double nominal = 0.0;//solver->nominU[i][j];
                solver->measurePointValues[j].z = DeltaFunction::lumpedPoint4(u, solver->measurePoints[j], spaceDimensionX(), spaceDimensionY());
                solver->measurePointValues[j].x = solver->measurePointValues[j].y = 0.0;

                const SpacePoint &mp = solver->measurePoints[j];
                double dist = sqrt((zi.x - mp.x)*(zi.x - mp.x) + (zi.y - mp.y)*(zi.y - mp.y));
                qi += (solver->alpha1[i][j]*dist*dist + solver->alpha2[i][j]*dist + solver->alpha3[i][j]) * (solver->measurePointValues[j].z - nominal);
                vi += (solver->betta1[i][j]*dist*dist + solver->betta2[i][j]*dist + solver->betta3[i][j]) * (solver->measurePointValues[j].z - nominal);
            }
            solver->sourceParams[i][ln].q = qi;
            solver->sourceParams[i][ln].v = vi;
        }
    }

    if (ln != timeDimension().max())
    {
        for (size_t i=0; i<solver->heatSourceNumber; i++)
        {
            const_this->i = i+1;
            double qi = 0.0;
            double vi = 0.0;

            z0[0] = solver->sourceParams[i][ln].z.x;
            z0[1] = solver->sourceParams[i][ln].z.y;
            z0[2] = solver->sourceParams[i][ln].zt.x;
            z0[3] = solver->sourceParams[i][ln].zt.y;
            const_this->next( z0, n0, z1, n1, ODESolverMethod::EULER );
            solver->sourceParams[i][ln+1].z = SpacePoint(z1[0], z1[1]);
            solver->sourceParams[i][ln+1].zt = SpacePoint(z1[2], z1[3]);
            const SpacePoint &zi = solver->sourceParams[i][ln+1].z;

            for (size_t j=0; j<solver->measrPointNumber; j++)
            {
                double nominal = 0.0;//solver->nominU[i][j];
                solver->measurePointValues[j].z = DeltaFunction::lumpedPoint3(u, solver->measurePoints[j], spaceDimensionX(), spaceDimensionY());
                solver->measurePointValues[j].x = solver->measurePointValues[j].y = 0.0;

                const SpacePoint &mp = solver->measurePoints[j];
                double dist = sqrt((zi.x - mp.x)*(zi.x - mp.x) + (zi.y - mp.y)*(zi.y - mp.y));
                qi += (solver->alpha1[i][j]*dist*dist + solver->alpha2[i][j]*dist + solver->alpha3[i][j]) * (solver->measurePointValues[j].z - nominal);
                vi += (solver->betta1[i][j]*dist*dist + solver->betta2[i][j]*dist + solver->betta3[i][j]) * (solver->measurePointValues[j].z - nominal);
            }

            solver->sourceParams[i][ln+1].q = qi;
            solver->sourceParams[i][ln+1].v = vi;
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
        const SpacePoint &z = sourceParams[i][ln].z;
        //if (0.05 <= z.x <= 0.95 && 0.05 <= z.y <= 0.95)
            fx += sourceParams[i][ln].q * DeltaFunction::gaussian(sn, z, SpacePoint(spaceDimensionX().step(), spaceDimensionY().step()));
    }
    return fx ;
}

void Solver1::frw_layerInfo(const DoubleMatrix &u, const TimeNodePDE &tn) const
{
    unsigned int ln = static_cast<unsigned int>(tn.i);
    std::string msg = "[ " + std::to_string(tn.i+1) + " ]";
    msg += " [ " + std::to_string(sourceParams[0][tn.i].z.x) + ", " + std::to_string(sourceParams[0][tn.i].z.y) + " ]";
    msg += " [ " + std::to_string(sourceParams[1][tn.i].z.x) + ", " + std::to_string(sourceParams[1][tn.i].z.y) + " ]";
    msg += " [ " + std::to_string(sourceParams[0][tn.i].v) + ", " + std::to_string(sourceParams[0][tn.i].v) + " ]";
    msg += " [ " + std::to_string(sourceParams[1][tn.i].v) + ", " + std::to_string(sourceParams[1][tn.i].v) + " ]";
    msg += " [ " + std::to_string(sourceParams[0][tn.i].q) + ", " + std::to_string(sourceParams[0][tn.i].q) + " ]";
    msg += " [ " + std::to_string(sourceParams[1][tn.i].q) + ", " + std::to_string(sourceParams[1][tn.i].q) + " ]";
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

auto HeatEquationIBVP::A(const PointNodeODE &, size_t r, size_t c) const -> double { return solver->A[i-1][r-1][c-1]; }

auto HeatEquationIBVP::B(const PointNodeODE &, size_t r, size_t c) const -> double
{
    return solver->B[i-1][r-1][c-1];
}

auto HeatEquationIBVP::C(const PointNodeODE &node, size_t r) const -> double
{
    size_t ln = static_cast<size_t>(node.i);
    return solver->C[i-1][r-1] * solver->sourceParams[i-1][ln].v;
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

auto HeatEquationFBVP::layerInfo(const DoubleMatrix &u, const TimeNodePDE &tn) const -> void {}

auto HeatEquationFBVP::timeDimension() const -> Dimension { return solver->timeDimension(); }

auto HeatEquationFBVP::spaceDimensionX() const -> Dimension { return solver->spaceDimensionX(); }

auto HeatEquationFBVP::spaceDimensionY() const -> Dimension { return solver->spaceDimensionY(); }

auto HeatEquationFBVP::spaceDimensionZ() const -> Dimension { throw std::runtime_error(""); }

//--------------------------------------------------------------------------------------------------------------//

auto HeatEquationFBVP::A(const PointNodeODE &, size_t r, size_t c) const -> double { return -solver->A[i-1][r-1][c-1]; }

auto HeatEquationFBVP::B(const PointNodeODE &, size_t r, size_t c) const -> double { return -solver->B[i-1][r-1][c-1]; }

auto HeatEquationFBVP::C(const PointNodeODE &node, size_t r) const -> double
{
    size_t ln = static_cast<size_t>(node.i);
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
