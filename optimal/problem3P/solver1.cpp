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

    puts("----");
    s.backward.setThermalDiffusivity(-0.01);
    s.backward.setThermalConductivity(0.0);
    s.backward.setThermalConvection(s.lambda0);
    s.backward.implicit_calculate_D2V1();

    DoubleVector x;
    s.parameterToVector(x);
    DoubleVector g;
    s.gradient(x, g);
    IPrinter::print(g, g.length());
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

    alpha1.resize(heatSourceNumber, measrPointNumber, 0.200);
    alpha2.resize(heatSourceNumber, measrPointNumber, 0.250);
//    alpha1.resize(heatSourceNumber, measrPointNumber, 0.000);
//    alpha2.resize(heatSourceNumber, measrPointNumber, 0.000);
    alpha3.resize(heatSourceNumber, measrPointNumber, 0.275);
    betta1.resize(heatSourceNumber, measrPointNumber, 0.200);
    betta2.resize(heatSourceNumber, measrPointNumber, 0.250);
//    betta1.resize(heatSourceNumber, measrPointNumber, 0.000);
//    betta2.resize(heatSourceNumber, measrPointNumber, 0.000);
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
        pp.ft = new SpacePoint[heatSourceNumber];
        pp.p = new SpacePoint[heatSourceNumber];

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

void Solver1::bcw_saveToImage(const DoubleMatrix &p UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM) const
{
#ifdef USE_LIB_IMAGING
    static double MIN = +10000.0;
    static double MAX = -10000.0;
    if (p.max()>MAX) MAX = p.max();
    if (p.min()<MIN) MIN = p.min();
    printf(">>> %6d %14.6f %14.6f %14.6f %14.6f\n", tn.i, p.min(), p.max(), MIN, MAX);

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

double Solver1::bcw_final(const SpaceNodePDE &sn, FinalCondition condition) const
{
    size_t i = static_cast<size_t>(sn.i);
    size_t j = static_cast<size_t>(sn.j);
    return -2.0*(U[j][i]-V[j][i]);
}

double Solver1::bcw_boundary(const SpaceNodePDE &, const TimeNodePDE &, BoundaryConditionPDE &cn) const
{
    //cn = BoundaryConditionPDE(BoundaryCondition::Neumann, 1.0, 0.0, 0.0);
    cn = BoundaryConditionPDE::Neumann(1.0, 0.0);
    return 0.0;
}

auto Solver1::bcw_f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double
{
    double fx = 0.0;
    unsigned int ln = static_cast<unsigned int>(tn.i);
    PointNodeODE node; node.i = tn.i; node.x = tn.t;
    const ProblemParams &pp = sourceParams[ln];

    for (size_t j=0; j<measrPointNumber; j++)
    {
        const SpacePoint &mp = measurePoints[j];
        double sum = 0.0;
        for (size_t i=0; i<heatSourceNumber; i++)
        {
            const SpacePoint &zi = pp.z[i];
            double dist = sqrt((zi.x - mp.x)*(zi.x - mp.x) + (zi.y - mp.y)*(zi.y - mp.y));

            sum += (alpha1[i][j]*dist*dist + alpha2[i][j]*dist + alpha3[i][j]) * ( A4(node, 1, i+1)*pp.f[i].x + A4(node, 2, i+1)*pp.f[i].y );
            sum += (betta1[i][j]*dist*dist + betta2[i][j]*dist + betta3[i][j]) * ( pp.p[i].z );
        }
        fx += sum * DeltaFunction::gaussian(sn, measurePoints[j],
                                            SpacePoint(spaceDimensionX().step(), spaceDimensionY().step()));
    }
    return fx;
}

void Solver1::bcw_layerInfo(const DoubleMatrix &p, const TimeNodePDE &tn) const
{
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
            ProblemParams &param = solver->sourceParams[ln];
            double qi = 0.0;
            double vi = 0.0;

            const_this->start( z0, n0 );
            param.z[i] = SpacePoint(z0[0], z0[1]);
            param.zt[i] = SpacePoint(z0[2], z0[3]);
            const SpacePoint &zi = param.z[i];

            for (size_t j=0; j<solver->measrPointNumber; j++)
            {
                double nominal = solver->nU.at(i, j);
                SpacePoint d;
                param.u[j].z = DeltaFunction::lumpedPoint4(u, solver->measurePoints[j], spaceDimensionX(), spaceDimensionY(), d);
                param.u[j].x = d.x; param.u[j].y = d.y;

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
                SpacePoint d;
                pp1.u[j].z = DeltaFunction::lumpedPoint4(u, solver->measurePoints[j], spaceDimensionX(), spaceDimensionY(), d);
                pp1.u[j].x = d.x; pp1.u[j].y = d.y;

                const SpacePoint &mp = solver->measurePoints[j];
                double dist = sqrt((zi.x - mp.x)*(zi.x - mp.x) + (zi.y - mp.y)*(zi.y - mp.y));
                qi += (solver->alpha1[i][j]*dist*dist + solver->alpha2[i][j]*dist + solver->alpha3[i][j]) * (pp1.u[j].z - nominal);
                vi += (solver->betta1[i][j]*dist*dist + solver->betta2[i][j]*dist + solver->betta3[i][j]) * (pp1.u[j].z - nominal);
            }
            pp1.q[i] = qi;
            pp1.v[i] = vi;
        }

        solver->U = u;
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
        if ( (0.0 <= z.x <= 1.0) && (0.0 <= z.y <= 1.0) )
        {
            fx += sourceParams[ln].q[i] * DeltaFunction::gaussian(sn, z,
                  SpacePoint(spaceDimensionX().step(), spaceDimensionY().step()));
        }
    }
    return fx ;
}

void Solver1::frw_layerInfo(const DoubleMatrix &u, const TimeNodePDE &tn) const
{
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
    //frw_saveToImage(u, tn);
}

void Solver1::gradient(const DoubleVector &x, DoubleVector &g) const
{
    const_cast<Solver1*>(this)->vectorToParameter(x);
    forward.implicit_calculate_D2V1();
    backward.implicit_calculate_D2V1();

    size_t L = timeDimension().size()-1;
    double ht = timeDimension().step();

    size_t size = heatSourceNumber*measrPointNumber;
    g.clear();
    g.resize(size*6);

    for (size_t i=0; i<heatSourceNumber; i++)
    {
        for (size_t j=0; j<measrPointNumber; j++)
        {
            const SpacePoint &mp = measurePoints[j];
            const SpacePoint &zi = sourceParams[0].z[i];
            double dist = sqrt((zi.x - mp.x)*(zi.x - mp.x) + (zi.y - mp.y)*(zi.y - mp.y));

            size_t ln = 0;
            PointNodeODE node; node.i  = ln; node.x = ln*ht;
            g[0*size + i*heatSourceNumber  + j] = 0.5*(A4(node, 1, i+1)*sourceParams[ln].f[i].x + A4(node, 2, i+1)*sourceParams[ln].f[i].y)*(dist*dist)*(sourceParams[ln].u[j].z-nU[i][j]);
            g[1*size + i*heatSourceNumber  + j] = 0.5*(A4(node, 1, i+1)*sourceParams[ln].f[i].x + A4(node, 2, i+1)*sourceParams[ln].f[i].y)*(dist)*(sourceParams[ln].u[j].z-nU[i][j]);
            g[2*size + i*heatSourceNumber  + j] = 0.5*(A4(node, 1, i+1)*sourceParams[ln].f[i].x + A4(node, 2, i+1)*sourceParams[ln].f[i].y)*(sourceParams[ln].u[j].z-nU[i][j]);
            g[3*size + i*heatSourceNumber  + j] = 0.5*(sourceParams[ln].p[i].z)*(dist*dist)*(sourceParams[ln].u[j].z-nU[i][j]);
            g[4*size + i*heatSourceNumber  + j] = 0.5*(sourceParams[ln].p[i].z)*(dist)*(sourceParams[ln].u[j].z-nU[i][j]);
            g[5*size + i*heatSourceNumber  + j] = 0.5*(sourceParams[ln].p[i].z)*(sourceParams[ln].u[j].z-nU[i][j]);

            for (size_t ln=1; ln<L; ln++)
            {
                node; node.i  = ln; node.x = ln*ht;
                g[0*size + i*heatSourceNumber  + j] += (A4(node, 1, i+1)*sourceParams[ln].f[i].x + A4(node, 2, i+1)*sourceParams[ln].f[i].y)*(dist*dist)*(sourceParams[ln].u[j].z-nU[i][j]);
                g[1*size + i*heatSourceNumber  + j] += (A4(node, 1, i+1)*sourceParams[ln].f[i].x + A4(node, 2, i+1)*sourceParams[ln].f[i].y)*(dist)*(sourceParams[ln].u[j].z-nU[i][j]);
                g[2*size + i*heatSourceNumber  + j] += (A4(node, 1, i+1)*sourceParams[ln].f[i].x + A4(node, 2, i+1)*sourceParams[ln].f[i].y)*(sourceParams[ln].u[j].z-nU[i][j]);
                g[3*size + i*heatSourceNumber  + j] += (sourceParams[ln].p[i].z)*(dist*dist)*(sourceParams[ln].u[j].z-nU[i][j]);
                g[4*size + i*heatSourceNumber  + j] += (sourceParams[ln].p[i].z)*(dist)*(sourceParams[ln].u[j].z-nU[i][j]);
                g[5*size + i*heatSourceNumber  + j] += (sourceParams[ln].p[i].z)*(sourceParams[ln].u[j].z-nU[i][j]);
            }

            ln = L;
            node; node.i  = ln; node.x = ln*ht;
            g[0*size + i*heatSourceNumber  + j] += 0.5*(A4(node, 1, i+1)*sourceParams[ln].f[i].x + A4(node, 2, i+1)*sourceParams[ln].f[i].y)*(dist*dist)*(sourceParams[ln].u[j].z-nU[i][j]);
            g[1*size + i*heatSourceNumber  + j] += 0.5*(A4(node, 1, i+1)*sourceParams[ln].f[i].x + A4(node, 2, i+1)*sourceParams[ln].f[i].y)*(dist)*(sourceParams[ln].u[j].z-nU[i][j]);
            g[2*size + i*heatSourceNumber  + j] += 0.5*(A4(node, 1, i+1)*sourceParams[ln].f[i].x + A4(node, 2, i+1)*sourceParams[ln].f[i].y)*(sourceParams[ln].u[j].z-nU[i][j]);
            g[3*size + i*heatSourceNumber  + j] += 0.5*(sourceParams[ln].p[i].z)*(dist*dist)*(sourceParams[ln].u[j].z-nU[i][j]);
            g[4*size + i*heatSourceNumber  + j] += 0.5*(sourceParams[ln].p[i].z)*(dist)*(sourceParams[ln].u[j].z-nU[i][j]);
            g[5*size + i*heatSourceNumber  + j] += 0.5*(sourceParams[ln].p[i].z)*(sourceParams[ln].u[j].z-nU[i][j]);

            g[0*size + i*heatSourceNumber  + j] *= ht;
            g[1*size + i*heatSourceNumber  + j] *= ht;
            g[2*size + i*heatSourceNumber  + j] *= ht;
            g[3*size + i*heatSourceNumber  + j] *= ht;
            g[4*size + i*heatSourceNumber  + j] *= ht;
            g[5*size + i*heatSourceNumber  + j] *= ht;
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
        for (size_t i=0; i<solver->heatSourceNumber; i++)
        {
            const_this->i = i+1;
            ProblemParams &param = solver->sourceParams[ln];
            const_this->start( f0, n0 );
            param.f[i]  = SpacePoint(f0[0], f0[1]);
            param.ft[i] = SpacePoint(f0[2], f0[3]);

            SpacePoint d;
            param.p[i].z = DeltaFunction::lumpedPoint4(p, param.z[i], spaceDimensionX(), spaceDimensionY(), d);
            param.p[i].x = d.x; param.p[i].y = d.y;
        }
    }

    if (ln != timeDimension().min())
    {
        for (size_t i=0; i<solver->heatSourceNumber; i++)
        {
            const_this->i = i+1;
            ProblemParams &pp0 = solver->sourceParams[ln];
            ProblemParams &pp1 = solver->sourceParams[ln-1];
            f0[0] = pp0.f[i].x; f0[1] = pp0.f[i].y; f0[2] = pp0.ft[i].x; f0[3] = pp0.ft[i].y;
            const_this->next( f0, n0, f1, n1, ODESolverMethod::EULER );
            pp1.f[i]  = SpacePoint(f1[0], f1[1]);
            pp1.ft[i] = SpacePoint(f1[2], f1[3]);

            SpacePoint d;
            pp1.p[i].z = DeltaFunction::lumpedPoint4(p, pp1.z[i], spaceDimensionX(), spaceDimensionY(), d);
            pp1.p[i].x = d.x; pp1.p[i].y = d.y;
        }
    }

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
    size_t ln = static_cast<size_t>(node.i);
    const ProblemParams &param = solver->sourceParams[ln];
    return (r==1) ? param.p[i-1].x * param.q[i-1] : param.p[i-1].y * param.q[i-1];
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
            alpha1[i][j] = x[0*size + i*heatSourceNumber  + j];
            alpha2[i][j] = x[1*size + i*heatSourceNumber  + j];
            alpha3[i][j] = x[2*size + i*heatSourceNumber  + j];
            betta1[i][j] = x[3*size + i*heatSourceNumber  + j];
            betta2[i][j] = x[4*size + i*heatSourceNumber  + j];
            betta3[i][j] = x[5*size + i*heatSourceNumber  + j];
        }
    }
}

void Solver1::parameterToVector(DoubleVector &x)
{
    size_t size = heatSourceNumber*measrPointNumber;
    x.clear();
    x.resize(size*6);

    for (size_t i=0; i<heatSourceNumber; i++)
    {
        for (size_t j=0; j<measrPointNumber; j++)
        {
            x[0*size + i*heatSourceNumber  + j] = alpha1[i][j];
            x[1*size + i*heatSourceNumber  + j] = alpha2[i][j];
            x[2*size + i*heatSourceNumber  + j] = alpha3[i][j];
            x[3*size + i*heatSourceNumber  + j] = betta1[i][j];
            x[4*size + i*heatSourceNumber  + j] = betta2[i][j];
            x[5*size + i*heatSourceNumber  + j] = betta3[i][j];
        }
    }
}


