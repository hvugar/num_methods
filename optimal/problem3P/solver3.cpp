#include "solver3.h"

//#define OMTIMZIE_V

using namespace p3p3;

void Functional::Main(int /*argc*/, char** /*argv*/)
{
    Functional *fn = new Functional();

    const size_t time_size = fn->timeDimension().size();
    const size_t dimX_size = fn->spaceDimensionX().size();
    const size_t dimY_size = fn->spaceDimensionY().size();0.0;//
    const double time_step = fn->timeDimension().step();
    //const double dimX_step = fn->spaceDimensionX().step();
    //const double dimY_step = fn->spaceDimensionX().step();

    fn->V.resize(dimY_size, dimX_size, 5.0);
    fn->U.resize(dimY_size, dimX_size, 0.0);

#ifdef OMTIMZIE_V
    DoubleVector x(4*time_size);
#else
    DoubleVector x(6*time_size);
#endif

    for (unsigned int ln=0; ln<time_size; ln++)
    {
        const double t = ln*time_step;
        const TimeNodePDE tn(ln, t);

        x[0*time_size + ln] = 50.00;
        x[1*time_size + ln] = 50.00;
#ifdef OMTIMZIE_V
        x[2*time_size + ln] = 2.0*t;
        x[3*time_size + ln] = 2.0*t;
#else
        x[2*time_size + ln] = /*+0.8*t+0.1;*/+fn->R[0]*sin(fn->v(tn, 0)*t) + 0.50;
        x[3*time_size + ln] = /*+0.8*t+0.1;*/+fn->R[0]*cos(fn->v(tn, 0)*t) + 0.50;
        x[4*time_size + ln] = /*+0.8*t+0.1;*/-fn->R[1]*sin(fn->v(tn, 1)*t) + 0.50;
        x[5*time_size + ln] = /*-0.8*t+0.9;*/+fn->R[1]*cos(fn->v(tn, 1)*t) + 0.50;
#endif
    }

    unsigned int w = 10, p = 4;
    {
        DoubleVector g0;
        fn->gradient(x, g0);
        IPrinter::printVector(w, p, g0.mid(0,   100).L2Normalize());
        IPrinter::printVector(w, p, g0.mid(101, 201).L2Normalize());
        IPrinter::printVector(w, p, g0.mid(202, 302).L2Normalize());
        IPrinter::printVector(w, p, g0.mid(303, 403).L2Normalize());
#ifndef OMTIMZIE_V
        IPrinter::printVector(w, p, g0.mid(404, 504).L2Normalize());
        IPrinter::printVector(w, p, g0.mid(505, 605).L2Normalize());
#endif
        IPrinter::printSeperatorLine();

//#ifndef OMTIMZIE_V
//        IPrinter::printVector(w, p, fn->pp[0]/*.L2Normalize()*/);
//        IPrinter::printVector(w, p, fn->px[0]/*.L2Normalize()*/);
//        IPrinter::printVector(w, p, fn->py[0]/*.L2Normalize()*/);

//        IPrinter::printVector(w, p, fn->pp[1]/*.L2Normalize()*/);
//        IPrinter::printVector(w, p, fn->px[1]/*.L2Normalize()*/);
//        IPrinter::printVector(w, p, fn->py[1]/*.L2Normalize()*/);
//#endif


        //        IPrinter::printVector(10, 6, g0.mid(0,   100));
        //        IPrinter::printVector(10, 6, g0.mid(101, 201));
        //        IPrinter::printVector(10, 6, g0.mid(202, 302));
        //        IPrinter::printVector(10, 6, g0.mid(303, 403));
        //        IPrinter::printSeperatorLine();
    }

    {
        DoubleVector g1(x.length());

        const double step = 0.01;
        IGradient::Gradient(fn, step, x, g1, 0, 100);   g1[  0] *= 2.0; g1[100] *= 2.0; IPrinter::printVector(w, p, g1.mid(  0, 100).L2Normalize());
        IGradient::Gradient(fn, step, x, g1, 101, 201); g1[101] *= 2.0; g1[201] *= 2.0; IPrinter::printVector(w, p, g1.mid(101, 201).L2Normalize());
        IGradient::Gradient(fn, step, x, g1, 202, 302); g1[202] *= 2.0; g1[302] *= 2.0; IPrinter::printVector(w, p, g1.mid(202, 302).L2Normalize());
        IGradient::Gradient(fn, step, x, g1, 303, 403); g1[303] *= 2.0; g1[403] *= 2.0; IPrinter::printVector(w, p, g1.mid(303, 403).L2Normalize());
#ifndef OMTIMZIE_V
        IGradient::Gradient(fn, step, x, g1, 404, 504); g1[404] *= 2.0; g1[504] *= 2.0; IPrinter::printVector(w, p, g1.mid(404, 504).L2Normalize());
        IGradient::Gradient(fn, step, x, g1, 505, 605); g1[505] *= 2.0; g1[605] *= 2.0; IPrinter::printVector(w, p, g1.mid(505, 605).L2Normalize());
#endif
        //        IPrinter::printSeperatorLine();

        //        IPrinter::printVector(10, 6, g1.mid(0, 100));
        //        IPrinter::printVector(10, 6, g1.mid(101, 201));
        //        IPrinter::printVector(10, 6, g1.mid(202, 302));
        //        IPrinter::printVector(10, 6, g1.mid(303, 403));
    }
}

Functional::Functional()
{
    ih = new HeatEquationIBVP(this);
    fh = new HeatEquationFBVP(this);

    const double a = 1.00;
    const double lambda0 = 0.0;

    ih->setThermalDiffusivity(a);
    ih->setThermalConductivity(0.0);
    ih->setThermalConvection(-lambda0);

    fh->setThermalDiffusivity(-a);
    fh->setThermalConductivity(0.0);
    fh->setThermalConvection(lambda0);
}

SpacePoint Functional::tr(const TimeNodePDE &tn, size_t i) const
{
    const double t = tn.t;
    const size_t time_size = timeDimension().size();
    SpacePoint sp;

#ifdef OMTIMZIE_V
    switch (i)
    {

    case 0:
    {
        sp.x = +R[0]*sin(v(tn, i)*t) + 0.50;
        sp.y = +R[0]*cos(v(tn, i)*t) + 0.50;
    } break;

    case 1:
    {
        sp.x = -R[1]*sin(v(tn,i)*t) + 0.50;
        sp.y = +R[1]*cos(v(tn,i)*t) + 0.50;
    } break;

    default:
    {
        //throw std::exception("i is not valid...");
    }

    }
#else
    if (i==0) { sp.x = x[2*time_size + tn.i]; sp.y = x[3*time_size + tn.i]; }
    if (i==1) { sp.x = x[4*time_size + tn.i]; sp.y = x[5*time_size + tn.i]; }
#endif

    return sp;
}

double Functional::v(const TimeNodePDE &tn, size_t i) const
{
#ifndef OPTIMIZE_V
    return 2.0*tn.t;
#else
    const size_t time_size = timeDimension().size();
    return x[(2+i)*time_size + tn.i];
#endif
}

double Functional::q(const TimeNodePDE &tn, size_t i) const
{
    const size_t time_size = timeDimension().size();
    return x[i*time_size + tn.i];
}

double Functional::fx(const DoubleVector &x) const
{
    setVector(x);

    ih->implicit_calculate_D2V1();

    return integral(U);
}

auto Functional::integral(const DoubleMatrix &U) const -> double
{
    const Dimension &dimensionX = spaceDimensionX();
    const Dimension &dimensionY = spaceDimensionY();
    const size_t N = static_cast<size_t> ( dimensionX.size()-1 );
    const size_t M = static_cast<size_t> ( dimensionY.size()-1 );
    const double hx = dimensionX.step();
    const double hy = dimensionY.step();

    double udiff = 0.0;
    double usum = 0.0;

    udiff = (U[0][0]-V[0][0]); usum += 0.25 * udiff * udiff;// * mu(0, 0);
    udiff = (U[0][N]-V[0][N]); usum += 0.25 * udiff * udiff;// * mu(N, 0);
    udiff = (U[M][0]-V[M][0]); usum += 0.25 * udiff * udiff;// * mu(0, M);
    udiff = (U[M][N]-V[M][N]); usum += 0.25 * udiff * udiff;// * mu(N, M);

    for (size_t n=1; n<=N-1; n++)
    {
        udiff = U[0][n]-V[0][n]; usum += 0.5 * udiff * udiff;// * mu(n, 0);
        udiff = U[M][n]-V[M][n]; usum += 0.5 * udiff * udiff;// * mu(n, M);
    }

    for (size_t m=1; m<=M-1; m++)
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

void Functional::gradient(const DoubleVector &x, DoubleVector &g) const
{
    setVector(x);

    const double time_step = timeDimension().step();
    const unsigned int time_size = timeDimension().size();

    ih->implicit_calculate_D2V1();
    fh->implicit_calculate_D2V1();

    g.clear();
    g.resize(x.length());

    for (unsigned int ln=0; ln<time_size; ln++)
    {
        const double t = ln*time_step;
        const TimeNodePDE tn = TimeNodePDE(ln, t);

        for (size_t i=0; i<heat_source_number; i++)
        {
            const double qi = q(tn, i);

            g[i*time_size+ln] = -pp[i][ln];


#ifdef OMTIMZIE_V
            if (i==0) g[(2+i)*time_size+ln] = -qi * R[0] * ( +px[0][ln]*cos(v(tn, 0)*t)*t - py[0][ln]*sin(v(tn, 0)*t)*t );
            if (i==1) g[(2+i)*time_size+ln] = -qi * R[1] * ( -px[1][ln]*cos(v(tn, 1)*t)*t - py[1][ln]*sin(v(tn, 1)*t)*t );
#else
            if (i==0)
            {
                g[2*time_size+ln] = -qi * px[0][ln];
                g[3*time_size+ln] = -qi * py[0][ln];
            }

            if (i==1)
            {
                g[4*time_size+ln] = -qi * px[1][ln];
                g[5*time_size+ln] = -qi * py[1][ln];
            }
#endif
        }
    }
}

void Functional::setVector(const DoubleVector &x) const
{
    const_cast<Functional*>(this)->x = x;
}

void Functional::frw_saveToImage(const DoubleMatrix &u UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM) const
{
#ifdef USE_LIB_IMAGING
    static double MIN = +10000.0;
    static double MAX = -10000.0;
    if (u.max()>MAX) MAX = u.max();
    if (u.min()<MIN) MIN = u.min();
    printf(">>> %6d %14.6f %14.6f %14.6f %14.6f\n", tn.i, u.min(), u.max(), MIN, MAX);

    //if (tn.i%10==0)
    {
        QString filename = QString("data/problem3P/f/png/%1.png").arg(tn.i, 8, 10, QChar('0'));
        QPixmap pixmap;
        //visualizeMatrixHeat(u, 0.0, 46.052, pixmap, 101, 101);
        //visualizeMatrixHeat(u, 1.50, 2.1807, pixmap, spaceDimensionX().size(), spaceDimensionY().size());
        //visualizeMatrixHeat(u, 0.50, 2.32, pixmap, spaceDimensionX().size(), spaceDimensionY().size());
        visualizeMatrixHeat(u, u.min(), u.max(), pixmap, spaceDimensionX().size(), spaceDimensionY().size());

        //        QPainter painter(&pixmap);
        //        painter.setFont(QFont("Consolas", 12));
        //        painter.setPen(Qt::white);

        //        painter.drawText(2, 30, QString("time: ")+QString::number(tn.t, 'f', 3));
        //        painter.drawText(2, 60, QString("temp: ")+QString::number(u.max(), 'f', 3)+QString(" max temp: ")+QString::number(MAX, 'f', 3));


        pixmap.save(filename);
    }
#endif
}

/********************************************************************************************************************************************************/

HeatEquationIBVP::HeatEquationIBVP(Functional *function) : _functional(function) {}

HeatEquationIBVP::HeatEquationIBVP(const HeatEquationIBVP &other) { *this = other; }

HeatEquationIBVP& HeatEquationIBVP::operator =(const HeatEquationIBVP &) { return *this; }

HeatEquationIBVP::~HeatEquationIBVP() {}

auto HeatEquationIBVP::initial(const SpaceNodePDE &/*sn*/, InitialCondition /*ic*/) const -> double { return _initial_temperature; }

auto HeatEquationIBVP::boundary(const SpaceNodePDE &/*sn*/, const TimeNodePDE &/*tn*/, BoundaryConditionPDE &bc) const -> double
{
    const double lambda1 = _functional->_lambda1;
    //bc = BoundaryConditionPDE::Dirichlet(1.0, 1.0); return 0.0;
    //bc = BoundaryConditionPDE::Neumann(1.0, 0.0); return 0.0;
    bc = BoundaryConditionPDE::Robin(lambda1, +1.0, lambda1); return _enviroment_temperature;
}

auto HeatEquationIBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double
{
    //size_t ln = static_cast<size_t>(tn.i);

    double fx = -thermalConvection() * _enviroment_temperature;
    const double dimX_step = spaceDimensionX().step();
    const double dimY_step = spaceDimensionY().step();
    //const unsigned dimX_size = spaceDimensionX().size();
    //const unsigned dimY_size = spaceDimensionY().size();


    double sum = 0.0;
    for (size_t i=0; i<_functional->heat_source_number; i++)
    {
        const SpacePoint &zi = _functional->tr(tn, i);
        const double qi = _functional->q(tn, i);
        sum += qi * DeltaFunction::gaussian(sn, zi, SpacePoint(dimX_step, dimY_step));
        //sum += qi * DeltaFunction::nearest(sn, zi, dimX_step, dimY_step, dimX_size, dimY_size);
    }

    return fx + sum;
}

auto HeatEquationIBVP::layerInfo(const DoubleMatrix &u, const TimeNodePDE &tn) const -> void
{
    if (static_cast<int>(tn.i) == timeDimension().max()) { _functional->U = u; }
    //if (tn.i%100!=0) return;

    //IPrinter::printSeperatorLine();
    //IPrinter::printMatrix(u);

    //frw_saveToImage(u, tn);
}

auto HeatEquationIBVP::timeDimension() const -> Dimension { return _functional->timeDimension(); }
auto HeatEquationIBVP::spaceDimensionX() const -> Dimension { return _functional->spaceDimensionX(); }
auto HeatEquationIBVP::spaceDimensionY() const -> Dimension { return _functional->spaceDimensionY(); }
auto HeatEquationIBVP::spaceDimensionZ() const -> Dimension { return _functional->spaceDimensionZ(); }

/********************************************************************************************************************************************************/

HeatEquationFBVP::HeatEquationFBVP(Functional *function) :_functional(function) {}

HeatEquationFBVP::HeatEquationFBVP(const HeatEquationFBVP &other) { *this = other; }

HeatEquationFBVP& HeatEquationFBVP::operator =(const HeatEquationFBVP &) { return *this; }

HeatEquationFBVP::~HeatEquationFBVP() {}

auto HeatEquationFBVP::final(const SpaceNodePDE &sn, FinalCondition /*condition*/) const -> double
{
    const size_t i = static_cast<size_t>(sn.i), j = static_cast<size_t>(sn.j);
    return -2.0*(_functional->U[j][i] - _functional->V[j][i]);
}

auto HeatEquationFBVP::boundary(const SpaceNodePDE &/*sn*/, const TimeNodePDE &/*tn*/, BoundaryConditionPDE &bc) const -> double
{
    const double lambda1 = _functional->_lambda1;
    //bc = BoundaryConditionPDE::Dirichlet(1.0, 1.0); return 0.0;
    //bc = BoundaryConditionPDE::Neumann(1.0, 0.0); return 0.0;
    bc = BoundaryConditionPDE::Robin(lambda1, +1.0, lambda1); return 0.0;
}

auto HeatEquationFBVP::f(const SpaceNodePDE &/*sn*/, const TimeNodePDE &/*tn*/) const -> double
{
    double fx = 0.0;
    return fx;
}

auto HeatEquationFBVP::layerInfo(const DoubleMatrix &p, const TimeNodePDE &tn) const -> void
{
    if (_functional->pp == nullptr)
    {
        _functional->pp = new DoubleVector[_functional->heat_source_number];
        _functional->px = new DoubleVector[_functional->heat_source_number];
        _functional->py = new DoubleVector[_functional->heat_source_number];

        for (size_t i=0; i<_functional->heat_source_number; i++)
        {
            _functional->pp[i].resize(timeDimension().size());
            _functional->px[i].resize(timeDimension().size());
            _functional->py[i].resize(timeDimension().size());
        }
    }

    const double dimX_step = _functional->spaceDimensionX().step();
    const double dimY_step = _functional->spaceDimensionX().step();
    const size_t dimX_size = _functional->spaceDimensionX().size();
    const size_t dimY_size = _functional->spaceDimensionY().size();

    for (size_t i=0; i<_functional->heat_source_number; i++)
    {
        SpacePoint spi = _functional->tr(tn, i);

        //const size_t rx = static_cast<unsigned int>(round(spi.x * (dimX_size-1)));
        //const size_t ry = static_cast<unsigned int>(round(spi.y * (dimY_size-1)));

        //_functional->pp[i][tn.i] = DeltaFunction::lumpedPoint2(p, spi, spaceDimensionX(), spaceDimensionY());
        double dx, dy;
        _functional->pp[i][tn.i] = DeltaFunction::lumpedPointG(p, spi, spaceDimensionX(), spaceDimensionY(), 1, 4, dx, dy);
        _functional->px[i][tn.i] = dx;
        _functional->py[i][tn.i] = dy;

        //_functional->px[i][tn.i] = (p[ry][rx+1]-p[ry][rx-1])/(2.0*dimX_step);
        //_functional->py[i][tn.i] = (p[ry+1][rx]-p[ry-1][rx])/(2.0*dimY_step);
    }
}

auto HeatEquationFBVP::timeDimension() const -> Dimension { return _functional->timeDimension(); }
auto HeatEquationFBVP::spaceDimensionX() const -> Dimension { return _functional->spaceDimensionX(); }
auto HeatEquationFBVP::spaceDimensionY() const -> Dimension { return _functional->spaceDimensionY(); }
auto HeatEquationFBVP::spaceDimensionZ() const -> Dimension { return _functional->spaceDimensionZ(); }
