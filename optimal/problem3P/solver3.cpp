#include "solver3.h"

using namespace p3p3;

void Functional::Main(int /*argc*/, char** /*argv*/)
{
    Functional *functional = new Functional();

    unsigned int time_size = functional->timeDimension().size();
    unsigned int dimX_size = functional->spaceDimensionX().size();
    unsigned int dimY_size = functional->spaceDimensionY().size();
    functional->U.resize(dimY_size, dimX_size, 2.0);
    functional->uT.resize(dimY_size, dimX_size, 0.0);
    DoubleVector x(2*time_size);

    //functional->fx(x);

    DoubleVector g;
    functional->gradient(x, g);

    DoubleVector g1(x.length());
    //IGradient::Gradient(functional, 0.01, x, g1);

    IPrinter::printSeperatorLine();

    g.L2Normalize();
    g1.L2Normalize();

    IPrinter::printVector(g);
    IPrinter::printVector(g1);


}

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
    size_t ln = static_cast<size_t>(tn.i);

    double fx = -thermalConvection() * _enviroment_temperature;

    double sum = 0.0;
    for (size_t i=0; i<_functional->heat_source_number; i++)
    {
        printf("time: %d\n", ln);
        const SpacePoint &zi = _functional->tr(tn, i);
        const double qi = _functional->q(tn, i);
        //sum += qi * DeltaFunction::gaussian(sn, zi, SpacePoint(spaceDimensionX().step(), spaceDimensionY().step()));
        sum += qi * DeltaFunction::nearest(sn, zi, spaceDimensionX().step(), spaceDimensionY().step(), spaceDimensionX().size(), spaceDimensionY().size());
    }

    return fx + sum;

    //if (sn.i==250 && sn.j == 250) { return  1.0/(0.002*0.002) * 0.1; } else { return  0.0; }
}



double Functional::q(const TimeNodePDE &tn, size_t i) const
{
    //    switch (i)
    //    {
    //    case 0: { return 0.05; }
    //    case 1: { return 0.05; }
    //    default: throw std::exception("error");
    //    }

    unsigned int size = timeDimension().size();
    return x[i*size + tn.i];
}



auto HeatEquationIBVP::layerInfo(const DoubleMatrix &u, const TimeNodePDE &tn) const -> void
{
    //if (tn.i%100!=0) return;

    //IPrinter::printSeperatorLine();
    //IPrinter::printMatrix(u);

    //frw_saveToImage(u, tn);

    if (tn.i == timeDimension().max()) { _functional->uT = u; }
}

auto HeatEquationIBVP::timeDimension() const -> Dimension { return _functional->timeDimension(); }
auto HeatEquationIBVP::spaceDimensionX() const -> Dimension { return _functional->spaceDimensionX(); }
auto HeatEquationIBVP::spaceDimensionY() const -> Dimension { return _functional->spaceDimensionY(); }
auto HeatEquationIBVP::spaceDimensionZ() const -> Dimension { return _functional->spaceDimensionZ(); }

Functional::Functional()
{
    ih = new HeatEquationIBVP(this);
    fh = new HeatEquationFBVP(this);

    double _a = 1.0;
    double _lambda0 = 1.0;

    ih->setThermalDiffusivity(_a);
    ih->setThermalConductivity(0.0);
    ih->setThermalConvection(-_lambda0);//-1.0 heating

    fh->setThermalDiffusivity(-_a);
    fh->setThermalConductivity(0.0);
    fh->setThermalConvection(_lambda0);//+1.0 heating
}

SpacePoint Functional::tr(const TimeNodePDE &tn, size_t i) const
{
    const double R = 0.40;
    SpacePoint sp;

    switch (i)
    {

    case 0:
    {
        sp.x = R*sin(v(tn, i)) + 0.50;
        sp.y = R*cos(v(tn, i)) + 0.50;
    } break;

    case 1:
    {
        sp.x = 0.5*R*sin(-v(tn,i)-2) + 0.50;
        sp.y = 0.5*R*cos(-v(tn,i)-2) + 0.50;
    } break;

    default:
    {
        throw std::exception("i is not valid...");
    }

    }

    return sp;
}

double Functional::v(const TimeNodePDE &tn, size_t /*i*/) const
{
    return 2.0*tn.t;
}

void HeatEquationIBVP::frw_saveToImage(const DoubleMatrix &u UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM) const
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


HeatEquationFBVP::HeatEquationFBVP(Functional *function) :_functional(function) {}

HeatEquationFBVP::HeatEquationFBVP(const HeatEquationFBVP &other) { *this = other; }

HeatEquationFBVP& HeatEquationFBVP::operator =(const HeatEquationFBVP &) { return *this; }

HeatEquationFBVP::~HeatEquationFBVP() {}

auto HeatEquationFBVP::final(const SpaceNodePDE &sn, FinalCondition /*condition*/) const -> double
{
    return -2.0*(_functional->uT[sn.j][sn.i] - _functional->U[sn.j][sn.i]);
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

    for (size_t i=0; i<_functional->heat_source_number; i++)
    {
        SpacePoint spi = _functional->tr(tn, i);
        _functional->pp[i][tn.i] = DeltaFunction::lumpedPoint2(p, spi, spaceDimensionX(), spaceDimensionY());
        //_functional->px[i][tn.i] = DeltaFunction::lumpedPoint2(p, spi, spaceDimensionX(), spaceDimensionY());
        //_functional->py[i][tn.i] = DeltaFunction::lumpedPoint2(p, spi, spaceDimensionX(), spaceDimensionY());
    }

}

auto HeatEquationFBVP::timeDimension() const -> Dimension { return _functional->timeDimension(); }
auto HeatEquationFBVP::spaceDimensionX() const -> Dimension { return _functional->spaceDimensionX(); }
auto HeatEquationFBVP::spaceDimensionY() const -> Dimension { return _functional->spaceDimensionY(); }
auto HeatEquationFBVP::spaceDimensionZ() const -> Dimension { return _functional->spaceDimensionZ(); }

double Functional::fx(const DoubleVector &x) const
{
    const_cast<Functional*>(this)->x = x;

    ih->implicit_calculate_D2V1();

    return integral(uT);
}

auto Functional::integral(const DoubleMatrix &uT) const -> double
{
    const Dimension &dimensionX = spaceDimensionX();
    const Dimension &dimensionY = spaceDimensionY();
    const unsigned int N = static_cast<unsigned int> ( dimensionX.size()-1 );
    const unsigned int M = static_cast<unsigned int> ( dimensionY.size()-1 );
    const double hx = dimensionX.step();
    const double hy = dimensionY.step();

    double udiff = 0.0;
    double usum = 0.0;

    udiff = (U[0][0]-uT[0][0]); usum += 0.25 * udiff * udiff;// * mu(0, 0);
    udiff = (U[0][N]-uT[0][N]); usum += 0.25 * udiff * udiff;// * mu(N, 0);
    udiff = (U[M][0]-uT[M][0]); usum += 0.25 * udiff * udiff;// * mu(0, M);
    udiff = (U[M][N]-uT[M][N]); usum += 0.25 * udiff * udiff;// * mu(N, M);

    for (unsigned int n=1; n<=N-1; n++)
    {
        udiff = U[0][n]-uT[0][n]; usum += 0.5 * udiff * udiff;// * mu(n, 0);
        udiff = U[M][n]-uT[M][n]; usum += 0.5 * udiff * udiff;// * mu(n, M);
    }

    for (unsigned int m=1; m<=M-1; m++)
    {
        udiff = U[m][0]-uT[m][0]; usum += 0.5 * udiff * udiff;// * mu(0, m);
        udiff = U[m][N]-uT[m][N]; usum += 0.5 * udiff * udiff;// * mu(N, m);
    }

    for (unsigned int m=1; m<=M-1; m++)
    {
        for (unsigned int n=1; n<=N-1; n++)
        {
            udiff = U[m][n]-uT[m][n]; usum += udiff * udiff;// * mu(n, m);
        }
    }

    return usum*(hx*hy);
}

void Functional::gradient(const DoubleVector &v, DoubleVector &g) const
{
    const_cast<Functional*>(this)->x = x;

    puts("OK1");
    ih->implicit_calculate_D2V1();
    puts("OK2");
    fh->implicit_calculate_D2V1();
    puts("OK3");

    g.clear();
    g.resize(x.length());

    unsigned int size = timeDimension().size();

    for (unsigned int ln=0; ln<size; ln++)
    {
        for (size_t i=0; i<heat_source_number; i++)
        {
            g[i*size+ln] = -pp[i][ln];
        }
    }

}

