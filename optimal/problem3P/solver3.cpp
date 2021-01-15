#include "solver3.h"

using namespace p3p3;

void Functional::Main(int /*argc*/, char** /*argv*/)
{
    HeatEquationIBVP h;
    h.setThermalDiffusivity(1.0);
    h.setThermalConductivity(0.0);
    h.setThermalConvection(-0.0);//-1.0 heating
    h.implicit_calculate_D2V1();
    //h.explicit_calculate_D2V1();
}

HeatEquationIBVP::HeatEquationIBVP() {}

HeatEquationIBVP::HeatEquationIBVP(const HeatEquationIBVP &other) { *this = other; }

HeatEquationIBVP& HeatEquationIBVP::operator =(const HeatEquationIBVP &) { return *this; }

/*virtual*/ HeatEquationIBVP::~HeatEquationIBVP() {}

/*virtual*/ auto HeatEquationIBVP::initial(const SpaceNodePDE &/*sn*/, InitialCondition /*ic*/) const -> double { return _initial_temperature; }

/*virtual*/ auto HeatEquationIBVP::boundary(const SpaceNodePDE &/*sn*/, const TimeNodePDE &/*tn*/, BoundaryConditionPDE &bc) const -> double
{
    //bc = BoundaryConditionPDE::Dirichlet(1.0, 1.0); return 0.0;
    //bc = BoundaryConditionPDE::Neumann(1.0, 0.0); return 0.0;
    bc = BoundaryConditionPDE::Robin(_lambda1, +1.0, _lambda1); return _enviroment_temperature;
}

/*virtual*/ auto HeatEquationIBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double
{
    //size_t ln = static_cast<size_t>(tn.i);

    double fx = -thermalConvection() * _enviroment_temperature;

    double sum = 0.0;
    for (size_t i=0; i<_heat_source_number; i++)
    {
        const SpacePoint &zi = z(tn.t, i);
        const double qi = q(tn.t, i);
        //sum += qi * DeltaFunction::gaussian(sn, zi, SpacePoint(spaceDimensionX().step(), spaceDimensionY().step()));
        sum += qi * DeltaFunction::nearest(sn, zi, spaceDimensionX().step(), spaceDimensionY().step(), spaceDimensionX().size(), spaceDimensionY().size());
    }

    return fx + sum;

    //if (sn.i==250 && sn.j == 250) { return  1.0/(0.002*0.002) * 0.1; } else { return  0.0; }
}

SpacePoint HeatEquationIBVP::z(double t, size_t i) const
{
    const double R = 0.40;
    SpacePoint sp;
    if (i==0)
    {
        sp.x = R*sin(v(t,i)) + 0.50;
        sp.y = R*cos(v(t,i)) + 0.50;

        //printf("0 %f %f %f\n", sp.x, sp.y, v(t,i));
    }
    else
    {
        sp.x = 0.5*R*sin(-v(t,i)-2) + 0.50;
        sp.y = 0.5*R*cos(-v(t,i)-2) + 0.50;

        //printf("1 %f %f %f\n", sp.x, sp.y, v(t,i));
    }
    return sp;
}

double HeatEquationIBVP::q(double /*t*/, size_t i) const
{
    if (i==0)
    {
        return 0.05;
    }
    else
    {
        return 0.05;
    }
}

double HeatEquationIBVP::v(double t, size_t /*i*/) const
{
    return 2.0*t;
}


auto HeatEquationIBVP::layerInfo(const DoubleMatrix &u, const TimeNodePDE &/*tn*/) const -> void
{
    //if (tn.i%1000!=0) return;

    IPrinter::printSeperatorLine();
    IPrinter::printMatrix(u);

    //frw_saveToImage(u, tn);
}

auto HeatEquationIBVP::timeDimension() const -> Dimension { return Dimension(0.01, 0, 1000) /*Dimension(0.0000005, 0, 20000000)*/; }

auto HeatEquationIBVP::spaceDimensionX() const -> Dimension { return Dimension(0.002, 0, 500); }

auto HeatEquationIBVP::spaceDimensionY() const -> Dimension { return Dimension(0.002, 0, 500); }

auto HeatEquationIBVP::spaceDimensionZ() const -> Dimension { return Dimension(0.002, 0, 500); }

Functional::Functional()
{

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


HeatEquationFBVP::HeatEquationFBVP() {}

HeatEquationFBVP::HeatEquationFBVP(const HeatEquationFBVP &other) { *this = other; }

HeatEquationFBVP& HeatEquationFBVP::operator =(const HeatEquationFBVP &) { return *this; }

HeatEquationFBVP::~HeatEquationFBVP() {}

auto HeatEquationFBVP::final(const SpaceNodePDE &/*sn*/, FinalCondition /*condition*/) const -> double { return 0.0; }

auto HeatEquationFBVP::boundary(const SpaceNodePDE &/*sn*/, const TimeNodePDE &/*tn*/, BoundaryConditionPDE &/*bc*/) const -> double
{
    //bc = BoundaryConditionPDE::Dirichlet(1.0, 1.0); return 0.0;
    //bc = BoundaryConditionPDE::Neumann(1.0, 0.0); return 0.0;
    //bc = BoundaryConditionPDE::Robin(_lambda1, +1.0, _lambda1); return _enviroment_temperature;
    return 0.0;
}

auto HeatEquationFBVP::f(const SpaceNodePDE &/*sn*/, const TimeNodePDE &/*tn*/) const -> double
{
    double fx = 0.0;
    return fx;
}

auto HeatEquationFBVP::layerInfo(const DoubleMatrix &/*p*/, const TimeNodePDE &/*tn*/) const -> void
{
}

auto HeatEquationFBVP::timeDimension() const -> Dimension { return Dimension(0.01, 0, 1000) /*Dimension(0.0000005, 0, 20000000)*/; }

auto HeatEquationFBVP::spaceDimensionX() const -> Dimension { return Dimension(0.002, 0, 500); }

auto HeatEquationFBVP::spaceDimensionY() const -> Dimension { return Dimension(0.002, 0, 500); }

auto HeatEquationFBVP::spaceDimensionZ() const -> Dimension { return Dimension(0.002, 0, 500); }
