#include "solver5.h"

using namespace p3p5;

void Functional::Main(int /*argc*/, char** /*argv*/)
{
    Functional f;

    LoadedHeatEquationIBVP lhe(&f);
    lhe.setThermalDiffusivity(1.0);
    lhe.setThermalConvection(-f.lambda0);
    lhe.setThermalConductivity(0.0);
    lhe.implicit_calculate_D2V1();
}

LoadedHeatEquationIBVP::LoadedHeatEquationIBVP(Functional*f) : _functional(f) {}

LoadedHeatEquationIBVP::LoadedHeatEquationIBVP(const LoadedHeatEquationIBVP &) {}

LoadedHeatEquationIBVP::~LoadedHeatEquationIBVP() {}

LoadedHeatEquationIBVP & LoadedHeatEquationIBVP::operator =(const LoadedHeatEquationIBVP &) {}

auto LoadedHeatEquationIBVP::initial(const SpaceNodePDE &/*sn*/, InitialCondition /*ic*/) const -> double
{
    const double initial_temperature = _functional->initial_temperature;
    return initial_temperature;
}

auto LoadedHeatEquationIBVP::boundary(const SpaceNodePDE &/*sn*/, const TimeNodePDE &/*tn*/, BoundaryConditionPDE &bc) const -> double
{
    //const double lambda = _functional->lambda1;
    //const double envrmnt_temperature = _functional->envrmnt_temperature;
    //bc = BoundaryConditionPDE::Dirichlet(1.0, 1.0); return 0.0;
    bc = BoundaryConditionPDE::Neumann(); return 0.0;
    //bc = BoundaryConditionPDE::Robin(lambda, +1.0); return lambda * envrmnt_temperature;
}

auto LoadedHeatEquationIBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double
{
    const size_t ln = static_cast<size_t>(tn.i);
    const double envrmnt_temperature = _functional->envrmnt_temperature;
    const double lambda0 = _functional->lambda0;
    const double dimX_step = spaceDimensionX().step();
    const double dimY_step = spaceDimensionY().step();
    const size_t heating_point_number = _functional->heating_point_number;

    //double fx = -thermalConvection() * envrmnt_temperature;
    double fx = lambda0 * envrmnt_temperature;

    double sum = 0.0;
    for (size_t i=0; i<heating_point_number; i++)
    {
        const SpacePoint &zi = _functional->z(tn, i);
        const double qi = _functional->q(tn, i);
        sum += qi * DeltaFunction::gaussian(sn, zi, SpacePoint(dimX_step, dimY_step));
        //sum += qi * DeltaFunction::nearest(sn, zi, dimX_step, dimY_step, dimX_size, dimY_size);
    }

    return fx + sum;
}

auto LoadedHeatEquationIBVP::layerInfo(const DoubleMatrix &u, const TimeNodePDE &tn) const -> void
{
    static double MIN = +10000.0;
    static double MAX = -10000.0;
    if (u.max()>MAX) MAX = u.max();
    if (u.min()<MIN) MIN = u.min();
    printf(">>> %6d %14.6f %14.6f %14.6f %14.6f\n", tn.i, u.min(), u.max(), MIN, MAX);

    //if (tn.i%10==0)
//    {
        QString filename = QString("data/problem3P/f/png/%1.png").arg(tn.i, 8, 10, QChar('0'));
        QPixmap pixmap;
        visualizeMatrixHeat(u, -5.0,  15.0, pixmap, 101, 101);
//        //visualizeMatrixHeat(u, 1.50, 2.1807, pixmap, spaceDimensionX().size(), spaceDimensionY().size());
//        //visualizeMatrixHeat(u, 0.50, 2.32, pixmap, spaceDimensionX().size(), spaceDimensionY().size());
//        visualizeMatrixHeat(u, u.min(), u.max(), pixmap, spaceDimensionX().size(), spaceDimensionY().size());

//                QPainter painter(&pixmap);
//                painter.setFont(QFont("Consolas", 12));
//                painter.setPen(Qt::white);

//                painter.drawText(2, 30, QString("time: ")+QString::number(tn.t, 'f', 3));
//                painter.drawText(2, 60, QString("temp: ")+QString::number(u.max(), 'f', 3)+QString(" max temp: ")+QString::number(MAX, 'f', 3));
        pixmap.save(filename);
//    }
}

auto LoadedHeatEquationIBVP::timeDimension() const -> Dimension { return _functional->timeDimension(); }

auto LoadedHeatEquationIBVP::spaceDimensionX() const -> Dimension { return _functional->spaceDimensionX(); }

auto LoadedHeatEquationIBVP::spaceDimensionY() const -> Dimension { return _functional->spaceDimensionY(); }

auto LoadedHeatEquationIBVP::spaceDimensionZ() const -> Dimension { return _functional->spaceDimensionZ(); }

auto Functional::q(const TimeNodePDE &tn, size_t i) const -> double
{
    return 1.0;
}

auto Functional::z(const TimeNodePDE &tn, size_t i) const -> SpacePoint
{
    SpacePoint sp;
    const double t = tn.t;
    const double v = 2.0*M_PI*t;
    switch (i) {
    case 0: {
        sp.x = 0.4*sin(v*t) + 0.5;
        sp.y = 0.4*cos(v*t) + 0.5;
    } break;
    case 1: {
        sp.x = -0.2*sin(v*t) + 0.5;
        sp.y = -0.2*cos(v*t) + 0.5;
    } break;
    }

    return sp;
}
