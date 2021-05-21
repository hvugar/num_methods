#include "solver5.h"

using namespace p3p5;

void Functional::Main(int /*argc*/, char** /*argv*/)
{
    Functional f(0.01, 0.10, 0.01);

    std::cout << f.LoadedHeatEquationIBVP::thermalDiffusivity() << " " << f.LoadedHeatEquationIBVP::thermalConvection() << " " << f.LoadedHeatEquationIBVP::thermalConductivity() << std::endl;
    std::cout << f.LoadedHeatEquationFBVP::thermalDiffusivity() << " " << f.LoadedHeatEquationFBVP::thermalConvection() << " " << f.LoadedHeatEquationFBVP::thermalConductivity() << std::endl;

    f.LoadedHeatEquationIBVP::implicit_calculate_D2V1();
    f.LoadedHeatEquationFBVP::implicit_calculate_D2V1();
}

Functional::Functional(double diffusivity, double convection, double conductivity, double lambda):
        LoadedHeatEquationIBVP(+diffusivity, -convection, conductivity),
        LoadedHeatEquationFBVP(-diffusivity, +convection, conductivity)
{
    this->lambda1 = lambda;
}

/*********************************************************************************************************************************************/

LoadedHeatEquationIBVP::LoadedHeatEquationIBVP(double a, double lambda0, double lambda1)
{
    setThermalDiffusivity(a);
    setThermalConvection(lambda0);
    setThermalConductivity(lambda1);
}

LoadedHeatEquationIBVP::LoadedHeatEquationIBVP(const LoadedHeatEquationIBVP &) {}

LoadedHeatEquationIBVP::~LoadedHeatEquationIBVP() {}

LoadedHeatEquationIBVP& LoadedHeatEquationIBVP::operator =(const LoadedHeatEquationIBVP &) {}

auto LoadedHeatEquationIBVP::initial(const SpaceNodePDE &/*sn*/, InitialCondition /*ic*/) const -> double { return initial_temperature; }

auto LoadedHeatEquationIBVP::boundary(const SpaceNodePDE &/*sn*/, const TimeNodePDE &/*tn*/, BoundaryConditionPDE &bc) const -> double
{
    //bc = BoundaryConditionPDE::Dirichlet(); return 0.0;
    //bc = BoundaryConditionPDE::Neumann(); return 0.0;
    {
        bc = BoundaryConditionPDE::Robin(lambda1, +1.0); return lambda1 * envrmnt_temperature;
    }
}

auto LoadedHeatEquationIBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double
{
    //const double envrmnt_temperature = _functional->envrmnt_temperature;
    //const double lambda0 = _functional->lambda0;
    const double dimX_step = spaceDimensionX().step();
    const double dimY_step = spaceDimensionY().step();
    //const size_t heating_point_number = _functional->heating_point_number;

    double fx = -thermalConvection() * envrmnt_temperature;

    double sum = 0.0;
    for (size_t i=0; i<heating_point_number; i++)
    {
        const SpacePoint &zi = z(tn, i);
        const double qi = q(tn, i);
        sum += qi * DeltaFunction::gaussian(sn, zi, SpacePoint(dimX_step, dimY_step));
        //sum += qi * DeltaFunction::nearest(sn, zi, dimX_step, dimY_step, dimX_size, dimY_size);
    }

    return fx + sum;
}

auto LoadedHeatEquationIBVP::layerInfo(const DoubleMatrix &u, const TimeNodePDE &tn) const -> void
{
    if (static_cast<int>(tn.i) == timeDimension().max())
    {
        LoadedHeatEquationIBVP* lhe = const_cast<LoadedHeatEquationIBVP*>(this);
        lhe->Shared::U = u;
    }

    static double MIN = +10000.0;
    static double MAX = -10000.0;
    if (u.max()>MAX) MAX = u.max();
    if (u.min()<MIN) MIN = u.min();
    printf(">>> %6d %14.6f %14.6f %14.6f %14.6f\n", tn.i, u.min(), u.max(), MIN, MAX);

    //if (tn.i%10==0)
    //    {
    QString filename = QString("data/problem3P/f/png/%1.png").arg(tn.i, 8, 10, QChar('0'));
    QPixmap pixmap;
    visualizeMatrixHeat(u, 0.0, 3.864, pixmap, 101, 101);
    //visualizeMatrixHeat(u, 1.50, 2.1807, pixmap, spaceDimensionX().size(), spaceDimensionY().size());
    //visualizeMatrixHeat(u, 0.50, 2.32, pixmap, spaceDimensionX().size(), spaceDimensionY().size());
    //visualizeMatrixHeat(u, u.min(), u.max(), pixmap, spaceDimensionX().size(), spaceDimensionY().size());

    //QPainter painter(&pixmap);
    //painter.setFont(QFont("Consolas", 12));
    //painter.setPen(Qt::white);

    //painter.drawText(2, 30, QString("time: ")+QString::number(tn.t, 'f', 3));
    //painter.drawText(2, 60, QString("temp: ")+QString::number(u.max(), 'f', 3)+QString(" max temp: ")+QString::number(MAX, 'f', 3));
    pixmap.save(filename);
    //    }
}

auto LoadedHeatEquationIBVP::timeDimension() const -> Dimension { return Shared::_timeDimension; }

auto LoadedHeatEquationIBVP::spaceDimensionX() const -> Dimension { return Shared::_spaceDimensionX; }

auto LoadedHeatEquationIBVP::spaceDimensionY() const -> Dimension { return Shared::_spaceDimensionY; }

auto LoadedHeatEquationIBVP::spaceDimensionZ() const -> Dimension { return Shared::_spaceDimensionZ; }

/*********************************************************************************************************************************************/

LoadedHeatEquationFBVP::LoadedHeatEquationFBVP(double a, double lambda0, double lambda1)
{
    setThermalDiffusivity(a);
    setThermalConvection(lambda0);
    setThermalConductivity(lambda1);
}

LoadedHeatEquationFBVP::LoadedHeatEquationFBVP(const LoadedHeatEquationFBVP &) {}

LoadedHeatEquationFBVP::~LoadedHeatEquationFBVP() {}

LoadedHeatEquationFBVP& LoadedHeatEquationFBVP::operator =(const LoadedHeatEquationFBVP &) {}

auto LoadedHeatEquationFBVP::final(const SpaceNodePDE &sn, FinalCondition /*ic*/) const -> double
{
    return -2.0 * ( U[sn.j][sn.i] /*- V[sn.j][sn.i]*/ );
}

auto LoadedHeatEquationFBVP::boundary(const SpaceNodePDE &/*sn*/, const TimeNodePDE &/*tn*/, BoundaryConditionPDE &bc) const -> double
{
    //bc = BoundaryConditionPDE::Dirichlet(); return 0.0;
    //bc = BoundaryConditionPDE::Neumann(); return 0.0;
    {
        //const double lambda = _functional->lambda1;
        bc = BoundaryConditionPDE::Robin(lambda1, +1.0); return 0.0;
    }
}

auto LoadedHeatEquationFBVP::f(const SpaceNodePDE &/*sn*/, const TimeNodePDE &tn) const -> double
{
    return 0.0;
}

auto LoadedHeatEquationFBVP::layerInfo(const DoubleMatrix &p, const TimeNodePDE &tn) const -> void
{
    static double MIN = +10000.0;
    static double MAX = -10000.0;
    if (p.max()>MAX) MAX = p.max();
    if (p.min()<MIN) MIN = p.min();
    printf(">>> %6d %14.6f %14.6f %14.6f %14.6f\n", tn.i, p.min(), p.max(), MIN, MAX);

    //if (tn.i%10==0)
    //    {
    QString filename = QString("data/problem3P/b/png/%1.png").arg(tn.i, 8, 10, QChar('0'));
    QPixmap pixmap;
    visualizeMatrixHeat(p, p.min(), p.max(), pixmap, 101, 101);
    //visualizeMatrixHeat(u, 1.50, 2.1807, pixmap, spaceDimensionX().size(), spaceDimensionY().size());
    //visualizeMatrixHeat(u, 0.50, 2.32, pixmap, spaceDimensionX().size(), spaceDimensionY().size());
    //visualizeMatrixHeat(u, u.min(), u.max(), pixmap, spaceDimensionX().size(), spaceDimensionY().size());

    //                QPainter painter(&pixmap);
    //                painter.setFont(QFont("Consolas", 12));
    //                painter.setPen(Qt::white);

    //                painter.drawText(2, 30, QString("time: ")+QString::number(tn.t, 'f', 3));
    //                painter.drawText(2, 60, QString("temp: ")+QString::number(u.max(), 'f', 3)+QString(" max temp: ")+QString::number(MAX, 'f', 3));
    pixmap.save(filename);
    //    }
}

auto LoadedHeatEquationFBVP::timeDimension() const -> Dimension { return Shared::_timeDimension; }

auto LoadedHeatEquationFBVP::spaceDimensionX() const -> Dimension { return Shared::_spaceDimensionX; }

auto LoadedHeatEquationFBVP::spaceDimensionY() const -> Dimension { return Shared::_spaceDimensionY; }

auto LoadedHeatEquationFBVP::spaceDimensionZ() const -> Dimension { return Shared::_spaceDimensionZ; }

/*********************************************************************************************************************************************/

auto Shared::q(const TimeNodePDE &tn, size_t /*i*/) const -> double
{
    return 0.05*tn.t * 0.0;
}

auto Shared::z(const TimeNodePDE &tn, size_t i) const -> SpacePoint
{
    SpacePoint sp;
    const double t = tn.t;
    const double v = 2.0*M_PI/**t*/;
    switch (i) {
    case 0: {
        sp.x = /*0.2*/0.4*sin(v*t) + 0.5;;
        sp.y = /*0.2*/0.4*cos(v*t) + 0.5;;
    } break;
    case 1: {
        sp.x = -0.2*sin(v*t) + 0.5;
        sp.y = -0.2*cos(v*t) + 0.5;
    } break;
    }

    return sp;
}
