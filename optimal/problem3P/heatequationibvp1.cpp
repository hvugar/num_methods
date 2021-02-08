#include "heatequationibvp1.h"

HeatEquationIBVP1::HeatEquationIBVP1() {}

HeatEquationIBVP1::HeatEquationIBVP1(const HeatEquationIBVP1 &) {}

HeatEquationIBVP1 & HeatEquationIBVP1::operator =(const HeatEquationIBVP1 &)
{
    return *this;
}

HeatEquationIBVP1::~HeatEquationIBVP1() {}

auto HeatEquationIBVP1::initial(const SpaceNodePDE &/*sn*/, InitialCondition /*condition*/) const -> double
{
    return 0.0;
}

auto HeatEquationIBVP1::boundary(const SpaceNodePDE &/*sn*/, const TimeNodePDE &tn, BoundaryConditionPDE &cn) const -> double
{
    //cn = BoundaryConditionPDE::Dirichlet(1.0, 1.0);
    //cn = BoundaryConditionPDE::Neumann(1.0, lambda1);
    cn = BoundaryConditionPDE::Robin(lambda1, -1.0);
    return lambda1 * theta(tn);
}

auto HeatEquationIBVP1::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double
{
    auto fx = 0.0;
    fx += q(0, tn)*DeltaFunction::gaussian(sn, SpacePoint(0.25+tn.t/2.0, 0.30), SpacePoint(0.01, 0.01));
    fx += q(1, tn)*DeltaFunction::gaussian(sn, SpacePoint(0.75-tn.t/2.0, 0.70), SpacePoint(0.01, 0.01));
    return fx;
}

auto HeatEquationIBVP1::q(size_t i, const TimeNodePDE &tn) const -> double
{
    if (i==0) return tn.t;
    if (i==1) return tn.t;
    return 0.0;
}

auto HeatEquationIBVP1::theta(const TimeNodePDE &tn) const -> double { return tn.t; }

auto HeatEquationIBVP1::layerInfo(const DoubleMatrix &u, const TimeNodePDE &tn) const -> void
{
    frw_saveToImage(u, tn);
}

auto HeatEquationIBVP1::timeDimension() const -> Dimension
{
    return Dimension(0.01, 0, 100);
}

auto HeatEquationIBVP1::spaceDimensionX() const -> Dimension
{
    return Dimension(0.01, 0, 100);
}

auto HeatEquationIBVP1::spaceDimensionY() const -> Dimension
{
    return Dimension(0.01, 0, 100);
}

auto HeatEquationIBVP1::spaceDimensionZ() const -> Dimension
{
    throw new std::exception;
}

void HeatEquationIBVP1::frw_saveToImage(const DoubleMatrix &u, const TimeNodePDE &tn) const
{
    static double MIN = +10000.0;
    static double MAX = -10000.0;
    if (u.max()>MAX) MAX = u.max();
    if (u.min()<MIN) MIN = u.min();

//    const auto _min = u.min();
//    const auto _max = u.max();

    const auto _min = 0.000000;
    const auto _max = 1.461945;

    printf(">>> %6d %14.6f %14.6f %14.6f %14.6f\n", tn.i, u.min(), u.max(), MIN, MAX);

    QString filename = QString("data/problem3P/f/png/2/%1.png").arg(tn.i, 8, 10, QChar('0'));
    QPixmap pixmap;
    visualizeMatrixHeat(u, _min, _max, pixmap, spaceDimensionX().size(), spaceDimensionY().size());
    pixmap.save(filename);
}

