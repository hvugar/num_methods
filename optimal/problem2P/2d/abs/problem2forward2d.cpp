#include "problem2forward2d.h"
#include "abstractproblem22d.h"

double Problem2Forward2D::initial(const SpaceNodePDE &) const
{
    return fi;
}

double Problem2Forward2D::boundary(const SpaceNodePDE &, const TimeNodePDE &, BoundaryType) const
{
    return NAN;
}

double Problem2Forward2D::f(const SpaceNodePDE &, const TimeNodePDE &) const
{
    return 0.0;
}

double Problem2Forward2D::g1(const SpaceNodePDE &, const TimeNodePDE &) const
{
    return 0.0;
}

double Problem2Forward2D::g2(const SpaceNodePDE &, const TimeNodePDE &) const
{
    return 0.0;
}

double Problem2Forward2D::g3(const SpaceNodePDE &, const TimeNodePDE &) const
{
    return 0.0;
}

double Problem2Forward2D::g4(const SpaceNodePDE &, const TimeNodePDE &) const
{
    return 0.0;
}

double Min = +10000000.0;
double Max = -10000000.0;
void Problem2Forward2D::layerInfo(const DoubleMatrix &u, unsigned int ln) const
{
    C_UNUSED(u);
    C_UNUSED(ln);

//    if (u.min()<Min) Min = u.min();
//    if (u.max()>Max) Max = u.max();
//    printf("%f %f\n", Min, Max);

//    QPixmap px;
//    visualizeMatrixHeat(u, 0.2, 15.344541, px);
//    px.save(QString("pics/u/u%1.png").arg(ln), "PNG");

    //    QPixmap pixmap;
    //    visualizeMatrixHeat(u, u.min(), u.max(), pixmap, u.cols(), u.rows());
    //    pixmap.save(QString("image%1.png").arg(n), "PNG");
}
