#include "problem2backward2d.h"
#include "abstractproblem22d.h"

double Problem2Backward2D::initial(const SpaceNodePDE & sn) const
{
    return -2.0 * ap22d->mu(sn.x, sn.y) * ((*u)[sn.j][sn.i] - (*U)[sn.j][sn.i]);
}

double Problem2Backward2D::boundary(const SpaceNodePDE &, const TimeNodePDE &, BoundaryType) const
{
    return NAN;
}

double Problem2Backward2D::f(const SpaceNodePDE &, const TimeNodePDE &) const
{
    return 0.0;
}

double Problem2Backward2D::penalty(unsigned int i, const TimeNodePDE &tn) const
{
    return ap22d->gpi(i, tn.i, *info);
}

double Problem2Backward2D::g1(const SpaceNodePDE &, const TimeNodePDE &) const
{
    return 0.0;
}

double Problem2Backward2D::g2(const SpaceNodePDE &, const TimeNodePDE &) const
{
    return 0.0;
}

double Problem2Backward2D::g3(const SpaceNodePDE &, const TimeNodePDE &) const
{
    return 0.0;
}

double Problem2Backward2D::g4(const SpaceNodePDE &, const TimeNodePDE &) const
{
    return 0.0;
}

void Problem2Backward2D::layerInfo(const DoubleMatrix &p, unsigned int ln) const
{
    C_UNUSED(p);
    C_UNUSED(ln);

    //    if (ln==timeDimension().sizeN())
    //    {
    //        IPrinter::printSeperatorLine("P100");
    //        IPrinter::printMatrix(p);
    //        IPrinter::printSeperatorLine();

    //        QPixmap px;
    //        visualizeMatrixHeat(p,  p.min(), p.max(), px);
    //        px.save("p100.png", "PNG");
    //    }
    //    if (ln==0)
    //    {
    //        IPrinter::printSeperatorLine("P0");
    //        IPrinter::printMatrix(p);
    //        IPrinter::printSeperatorLine();

    //        QPixmap px;
    //        visualizeMatrixHeat(p,  p.min(), p.max(), px);
    //        px.save("p0.png", "PNG");
    //    }
}
