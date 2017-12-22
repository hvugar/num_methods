#include "problem2backward2dex4.h"
#include "abstractproblem22d.h"

double Problem2Backward2DEx4::initial(const SpaceNodePDE & sn) const
{
    return -2.0 * ap22d->mu(sn.x, sn.y) * ((*u)[sn.j][sn.i] - (*U)[sn.j][sn.i]);
}

double Problem2Backward2DEx4::boundary(const SpaceNodePDE &, const TimeNodePDE &, BoundaryType) const
{
    return NAN;
}

double Problem2Backward2DEx4::f(const SpaceNodePDE &, const TimeNodePDE &) const
{
    return 0.0;
}

double Problem2Backward2DEx4::g1(const SpaceNodePDE &, const TimeNodePDE &) const
{
    return 0.0;
}

double Problem2Backward2DEx4::g2(const SpaceNodePDE &, const TimeNodePDE &) const
{
    return 0.0;
}

double Problem2Backward2DEx4::g3(const SpaceNodePDE &, const TimeNodePDE &) const
{
    return 0.0;
}

double Problem2Backward2DEx4::g4(const SpaceNodePDE &, const TimeNodePDE &) const
{
    return 0.0;
}

void Problem2Backward2DEx4::layerInfo(const DoubleMatrix &p, unsigned int ln) const
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
