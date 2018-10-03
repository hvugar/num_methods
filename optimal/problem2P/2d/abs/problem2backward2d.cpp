#include "problem2backward2d.h"
#include "ifunctional.h"

double Problem2Backward2D::initial(const SpaceNodePDE & sn) const
{
    return -2.0 * func->mu(sn.x, sn.y) * ((*u)[sn.j][sn.i] - (*U)[sn.j][sn.i]);
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
    return func->gpi(i, tn.i, *info)*sgn(func->g0i(i, tn.i, *info));
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

//    if (p.min()<MinB) MinB = p.min();
//    if (p.max()>MaxB) MaxB = p.max();
//    printf("%f %f\n", MinB, MaxB);

//    QPixmap px;
//    visualizeMatrixHeat(p, p.min(), p.max(), px);
//    px.save(QString("pics/p/p%1.png").arg(ln), "PNG");


//        if (ln==timeDimension().sizeN())
//        {
////            IPrinter::printSeperatorLine("P100");
////            IPrinter::printMatrix(p);
////            IPrinter::printSeperatorLine();

//            QPixmap px;
//            visualizeMatrixHeat(p,  p.min(), p.max(), px);
//            px.save(QString("p%1").arg(ln), "PNG");
//        }
//        if (ln==0)
//        {
////            IPrinter::printSeperatorLine("P0");
////            IPrinter::printMatrix(p);
////            IPrinter::printSeperatorLine();

//            QPixmap px;
//            visualizeMatrixHeat(p,  p.min(), p.max(), px);
//            px.save("p0.png", "PNG");
//        }
}
