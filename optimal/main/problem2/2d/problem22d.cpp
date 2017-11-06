#include "problem22d.h"

void Problem22D::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    Problem22D p22d;
    p22d.mTimeDimension = Dimension(0.1, 0, 10);
    p22d.mSpaceDimensionX = Dimension(0.01, 0, 100);
    p22d.mSpaceDimensionY = Dimension(0.01, 0, 100);

    IProblem2Forward2D ipf2d(1.0, 1.0, 1.0, 10.0, 0, 0);
    ipf2d.setTimeDimension(p22d.mTimeDimension);
    ipf2d.addSpaceDimension(p22d.mSpaceDimensionX);
    ipf2d.addSpaceDimension(p22d.mSpaceDimensionY);

    DoubleMatrix u;
    ipf2d.calculateMVD(u);
    IPrinter::printMatrix(u);
    IPrinter::printSeperatorLine();

    IProblem2Backward2D ipb2d(1.0, 1.0, 1.0, 10.0, 0, 0);
    ipb2d.setTimeDimension(p22d.mTimeDimension);
    ipb2d.addSpaceDimension(p22d.mSpaceDimensionX);
    ipb2d.addSpaceDimension(p22d.mSpaceDimensionY);

    ipb2d.U = u;
    ipb2d.uT = u;

    DoubleMatrix p;
    ipb2d.calculateMVD(p);
    IPrinter::printMatrix(p);
    IPrinter::printSeperatorLine();
}
