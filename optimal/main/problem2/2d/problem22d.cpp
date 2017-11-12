#include "problem22d.h"

void Problem22D::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    Problem22D p22d(1.0, 0.1, 1.0, 10.0, 1, 1);
    p22d.setGridParameters(Dimension(0.01, 0, 100),
                           Dimension(0.1, 0, 10),
                           Dimension(0.1, 0, 10));

    DoubleMatrix u;
    p22d.forward.calculateMVD1(u);
    //IPrinter::printMatrix(u);
    //IPrinter::printSeperatorLine();

    //    DoubleMatrix p;
    //    p22d.backward.U = u;
    //    p22d.backward.uT = u;
    //    p22d.backward.calculateMVD(p);
    //    IPrinter::printMatrix(p);
    //    IPrinter::printSeperatorLine();
}

Problem22D::Problem22D(double a, double lambda0, double lambda, double theta, double Lc, double Lo)
{
    this->a = a;
    this->lambda0 = lambda0;
    this->lambda = lambda;
    this->theta = theta;
    this->Lc = Lc;
    this->Lo = Lo;

    forward.setSettings(a, lambda0, lambda, theta, Lc, Lo);
    //backward.setSettings(a, lambda0, lambda, theta, Lc, Lo);
}

void Problem22D::setGridParameters(Dimension timeDimension, Dimension spaceDimensionX, Dimension spaceDimensionY)
{
    mTimeDimension = timeDimension;
    mSpaceDimensionX = spaceDimensionX;
    mSpaceDimensionY = spaceDimensionY;

    forward.setTimeDimension(mTimeDimension);
    forward.addSpaceDimension(mSpaceDimensionX);
    forward.addSpaceDimension(mSpaceDimensionY);

    backward.setTimeDimension(mTimeDimension);
    backward.addSpaceDimension(mSpaceDimensionX);
    backward.addSpaceDimension(mSpaceDimensionY);

    unsigned int N1 = mSpaceDimensionX.sizeN();
    unsigned int N2 = mSpaceDimensionY.sizeN();
    U.resize(N2+1, N1+1, 10.0);


//    double hx = mSpaceDimensionX.step();
//    double hy = mSpaceDimensionY.step();
//    double sum = 0.0;
//    SpaceNodePDE eta;
//    eta.x = 0.42; eta.i = 5;
//    eta.y = 0.48; eta.j = 5;

//    //int X = (int)(round(N1*eta.x));
//    //int Y = (int)(round(N2*eta.y));
//    //printf("%d %f %d %f\n", X, round(N1*mux), Y, round(N2*muy));
//    SpaceNodePDE sn;
//    for (unsigned int m=0; m<=N2; m++)
//    {
//        double y = m*hy;
//        sn.j = m; sn.y = y;
//        for (unsigned int n=0; n<=N1; n++)
//        {
//            double x = n*hx;
//            sn.i = n; sn.x = x;
//            sum += forward.delta2(sn, eta, 0);
//        }
//    }
//    printf("%f\n", sum);
//    exit(-1);
}

double Problem22D::fx(const DoubleVector &x) const
{
    DoubleMatrix u;
    forward.calculateMVD(u);

    double intgrl = integral(u);

    return intgrl;
}

void Problem22D::gradient(const DoubleVector &x UNUSED_PARAM, DoubleVector &g UNUSED_PARAM)
{}

double Problem22D::integral(const DoubleMatrix &u) const
{
    double hx = mSpaceDimensionX.step();
    double hy = mSpaceDimensionY.step();
    unsigned int N1 = mSpaceDimensionX.sizeN();
    unsigned int N2 = mSpaceDimensionY.sizeN();

    double sum = 0.0;

    sum += 0.25* (u[0][0] - U[0][0])*(u[0][0] - U[0][0]);
    sum += 0.25*(u[0][N1] - U[0][N1])*(u[0][N1] - U[0][N1]);
    sum += 0.25*(u[N2][0] - U[N2][0])*(u[N2][0] - U[N2][0]);
    sum += 0.25*(u[N2][N1] - U[N2][N1])*(u[N2][N1] - U[N2][N1]);

    for (unsigned int n1=1; n1<=N1-1; n1++)
    {
        sum += 0.5*(u[0][n1] - U[0][n1])*(u[0][n1] - U[0][n1]);
        sum += 0.5*(u[N2][n1] - U[N2][n1])*(u[N2][n1] - U[N2][n1]);
    }

    for (unsigned int n2=1; n2<=N2-1; n2++)
    {
        sum += 0.5*(u[n2][0] - U[n2][0])*(u[n2][0] - U[n2][0]);
        sum += 0.5*(u[n2][N1] - U[n2][N1])*(u[n2][N1] - U[n2][N1]);
    }

    for (unsigned int n2 = 1; n2 <= N2-1; n2++)
    {
        for (unsigned int n1 = 1; n1 <= N1-1; n1++)
        {
            sum += (u[n2][n1] - U[n2][n1])*(u[n2][n1] - U[n2][n1]);
        }
    }
    sum *= (hx*hy);

    return sum;
}

double Problem22D::mu(double x UNUSED_PARAM, double y UNUSED_PARAM) const
{
    return 1.0;
}


