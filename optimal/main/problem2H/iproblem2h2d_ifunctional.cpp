#include "iproblem2h2d_ifunctional.h"

using namespace IProblem2H2D_NS;

double IFunctional::fx(const DoubleVector &prms) const
{
    IFunctional* ifunc = const_cast<IFunctional*>(this);

    //ifunc->fromVector(pv, ifunc->mParameter);
    //functional->mParameter.fromVector(prms);
    //ifunc->forward->setParameter(mParameter);
    //ifunc->backward->setParameter(mParameter);

    DoubleMatrix u;
    DoubleMatrix ut;
    //vector<ExtendedSpaceNode2D> info;
    ifunc->forward->calculateMVD(u, ut);

    double intgrl = integral(u, ut);
    //u.clear();

    return intgrl + regEpsilon*norm();// + r*penalty(info);
}

double IFunctional::integral(const DoubleMatrix &u, const DoubleMatrix ut) const
{
    double hx = mSpaceDimensionX.step();
    double hy = mSpaceDimensionY.step();
    unsigned int N1 = mSpaceDimensionX.sizeN();
    unsigned int N2 = mSpaceDimensionY.sizeN();

    double sum0 = 0.0;

    sum0 += 0.25*(u[0][0]   - U0[0][0])   * (u[0][0]   - U0[0][0]);
    sum0 += 0.25*(u[0][N1]  - U0[0][N1])  * (u[0][N1]  - U0[0][N1]);
    sum0 += 0.25*(u[N2][0]  - U0[N2][0])  * (u[N2][0]  - U0[N2][0]);
    sum0 += 0.25*(u[N2][N1] - U0[N2][N1]) * (u[N2][N1] - U0[N2][N1]);

    for (unsigned int n1=1; n1<=N1-1; n1++)
    {
        sum0 += 0.5*(u[0][n1]  - U0[0][n1]) *(u[0][n1]  - U0[0][n1]);
        sum0 += 0.5*(u[N2][n1] - U0[N2][n1])*(u[N2][n1] - U0[N2][n1]);
    }

    for (unsigned int n2=1; n2<=N2-1; n2++)
    {
        sum0 += 0.5*(u[n2][0]  - U0[n2][0]) *(u[n2][0]  - U0[n2][0]);
        sum0 += 0.5*(u[n2][N1] - U0[n2][N1])*(u[n2][N1] - U0[n2][N1]);
    }

    for (unsigned int n2 = 1; n2 <= N2-1; n2++)
    {
        for (unsigned int n1 = 1; n1 <= N1-1; n1++)
        {
            sum0 += (u[n2][n1] - U0[n2][n1])*(u[n2][n1] - U0[n2][n1]);
        }
    }
    sum0 *= (hx*hy);

    double sum1 = 0.0;

    sum1 += 0.25*(ut[0][0]   - U1[0][0])   * (ut[0][0]   - U1[0][0]);
    sum1 += 0.25*(ut[0][N1]  - U1[0][N1])  * (ut[0][N1]  - U1[0][N1]);
    sum1 += 0.25*(ut[N2][0]  - U1[N2][0])  * (ut[N2][0]  - U1[N2][0]);
    sum1 += 0.25*(ut[N2][N1] - U1[N2][N1]) * (ut[N2][N1] - U1[N2][N1]);

    for (unsigned int n1=1; n1<=N1-1; n1++)
    {
        sum1 += 0.5*(ut[0][n1]  - U1[0][n1]) *(ut[0][n1]  - U1[0][n1]);
        sum1 += 0.5*(ut[N2][n1] - U1[N2][n1])*(ut[N2][n1] - U1[N2][n1]);
    }

    for (unsigned int n2=1; n2<=N2-1; n2++)
    {
        sum1 += 0.5*(ut[n2][0]  - U1[n2][0]) *(ut[n2][0]  - U1[n2][0]);
        sum1 += 0.5*(ut[n2][N1] - U1[n2][N1])*(ut[n2][N1] - U1[n2][N1]);
    }

    for (unsigned int n2 = 1; n2 <= N2-1; n2++)
    {
        for (unsigned int n1 = 1; n1 <= N1-1; n1++)
        {
            sum1 += (ut[n2][n1] - U1[n2][n1])*(ut[n2][n1] - U1[n2][n1]);
        }
    }
    sum1 *= (hx*hy);

    return alpha0*sum0 + alpha1*sum1;
}

double IFunctional::norm() const
{
    return 0.0;
}
