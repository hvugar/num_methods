#include "iproblem2h2d_ifunctional.h"
#include "iproblem2hforward2d.h"
#include "iproblem2hbackward2d.h"

using namespace IProblem2H2D_NS;

IFunctional::IFunctional()
{

}

double IFunctional::fx(const DoubleVector &pv) const
{
    IProblem2H2D::Parameter parameter;
    parameter.a = mParameter.a;
    parameter.lambda = mParameter.lambda;
    parameter.lambda1 = mParameter.lambda1;
    parameter.Nc = mParameter.Nc;
    parameter.No = mParameter.No;
    parameter.Ns = mParameter.Ns;
    parameter.q = mParameter.q;
    parameter.theta = mParameter.theta;
    fromVector(pv, parameter);

    IProblem2HForward2D forward;
    forward.setTimeDimension(mTimeDimension);
    forward.addSpaceDimension(mSpaceDimensionX);
    forward.addSpaceDimension(mSpaceDimensionY);
    forward.mParameter = parameter;

    DoubleMatrix u;
    DoubleMatrix ut;
    vector<ExtendedSpaceNode2DH> info;
    forward.calculateMVD(u, ut, info, false);

    double intgrl = integral(u, ut);
    //u.clear();

    return intgrl;// + regEpsilon*norm();// + r*penalty(info);
}

double IFunctional::integral(const DoubleMatrix &u, const DoubleMatrix ut) const
{
    double hx = mSpaceDimensionX.step();
    double hy = mSpaceDimensionY.step();
    unsigned int N1 = mSpaceDimensionX.sizeN();
    unsigned int N2 = mSpaceDimensionY.sizeN();

    double sum0 = 0.0;

    sum0 += 0.25*(u[0][0]   - V0[0][0])   * (u[0][0]   - V0[0][0]);
    sum0 += 0.25*(u[0][N1]  - V0[0][N1])  * (u[0][N1]  - V0[0][N1]);
    sum0 += 0.25*(u[N2][0]  - V0[N2][0])  * (u[N2][0]  - V0[N2][0]);
    sum0 += 0.25*(u[N2][N1] - V0[N2][N1]) * (u[N2][N1] - V0[N2][N1]);

    for (unsigned int n1=1; n1<=N1-1; n1++)
    {
        sum0 += 0.5*(u[0][n1]  - V0[0][n1]) *(u[0][n1]  - V0[0][n1]);
        sum0 += 0.5*(u[N2][n1] - V0[N2][n1])*(u[N2][n1] - V0[N2][n1]);
    }

    for (unsigned int n2=1; n2<=N2-1; n2++)
    {
        sum0 += 0.5*(u[n2][0]  - V0[n2][0]) *(u[n2][0]  - V0[n2][0]);
        sum0 += 0.5*(u[n2][N1] - V0[n2][N1])*(u[n2][N1] - V0[n2][N1]);
    }

    for (unsigned int n2 = 1; n2 <= N2-1; n2++)
    {
        for (unsigned int n1 = 1; n1 <= N1-1; n1++)
        {
            sum0 += (u[n2][n1] - V0[n2][n1])*(u[n2][n1] - V0[n2][n1]);
        }
    }
    sum0 *= (hx*hy);

    double sum1 = 0.0;

    sum1 += 0.25*(ut[0][0]   - V1[0][0])   * (ut[0][0]   - V1[0][0]);
    sum1 += 0.25*(ut[0][N1]  - V1[0][N1])  * (ut[0][N1]  - V1[0][N1]);
    sum1 += 0.25*(ut[N2][0]  - V1[N2][0])  * (ut[N2][0]  - V1[N2][0]);
    sum1 += 0.25*(ut[N2][N1] - V1[N2][N1]) * (ut[N2][N1] - V1[N2][N1]);

    for (unsigned int n1=1; n1<=N1-1; n1++)
    {
        sum1 += 0.5*(ut[0][n1]  - V1[0][n1]) *(ut[0][n1]  - V1[0][n1]);
        sum1 += 0.5*(ut[N2][n1] - V1[N2][n1])*(ut[N2][n1] - V1[N2][n1]);
    }

    for (unsigned int n2=1; n2<=N2-1; n2++)
    {
        sum1 += 0.5*(ut[n2][0]  - V1[n2][0]) *(ut[n2][0]  - V1[n2][0]);
        sum1 += 0.5*(ut[n2][N1] - V1[n2][N1])*(ut[n2][N1] - V1[n2][N1]);
    }

    for (unsigned int n2 = 1; n2 <= N2-1; n2++)
    {
        for (unsigned int n1 = 1; n1 <= N1-1; n1++)
        {
            sum1 += (ut[n2][n1] - V1[n2][n1])*(ut[n2][n1] - V1[n2][n1]);
        }
    }
    sum1 *= (hx*hy);

    return alpha0*sum0 + alpha1*sum1;
}

double IFunctional::norm() const
{
    double norm = 0.0;

    for (unsigned int i=0; i<mParameter.Nc; i++)
    {
        norm += (mParameter.eta[i].x - mParameter0.eta[i].x)*(mParameter.eta[i].x - mParameter0.eta[i].x);
        norm += (mParameter.eta[i].y - mParameter0.eta[i].y)*(mParameter.eta[i].y - mParameter0.eta[i].y);

        for (unsigned int j=0; j<mParameter.No; j++)
        {
            norm += (mParameter.k[i][j] - mParameter0.k[i][j])*(mParameter.k[i][j] - mParameter0.k[i][j]);
            norm += (mParameter.z[i][j] - mParameter0.z[i][j])*(mParameter.z[i][j] - mParameter0.z[i][j]);

            norm += (mParameter.xi[j].x - mParameter0.xi[j].x)*(mParameter.xi[j].x - mParameter0.xi[j].x);
            norm += (mParameter.xi[j].y - mParameter0.xi[j].y)*(mParameter.xi[j].y - mParameter0.xi[j].y);
        }
    }

    return norm;
}

void IFunctional::gradient(const DoubleVector &pv, DoubleVector &g) const
{
    unsigned int L = mTimeDimension.sizeN();
    double ht = mTimeDimension.step();

    IProblem2H2D::Parameter parameter;
    parameter.a = mParameter.a;
    parameter.lambda = mParameter.lambda;
    parameter.lambda1 = mParameter.lambda1;
    parameter.Nc = mParameter.Nc;
    parameter.No = mParameter.No;
    parameter.Ns = mParameter.Ns;
    parameter.q = mParameter.q;
    parameter.theta = mParameter.theta;
    fromVector(pv, parameter);

    IProblem2HForward2D forward;
    forward.setTimeDimension(mTimeDimension);
    forward.addSpaceDimension(mSpaceDimensionX);
    forward.addSpaceDimension(mSpaceDimensionY);
    forward.mParameter = parameter;

    IProblem2HBackward2D backward;
    backward.setTimeDimension(mTimeDimension);
    backward.addSpaceDimension(mSpaceDimensionX);
    backward.addSpaceDimension(mSpaceDimensionY);
    backward.mParameter = parameter;
    backward.ifunc = const_cast<IProblem2H2D_NS::IFunctional*>(this);

    DoubleMatrix u;
    DoubleMatrix ut;
    DoubleMatrix p;

    vector<ExtendedSpaceNode2DH> u_info;
    forward.calculateMVD(u, ut, u_info, true);

    backward.UT = u;
    backward.UTt = ut;
    vector<ExtendedSpaceNode2DH> p_info;
    backward.calculateMVD(p, p_info, true);

    g.clear();
    g.resize(pv.length(), 0.0);
    unsigned int gi = 0;

    // k
    if (optimizeK)
    {
        for (unsigned int i=0; i<mParameter.Nc; i++)
        {
            ExtendedSpaceNode2DH &pi = p_info[i];

            for (unsigned int j=0; j<mParameter.No; j++)
            {
                ExtendedSpaceNode2DH &uj = u_info[j];

                double grad_Kij = 0.0;

                grad_Kij += 0.5 * (pi.value(0)/*+2.0*r*gpi(i,0,u_info)*sgn(g0i(i,0,u_info))*/) * (uj.value(0) - mParameter.z.at(i,j));
                for (unsigned int m=1; m<=L-1; m++)
                {
                    grad_Kij += (pi.value(m)/*+2.0*r*gpi(i,m,u_info)*sgn(g0i(i,m,u_info))*/) * (uj.value(m) - mParameter.z.at(i,j));
                }
                grad_Kij += 0.5 * (pi.value(L)/*+2.0*r*gpi(i,L,u_info)*sgn(g0i(i,L,u_info))*/) * (uj.value(L) - mParameter.z.at(i,j));
                grad_Kij *= -ht;

                g[gi++] = grad_Kij;// + 2.0*regEpsilon*(mParameter.k.at(i,j) - mParameter0.k.at(i,j));
            }
        }
    }
    else
    {
        for (unsigned int i=0; i<mParameter.Nc; i++)
        {
            for (unsigned int j=0; j<mParameter.No; j++)
            {
                g[gi++] = 0.0;
            }
        }
    }

    // z
    if (optimizeZ)
    {
        for (unsigned int i=0; i<mParameter.Nc; i++)
        {
            ExtendedSpaceNode2DH &pi = p_info[i];

            for (unsigned int j=0; j<mParameter.No; j++)
            {
                double grad_Zij = 0.0;

                grad_Zij += 0.5 * (pi.value(0)/*+2.0*r*gpi(i,0,u_info)*sgn(g0i(i,0,u_info))*/) * mParameter.k.at(i,j);
                for (unsigned int m=1; m<=L-1; m++)
                {
                    grad_Zij += (pi.value(m)/*+2.0*r*gpi(i,m,u_info)*sgn(g0i(i,m,u_info))*/)  * mParameter.k.at(i,j);
                }
                grad_Zij += 0.5 * (pi.value(L)/*+2.0*r*gpi(i,L,u_info)*sgn(g0i(i,L,u_info))*/) * mParameter.k.at(i,j);
                grad_Zij *= ht;

                g[gi++] = grad_Zij;// + 2.0*regEpsilon*(mParameter.z[i][j] - mParameter0.z[i][j]);
            }
        }
    }
    else
    {
        for (unsigned int i=0; i<mParameter.Nc; i++)
        {
            for (unsigned int j=0; j<mParameter.No; j++)
            {
                g[gi++] = 0.0;
            }
        }
    }

    // eta
    if (optimizeC)
    {
        for (unsigned int i=0; i<mParameter.Nc; i++)
        {
            ExtendedSpaceNode2DH &pi = p_info[i];

            double grad_EtaiX = 0.0;
            double grad_EtaiY = 0.0;
            double vi = 0.0;

            vi = 0.0;
            for (unsigned int j=0; j<mParameter.No; j++) vi += mParameter.k.at(i,j) * (u_info[j].value(0) - mParameter.z.at(i,j));
            grad_EtaiX += 0.5 * pi.valueDx(0) * vi;
            grad_EtaiY += 0.5 * pi.valueDy(0) * vi;

            for (unsigned int m=1; m<=L-1; m++)
            {
                vi = 0.0;
                for (unsigned int j=0; j<mParameter.No; j++) vi += mParameter.k.at(i,j) * (u_info[j].value(m) - mParameter.z.at(i,j));
                grad_EtaiX += pi.valueDx(m) * vi;
                grad_EtaiY += pi.valueDy(m) * vi;
            }

            vi = 0.0;
            for (unsigned int j=0; j<mParameter.No; j++) vi += mParameter.k.at(i,j) * (u_info[j].value(L) - mParameter.z.at(i,j));
            grad_EtaiX += 0.5 * pi.valueDx(L) * vi;
            grad_EtaiY += 0.5 * pi.valueDy(L) * vi;

            grad_EtaiX *= -ht;
            grad_EtaiY *= -ht;

            g[gi++] = grad_EtaiX;// + 2.0*regEpsilon*(mParameter.eta[i].x - mParameter0.eta[i].x);
            g[gi++] = grad_EtaiY;// + 2.0*regEpsilon*(mParameter.eta[i].y - mParameter0.eta[i].y);
        }
    }
    else
    {
        for (unsigned int i=0; i<mParameter.Nc; i++)
        {
            g[gi++] = 0.0;
            g[gi++] = 0.0;
        }
    }

    // xi
    if (optimizeO)
    {
        for (unsigned int j=0; j<mParameter.Nc; j++)
        {
            ExtendedSpaceNode2DH &uj = u_info[j];

            double gradXijX = 0.0;
            double gradXijY = 0.0;
            double vi = 0.0;

            vi = 0.0;
            for (unsigned int i=0; i<mParameter.Nc; i++)
                vi += mParameter.k.at(i,j) * (p_info[i].value(0)/*+2.0*r*gpi(i,0,u_info)*sgn(g0i(i,0,u_info))*/);
            gradXijX += 0.5 * uj.valueDx(0) * vi;
            gradXijY += 0.5 * uj.valueDy(0) * vi;

            for (unsigned int m=1; m<=L-1; m++)
            {
                vi = 0.0;
                for (unsigned int i=0; i<mParameter.Nc; i++)
                    vi += mParameter.k.at(i,j)*(p_info[i].value(m)/*+2.0*r*gpi(i,m,u_info)*sgn(g0i(i,m,u_info))*/);
                gradXijX += uj.valueDx(m) * vi;
                gradXijY += uj.valueDy(m) * vi;
            }

            vi = 0.0;
            for (unsigned int i=0; i<mParameter.Nc; i++)
                vi += mParameter.k.at(i,j)*(p_info[i].value(L)/*+2.0*r*gpi(i,L,u_info)*sgn(g0i(i,L,u_info))*/);
            gradXijX += 0.5 * uj.valueDx(L) * vi;
            gradXijY += 0.5 * uj.valueDy(L) * vi;

            gradXijX *= -ht;
            gradXijY *= -ht;

            g[gi++] = gradXijX;// + 2.0*regEpsilon*(mParameter.xi[j].x - mParameter0.xi[j].x);
            g[gi++] = gradXijY;// + 2.0*regEpsilon*(mParameter.xi[j].x - mParameter0.xi[j].x);
        }
    }
    else
    {
        for (unsigned int j=0; j<mParameter.No; j++)
        {
            g[gi++] = 0.0;
            g[gi++] = 0.0;
        }
    }

    u_info.clear();
    p_info.clear();
}

void IFunctional::toVector(const IProblem2H2D::Parameter &prm, DoubleVector &pv) const
{
    unsigned int Nc = prm.Nc;
    unsigned int No = prm.No;

    //unsigned int length = 0;
    //if (optimizeK) length += Lc*Lo;
    //if (optimizeZ) length += Lc*Lo;
    //if (optimizeC) length += 2*Lc;
    //if (optimizeO) length += 2*Lo;

    pv.clear();
    pv.resize(2*Nc*No+2*Nc+2*No);

    //if (optimizeK)
    {
        for (unsigned int i=0; i<Nc; i++)
        {
            for (unsigned int j=0; j<No; j++)
            {
                pv[i*No + j] = prm.k[i][j];
            }
        }
    }

    //if (optimizeZ)
    {
        for (unsigned int i=0; i<Nc; i++)
        {
            for (unsigned int j=0; j<No; j++)
            {
                pv[i*No + j + Nc*No] = prm.z[i][j];
            }
        }
    }

    //if (optimizeC)
    {
        for (unsigned int i=0; i<Nc; i++)
        {
            pv[2*i + 0 + 2*Nc*No] = prm.eta[i].x;
            pv[2*i + 1 + 2*Nc*No] = prm.eta[i].y;
        }
    }

    //if (optimizeO)
    {
        for (unsigned int j=0; j<No; j++)
        {
            pv[2*j + 0 + 2*Nc + 2*Nc*No] = prm.xi[j].x;
            pv[2*j + 1 + 2*Nc + 2*Nc*No] = prm.xi[j].y;
        }
    }

    //    for (unsigned int i=0; i<Lc; i++)
    //    {
    //        x[2*i + 0 + 2*Lc*Lo] = prm.eta[i].x;
    //        x[2*i + 1 + 2*Lc*Lo] = prm.eta[i].y;

    //        for (unsigned int j=0; j<Lo; j++)
    //        {
    //            x[i*Lo + j] = prm.k[i][j];
    //            x[i*Lo + j + Lc*Lo] = prm.z[i][j];

    //            x[2*j + 0 + 2*Lc + 2*Lc*Lo] = prm.xi[j].x;
    //            x[2*j + 1 + 2*Lc + 2*Lc*Lo] = prm.xi[j].y;
    //        }
    //    }
}

void IFunctional::fromVector(const DoubleVector &pv, IProblem2H2D::Parameter &prm) const
{
    unsigned int Nc = prm.Nc;
    unsigned int No = prm.No;

    unsigned int index = 0;

    //if (optimizeK)
    {
        prm.k.clear();
        prm.k.resize(Nc, No);

        for (unsigned int i=0; i<Nc; i++)
        {
            for (unsigned int j=0; j<No; j++)
            {
                prm.k[i][j] = pv[index]; index++;
            }
        }
    }

    //if (optimizeZ)
    {
        prm.z.clear();
        prm.z.resize(Nc, No);

        for (unsigned int i=0; i<Nc; i++)
        {
            for (unsigned int j=0; j<No; j++)
            {
                prm.z[i][j] = pv[index]; index++;
            }
        }
    }

    //if (optimizeC)
    {
        prm.eta.clear();
        prm.eta.resize(Nc);

        for (unsigned int i=0; i<Nc; i++)
        {
            prm.eta[i].x = pv[index]; index++;
            prm.eta[i].y = pv[index]; index++;
        }
    }

    //if (optimizeO)
    {
        prm.xi.clear();
        prm.xi.resize(No);

        for (unsigned int j=0; j<No; j++)
        {
            prm.xi[j].x = pv[index]; index++;
            prm.xi[j].y = pv[index]; index++;
        }
    }

    //    prm.k.clear();
    //    prm.k.resize(Lc, Lo);

    //    prm.z.clear();
    //    prm.z.resize(Lc, Lo);

    //    prm.eta.clear();
    //    prm.eta.resize(Lc);

    //    prm.xi.clear();
    //    prm.xi.resize(Lo);

    //    for (unsigned int i=0; i<Lc; i++)
    //    {
    //        prm.eta[i].x = pv[2*Lc*Lo + 2*i+0];
    //        prm.eta[i].y = pv[2*Lc*Lo + 2*i+1];

    //        for (unsigned int j=0; j<Lo; j++)
    //        {
    //            prm.k[i][j] = pv[i*Lo + j];
    //            prm.z[i][j] = pv[i*Lo + j + Lc*Lo];

    //            prm.xi[j].x = pv[2*Lc*Lo + 2*Lc + 2*j+0];
    //            prm.xi[j].y = pv[2*Lc*Lo + 2*Lc + 2*j+1];
    //        }
    //    }
}

