#include "iproblem2h2d_ifunctional.h"
#include "iproblem2hforward2d.h"
#include "iproblem2hbackward2d.h"

using namespace IProblem2H;

IFunctional::IFunctional()
{
    optimizeK = optimizeZ = optimizeO = optimizeC = true;
}

double IFunctional::fx(const DoubleVector &pv) const
{
    OptimizeParameter o_prm;
    fromVector(pv, o_prm);

    IProblem2HForward2D forward;
    forward.setTimeDimension(mTimeDimension);
    forward.addSpaceDimension(mSpaceDimensionX);
    forward.addSpaceDimension(mSpaceDimensionY);
    forward.mEquParameter = mEquParameter;
    forward.mOptParameter = o_prm;

    DoubleMatrix u;
    DoubleMatrix ut;
    vector<SpacePointInfo> info;
    forward.calculateMVD(u, ut, info, true);

    double intgrl = integral(u, ut);
    u.clear();
    ut.clear();

    double sum = intgrl + regEpsilon*norm(mEquParameter, o_prm, mOptParameter0) + r*penalty(info, o_prm);

    for (unsigned int i=0; i<info.size(); i++)
    {
        info[i].clearWeights();
    }

    return sum;
}

double IFunctional::integral(const DoubleMatrix &u, const DoubleMatrix &ut) const
{
    double hx = mSpaceDimensionX.step();
    double hy = mSpaceDimensionY.step();
    unsigned int N1 = mSpaceDimensionX.sizeN();
    unsigned int N2 = mSpaceDimensionY.sizeN();

    double sum0 = 0.0;
    double sum1 = 0.0;

    sum0 += 0.25*(u[0][0]   - V0[0][0])   * (u[0][0]   - V0[0][0]);
    sum0 += 0.25*(u[0][N1]  - V0[0][N1])  * (u[0][N1]  - V0[0][N1]);
    sum0 += 0.25*(u[N2][0]  - V0[N2][0])  * (u[N2][0]  - V0[N2][0]);
    sum0 += 0.25*(u[N2][N1] - V0[N2][N1]) * (u[N2][N1] - V0[N2][N1]);

    sum1 += 0.25*(ut[0][0]   - V1[0][0])   * (ut[0][0]   - V1[0][0]);
    sum1 += 0.25*(ut[0][N1]  - V1[0][N1])  * (ut[0][N1]  - V1[0][N1]);
    sum1 += 0.25*(ut[N2][0]  - V1[N2][0])  * (ut[N2][0]  - V1[N2][0]);
    sum1 += 0.25*(ut[N2][N1] - V1[N2][N1]) * (ut[N2][N1] - V1[N2][N1]);

    for (unsigned int n1=1; n1<=N1-1; n1++)
    {
        sum0 += 0.5*(u[0][n1]  - V0[0][n1]) *(u[0][n1]  - V0[0][n1]);
        sum0 += 0.5*(u[N2][n1] - V0[N2][n1])*(u[N2][n1] - V0[N2][n1]);
        sum1 += 0.5*(ut[0][n1]  - V1[0][n1]) *(ut[0][n1]  - V1[0][n1]);
        sum1 += 0.5*(ut[N2][n1] - V1[N2][n1])*(ut[N2][n1] - V1[N2][n1]);
    }

    for (unsigned int n2=1; n2<=N2-1; n2++)
    {
        sum0 += 0.5*(u[n2][0]  - V0[n2][0]) *(u[n2][0]  - V0[n2][0]);
        sum0 += 0.5*(u[n2][N1] - V0[n2][N1])*(u[n2][N1] - V0[n2][N1]);
        sum1 += 0.5*(ut[n2][0]  - V1[n2][0]) *(ut[n2][0]  - V1[n2][0]);
        sum1 += 0.5*(ut[n2][N1] - V1[n2][N1])*(ut[n2][N1] - V1[n2][N1]);
    }

    for (unsigned int n2=1; n2<=N2-1; n2++)
    {
        for (unsigned int n1=1; n1<=N1-1; n1++)
        {
            sum0 += (u[n2][n1] - V0[n2][n1])*(u[n2][n1] - V0[n2][n1]);
            sum1 += (ut[n2][n1] - V1[n2][n1])*(ut[n2][n1] - V1[n2][n1]);
        }
    }

    sum0 *= (hx*hy);
    sum1 *= (hx*hy);

    return alpha0*sum0 + alpha1*sum1;
}

double IFunctional::integral1(const DoubleMatrix &u, const DoubleMatrix &) const
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

    return sum0;
}

double IFunctional::integral2(const DoubleMatrix &, const DoubleMatrix &ut) const
{
    double hx = mSpaceDimensionX.step();
    double hy = mSpaceDimensionY.step();
    unsigned int N1 = mSpaceDimensionX.sizeN();
    unsigned int N2 = mSpaceDimensionY.sizeN();

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

    return sum1;
}

double IFunctional::norm(const EquationParameter& e_prm, const OptimizeParameter &o_prm, const OptimizeParameter &o_prm0) const
{
    double norm = 0.0;

    for (unsigned int i=0; i<e_prm.Nc; i++)
    {
        norm += (o_prm.eta[i].x - o_prm0.eta[i].x)*(o_prm.eta[i].x - o_prm0.eta[i].x);
        norm += (o_prm.eta[i].y - o_prm0.eta[i].y)*(o_prm.eta[i].y - o_prm0.eta[i].y);

        for (unsigned int j=0; j<e_prm.No; j++)
        {
            norm += (o_prm.k[i][j] - o_prm0.k[i][j])*(o_prm.k[i][j] - o_prm0.k[i][j]);
            norm += (o_prm.z[i][j] - o_prm0.z[i][j])*(o_prm.z[i][j] - o_prm0.z[i][j]);

            norm += (o_prm.xi[j].x - o_prm0.xi[j].x)*(o_prm.xi[j].x - o_prm0.xi[j].x);
            norm += (o_prm.xi[j].y - o_prm0.xi[j].y)*(o_prm.xi[j].y - o_prm0.xi[j].y);
        }
    }

    return norm;
}

double IFunctional::penalty(const vector<SpacePointInfo> &info, const OptimizeParameter &o_prm) const
{
    double ht = mTimeDimension.step();
    unsigned int L = mTimeDimension.sizeN();

    double p_sum = 0.0;

    for (unsigned int i=0; i<mEquParameter.Nc; i++)
    {
        double _gp0 = gpi(i, 0, info, o_prm);
        p_sum += 0.5*_gp0*_gp0;
        for (unsigned int l=1; l<=L-1; l++)
        {
            double _gpi = gpi(i, l, info, o_prm);
            p_sum += _gpi*_gpi;
        }
        double _gpL = gpi(i, L, info, o_prm);
        p_sum += 0.5*_gpL*_gpL;
    }

    return p_sum*ht;
}

double IFunctional::gpi(unsigned int i, unsigned int layer, const vector<SpacePointInfo> &info, const OptimizeParameter &o_prm) const
{
    double p = fabs(g0i(i, layer, info, o_prm)) - (vmax.at(i) - vmin.at(i))/2.0;
    return p > 0.0 ? p : 0.0;
}

double IFunctional::g0i(unsigned int i, unsigned int layer, const vector<SpacePointInfo> &info, const OptimizeParameter &o_prm) const
{
    double vi = 0.0;
    for (unsigned int j=0; j<mEquParameter.No; j++)
    {
        const SpacePointInfo &node = info[j];
        vi += o_prm.k[i][j] * (node.value(layer)-o_prm.z[i][j]);
    }
    return (vmax.at(i) + vmin.at(i))/2.0 - vi;
}

void IFunctional::gradient(const DoubleVector &pv, DoubleVector &g) const
{
    unsigned int L = mTimeDimension.sizeN();
    double ht = mTimeDimension.step();
    unsigned int Nc = mEquParameter.Nc;
    unsigned int No = mEquParameter.No;

    OptimizeParameter o_prm;
    fromVector(pv, o_prm);

    IProblem2HForward2D forward;
    forward.setTimeDimension(mTimeDimension);
    forward.addSpaceDimension(mSpaceDimensionX);
    forward.addSpaceDimension(mSpaceDimensionY);
    forward.mEquParameter = mEquParameter;
    forward.mOptParameter = o_prm;

    IProblem2HBackward2D backward;
    backward.setTimeDimension(mTimeDimension);
    backward.addSpaceDimension(mSpaceDimensionX);
    backward.addSpaceDimension(mSpaceDimensionY);
    backward.mEquParameter = mEquParameter;
    backward.mOptParameter = o_prm;
    backward.ifunc = const_cast<IFunctional*>(this);

    DoubleMatrix u;
    DoubleMatrix ut;
    DoubleMatrix p;

    vector<SpacePointInfo> u_info;
    forward.calculateMVD(u, ut, u_info, true);

    backward.UT = u;
    backward.UTt = ut;
    vector<SpacePointInfo> p_info;
    backward.calculateMVD(p, p_info, true, u_info);

    g.clear();
    g.resize(pv.length(), 0.0);
    unsigned int gi = 0;

    // k
    if (optimizeK)
    {
        for (unsigned int i=0; i<Nc; i++)
        {
            const SpacePointInfo &pi = p_info[i];

            for (unsigned int j=0; j<No; j++)
            {
                const SpacePointInfo &uj = u_info[j];

                double grad_Kij = 0.0;
                double zij = o_prm.z[i][j];

                grad_Kij += 0.5 * (pi.value(0) + 2.0*r*gpi(i,0,u_info,o_prm)*sgn(g0i(i,0,u_info,o_prm))) * (uj.value(0) - zij);
                for (unsigned int m=1; m<=L-1; m++)
                {
                    grad_Kij += (pi.value(m) + 2.0*r*gpi(i,m,u_info,o_prm)*sgn(g0i(i,m,u_info,o_prm))) * (uj.value(m) - zij);
                }
                grad_Kij += 0.5 * (pi.value(L) + 2.0*r*gpi(i,L,u_info,o_prm)*sgn(g0i(i,L,u_info,o_prm))) * (uj.value(L) - zij);
                grad_Kij *= -ht;

                g[gi++] = grad_Kij + 2.0*regEpsilon*(o_prm.k[i][j] - mOptParameter0.k[i][j]);
            }
        }
    }
    else
    {
        for (unsigned int i=0; i<Nc; i++)
        {
            for (unsigned int j=0; j<No; j++)
            {
                g[gi++] = 0.0;
            }
        }
    }

    // z
    if (optimizeZ)
    {
        for (unsigned int i=0; i<Nc; i++)
        {
            const SpacePointInfo &pi = p_info[i];

            for (unsigned int j=0; j<No; j++)
            {
                double grad_Zij = 0.0;
                double kij = o_prm.k[i][j];

                grad_Zij += 0.5 * (pi.value(0) + 2.0*r*gpi(i,0,u_info,o_prm)*sgn(g0i(i,0,u_info,o_prm))) * kij;
                for (unsigned int m=1; m<=L-1; m++)
                {
                    grad_Zij += (pi.value(m) + 2.0*r*gpi(i,m,u_info,o_prm)*sgn(g0i(i,m,u_info,o_prm))) * kij;
                }
                grad_Zij += 0.5 * (pi.value(L) + 2.0*r*gpi(i,L,u_info,o_prm)*sgn(g0i(i,L,u_info,o_prm))) * kij;
                grad_Zij *= ht;

                g[gi++] = grad_Zij + 2.0*regEpsilon*(o_prm.z[i][j] - mOptParameter0.z[i][j]);
            }
        }
    }
    else
    {
        for (unsigned int i=0; i<Nc; i++)
        {
            for (unsigned int j=0; j<No; j++)
            {
                g[gi++] = 0.0;
            }
        }
    }

    // xi
    if (optimizeO)
    {
        for (unsigned int j=0; j<No; j++)
        {
            const SpacePointInfo &uj = u_info[j];

            double gradXijX = 0.0;
            double gradXijY = 0.0;
            double vi = 0.0;

            vi = 0.0;
            for (unsigned int i=0; i<Nc; i++) vi += o_prm.k[i][j] * (p_info[i].value(0) + 2.0*r*gpi(i,0,u_info,o_prm)*sgn(g0i(i,0,u_info,o_prm)));
            gradXijX += 0.5 * uj.valueDx(0) * vi;
            gradXijY += 0.5 * uj.valueDy(0) * vi;

            for (unsigned int m=1; m<=L-1; m++)
            {
                vi = 0.0;
                for (unsigned int i=0; i<Nc; i++) vi += o_prm.k[i][j]*(p_info[i].value(m) + 2.0*r*gpi(i,m,u_info,o_prm)*sgn(g0i(i,m,u_info,o_prm)));
                gradXijX += uj.valueDx(m) * vi;
                gradXijY += uj.valueDy(m) * vi;
            }

            vi = 0.0;
            for (unsigned int i=0; i<Nc; i++) vi += o_prm.k[i][j]*(p_info[i].value(L) + 2.0*r*gpi(i,L,u_info,o_prm)*sgn(g0i(i,L,u_info,o_prm)));
            gradXijX += 0.5 * uj.valueDx(L) * vi;
            gradXijY += 0.5 * uj.valueDy(L) * vi;

            gradXijX *= -ht;
            gradXijY *= -ht;

            g[gi++] = gradXijX + 2.0*regEpsilon*(o_prm.xi[j].x - mOptParameter0.xi[j].x);
            g[gi++] = gradXijY + 2.0*regEpsilon*(o_prm.xi[j].y - mOptParameter0.xi[j].y);
        }
    }
    else
    {
        for (unsigned int j=0; j<No; j++)
        {
            g[gi++] = 0.0;
            g[gi++] = 0.0;
        }
    }

    // eta
    if (optimizeC)
    {
        for (unsigned int i=0; i<Nc; i++)
        {
            const SpacePointInfo &pi = p_info[i];

            double gradEtaiX = 0.0;
            double gradEtaiY = 0.0;
            double vi = 0.0;

            vi = 0.0;
            for (unsigned int j=0; j<No; j++) vi += o_prm.k[i][j] * (u_info[j].value(0) - o_prm.z[i][j]);
            gradEtaiX += 0.5 * pi.valueDx(0) * vi;
            gradEtaiY += 0.5 * pi.valueDy(0) * vi;

            for (unsigned int m=1; m<=L-1; m++)
            {
                vi = 0.0;
                for (unsigned int j=0; j<No; j++) vi += o_prm.k[i][j] * (u_info[j].value(m) - o_prm.z[i][j]);
                gradEtaiX += pi.valueDx(m) * vi;
                gradEtaiY += pi.valueDy(m) * vi;
            }

            vi = 0.0;
            for (unsigned int j=0; j<No; j++) vi += o_prm.k[i][j] * (u_info[j].value(L) - o_prm.z[i][j]);
            gradEtaiX += 0.5 * pi.valueDx(L) * vi;
            gradEtaiY += 0.5 * pi.valueDy(L) * vi;

            gradEtaiX *= -ht;
            gradEtaiY *= -ht;

            g[gi++] = gradEtaiX + 2.0*regEpsilon*(o_prm.eta[i].x - mOptParameter0.eta[i].x);
            g[gi++] = gradEtaiY + 2.0*regEpsilon*(o_prm.eta[i].y - mOptParameter0.eta[i].y);
        }
    }
    else
    {
        for (unsigned int i=0; i<Nc; i++)
        {
            g[gi++] = 0.0;
            g[gi++] = 0.0;
        }
    }

    for (unsigned int i=0; i<u_info.size(); i++)
    {
        u_info[i].clearWeights();
    }

    for (unsigned int i=0; i<p_info.size(); i++)
    {
        p_info[i].clearWeights();
    }

    u_info.clear();
    p_info.clear();

    DoubleVector gk = g.mid(0, Nc*No-1);                         gk.L2Normalize(); for (unsigned int i=0; i<Nc*No; i++) g[i] = gk[i];
    DoubleVector gz = g.mid(Nc*No, 2*Nc*No-1);                   gz.L2Normalize(); for (unsigned int i=0; i<Nc*No; i++) g[Nc*No+i] = gz[i];
    DoubleVector gx = g.mid(2*Nc*No, 2*Nc*No+2*No-1);            gx.L2Normalize(); for (unsigned int i=0; i<2*No;  i++) g[2*Nc*No+i] = gx[i];
    DoubleVector ge = g.mid(2*Nc*No+2*No, 2*Nc*No+2*No+2*Nc-1);  ge.L2Normalize(); for (unsigned int i=0; i<2*Nc; i++)  g[2*Nc*No+2*No+i] = ge[i];
}

void IFunctional::project(DoubleVector &pv, unsigned int index)
{
    unsigned int Nc = mEquParameter.Nc;
    unsigned int No = mEquParameter.No;

    unsigned int offset = 2*Nc*No;

    // xi
    if ( offset <= index && index <= offset + 2*No - 1 )
    {
        if (pv[index] < 0.05) pv[index] = 0.05;
        if (pv[index] > 0.95) pv[index] = 0.95;
    }

    // eta
    if ( offset + 2*No <= index && index <= offset + 2*No + 2*Nc - 1 )
    {
        if (pv[index] < 0.05) pv[index] = 0.05;
        if (pv[index] > 0.95) pv[index] = 0.95;
    }

//    double dx = 0.09;

//    if (index == 12 && fabs(pv[8] - pv[12])<dx)
//    {
//        pv[12] = pv[8] + dx;
//        if (pv[12] > 0.95) pv[12] = pv[8] - dx;
//    }

//    if (index == 12 && fabs(pv[10] - pv[12])<dx)
//    {
//        pv[12] = pv[10] + dx;
//        if (pv[12] > 0.95) pv[12] = pv[10] - dx;
//    }

//    if (index == 14 && fabs(pv[8] - pv[14])<dx)
//    {
//        pv[14] = pv[8] + dx;
//        if (pv[14] > 0.95) pv[14] = pv[8] - dx;
//    }
//    if (index == 14 && fabs(pv[10] - pv[14])<dx)
//    {
//        pv[14] = pv[10] + dx;
//        if (pv[14] > 0.95) pv[14] = pv[10] - dx;
//    }

//    if (index == 13 && fabs(pv[9] - pv[13])<dx)
//    {
//        pv[13] = pv[9] + dx;
//        if (pv[13] > 0.95) pv[13] = pv[9] - dx;
//    }
//    if (index == 13 && fabs(pv[11] - pv[13])<dx)
//    {
//        pv[13] = pv[11] + dx;
//        if (pv[13] > 0.95) pv[13] = pv[11] - dx;
//    }

//    if (index == 15 && fabs(pv[9] - pv[15])<dx)
//    {
//        pv[15] = pv[9] + dx;
//        if (pv[15] > 0.95) pv[15] = pv[9] - dx;
//    }
//    if (index == 15 && fabs(pv[11] - pv[15])<dx)
//    {
//        pv[15] = pv[11] + dx;
//        if (pv[15] > 0.95) pv[15] = pv[11] - dx;
//    }
}

void IFunctional::print(unsigned int i, const DoubleVector &x, const DoubleVector &g, double f, GradientMethod::MethodResult result) const
{
    IFunctional* ifunc = const_cast<IFunctional*>(this);
    OptimizeParameter o_prm;
    fromVector(x, o_prm);

    IProblem2HForward2D forward;
    forward.setTimeDimension(mTimeDimension);
    forward.addSpaceDimension(mSpaceDimensionX);
    forward.addSpaceDimension(mSpaceDimensionY);
    forward.mEquParameter = mEquParameter;
    forward.mOptParameter = o_prm;

    DoubleMatrix u;
    DoubleMatrix ut;

    std::vector<SpacePointInfo> info;
    forward.calculateMVD(u, ut, info, false);

    //printf("optimizeK:%d optimizeZ:%d optimizeC:%d optimizeO:%d\n", optimizeK, optimizeZ, optimizeC, optimizeO);
    if (result == GradientMethod::FIRST_ITERATION)
    {
        //printf("Nt:%d Nx:%d Ny:%d optimizeK:%d optimizeZ:%d optimizeC:%d optimizeO:%d\n", mTimeDimension.sizeN(), mSpaceDimensionX.sizeN(), mSpaceDimensionY.sizeN(), optimizeK, optimizeZ, optimizeC, optimizeO);
        IPrinter::printSeperatorLine();
        printf("I[%3d]: %8.6f %8.6f %8.6f R:%.2f e:%.3f  \n", i, integral1(u, ut), integral2(u, ut), f, r, regEpsilon);
        printf("k:%8.4f %8.4f %8.4f %8.4f z:%8.4f %8.4f %8.4f %8.4f o:%9.4f %9.4f %9.4f %9.4f c:%9.4f %9.4f %9.4f %9.4f\n", x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11], x[12], x[13], x[14], x[15]);
        printf("k:%8.4f %8.4f %8.4f %8.4f z:%8.4f %8.4f %8.4f %8.4f o:%9.4f %9.4f %9.4f %9.4f c:%9.4f %9.4f %9.4f %9.4f\n", g[0], g[1], g[2], g[3], g[4], g[5], g[6], g[7], g[8], g[9], g[10], g[11], g[12], g[13], g[14], g[15]);
    }

    if (result == GradientMethod::NEXT_ITERATION)
    {
        IPrinter::printSeperatorLine();
        printf("I[%3d]: %8.6f %8.6f %8.6f R:%.2f e:%.3f  \n", i, integral1(u, ut), integral2(u, ut), f, r, regEpsilon);
        printf("k:%8.4f %8.4f %8.4f %8.4f z:%8.4f %8.4f %8.4f %8.4f o:%9.4f %9.4f %9.4f %9.4f c:%9.4f %9.4f %9.4f %9.4f\n", x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11], x[12], x[13], x[14], x[15]);
        printf("k:%8.4f %8.4f %8.4f %8.4f z:%8.4f %8.4f %8.4f %8.4f o:%9.4f %9.4f %9.4f %9.4f c:%9.4f %9.4f %9.4f %9.4f\n", g[0], g[1], g[2], g[3], g[4], g[5], g[6], g[7], g[8], g[9], g[10], g[11], g[12], g[13], g[14], g[15]);
    }

//    ifunc->optimizeK = (i%4==3);
//    ifunc->optimizeZ = (i%4==0);
//    ifunc->optimizeO = (i%4==1);
//    ifunc->optimizeC = (i%4==2);


//    printf("I[%3d]: %8.6f %8.6f %8.6f R:%.2f e:%.3f  \n", i, integral1(u, ut), integral2(u, ut), f, r, regEpsilon);
//    printf("k: %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n", x[ 0], x[ 1], x[ 2], x[ 3], x[ 4], x[ 5], x[ 6], x[ 7]);
//    printf("z: %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n", x[ 8], x[ 9], x[10], x[11], x[12], x[13], x[14], x[15]);
//    printf("o: %8.4f %8.4f %8.4f %8.4f\n",                         x[16], x[17], x[18], x[19]);
//    printf("c: %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n", x[20], x[21], x[22], x[23], x[24], x[25], x[26], x[27]);

    //IPrinter::print(x,x.length(),8,4);
    //printf("k:%8.4f %8.4f %8.4f %8.4f z:%8.4f %8.4f %8.4f %8.4f   o:%.4f %.4f %.4f %.4f   c:%.4f %.4f %.4f %.4f\n",
    //       x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11], x[12], x[13], x[14], x[15]);
    //printf("k:%8.4f %8.4f %8.4f %8.4f z:%8.4f %8.4f %8.4f %8.4f   o:%.4f %.4f %.4f %.4f   c:%.4f %.4f %.4f %.4f\n",
    //       g[0], g[1], g[2], g[3], g[4], g[5], g[6], g[7], g[8], g[9], g[10], g[11], g[12], g[13], g[14], g[15]);
    //IPrinter::print(g,g.length(),10,4);
    //DoubleVector px(x.length());
    //IGradient::Gradient(ifunc, 0.01, x, px);
    //IPrinter::print(px,px.length(),10,4);
    u.clear();
    ut.clear();

//    DoubleVector v1;
//    DoubleVector v2;

//    for (int i=0; i<=mTimeDimension.sizeN(); i++)
//    {
//        v1 << mParameter.k[0][0]*(info[0].value(i) - mParameter.z[0][0]) + mParameter.k[0][1]*(info[1].value(i) - mParameter.z[0][1]);
//        v2 << mParameter.k[1][0]*(info[0].value(i) - mParameter.z[1][0]) + mParameter.k[1][1]*(info[1].value(i) - mParameter.z[1][1]);
//    }

//    IPrinter::printVector(v1,"v1: "); v1.clear();
//    IPrinter::printVector(v2,"v2: "); v2.clear();

//    if (result == GradientMethod::BREAK_DISTANCE_LESS || result == GradientMethod::BREAK_GRADIENT_NORM_LESS)
//    {
//        IPrinter::printSeperatorLine();
//    }

//    ifunc->optimizeK = i%2==1;
//    ifunc->optimizeZ = i%2==1;
//    ifunc->optimizeC = i%2==0;
//    ifunc->optimizeO = i%2==0;

//    ifunc->optimizeK = i%4==0;
//    ifunc->optimizeZ = i%4==1;
//    ifunc->optimizeC = i%4==2;
//    ifunc->optimizeO = i%4==3;

//    ifunc->optimizeK = i%4==3;
//    ifunc->optimizeZ = i%4==0;
//    ifunc->optimizeC = i%4==1;
//    ifunc->optimizeO = i%4==2;
}

void IFunctional::toVector(const OptimizeParameter &prm, DoubleVector &pv) const
{
    unsigned int Nc = mEquParameter.Nc;
    unsigned int No = mEquParameter.No;

    pv.clear();
    pv.resize(2*Nc*No+2*No+2*Nc);

    for (unsigned int i=0; i<Nc; i++)
    {
        for (unsigned int j=0; j<No; j++)
        {
            pv[i*No + j] = prm.k[i][j];
        }
    }

    for (unsigned int i=0; i<Nc; i++)
    {
        for (unsigned int j=0; j<No; j++)
        {
            pv[i*No + j + Nc*No] = prm.z[i][j];
        }
    }


    for (unsigned int j=0; j<No; j++)
    {
        pv[2*j + 0 + 2*Nc*No] = prm.xi[j].x;
        pv[2*j + 1 + 2*Nc*No] = prm.xi[j].y;
    }



    for (unsigned int i=0; i<Nc; i++)
    {
        pv[2*i + 0 + 2*No + 2*Nc*No] = prm.eta[i].x;
        pv[2*i + 1 + 2*No + 2*Nc*No] = prm.eta[i].y;
    }
}

void IFunctional::fromVector(const DoubleVector &pv, OptimizeParameter &prm) const
{
    unsigned int Nc = mEquParameter.Nc;
    unsigned int No = mEquParameter.No;

    unsigned int index = 0;

    prm.k.clear();
    prm.k.resize(Nc, No);

    for (unsigned int i=0; i<Nc; i++)
    {
        for (unsigned int j=0; j<No; j++)
        {
            prm.k[i][j] = pv[index]; index++;
        }
    }

    prm.z.clear();
    prm.z.resize(Nc, No);

    for (unsigned int i=0; i<Nc; i++)
    {
        for (unsigned int j=0; j<No; j++)
        {
            prm.z[i][j] = pv[index]; index++;
        }
    }

    prm.xi.clear();
    prm.xi.resize(No);

    for (unsigned int j=0; j<No; j++)
    {
        prm.xi[j].x = pv[index]; index++;
        prm.xi[j].y = pv[index]; index++;
    }

    prm.eta.clear();
    prm.eta.resize(Nc);

    for (unsigned int i=0; i<Nc; i++)
    {
        prm.eta[i].x = pv[index]; index++;
        prm.eta[i].y = pv[index]; index++;
    }
}

