#include "ifunctional.h"
#include "problem2forward2d.h"
#include "problem2backward2d.h"
#include <math.h>

IFunctional::IFunctional()
{
    forward = new Problem2Forward2D();
    backward = new Problem2Backward2D();
    backward->func = this;
}

double IFunctional::fx(const DoubleVector &prms) const
{
    IFunctional* functional = const_cast<IFunctional*>(this);

    functional->mParameter.fromVector(prms);
    functional->forward->setParameter(mParameter);
    functional->backward->setParameter(mParameter);

    DoubleMatrix u;
    vector<ExtendedSpaceNode2D> info;
    functional->forward->calculateMVD(u, info, true);

    double intgrl = integral(u);
    u.clear();

    double v1 = mParameter.k[0][0]*(info[0].value(mTimeDimension.sizeN())-mParameter.z[0][0])
            + mParameter.k[0][1]*(info[1].value(mTimeDimension.sizeN())-mParameter.z[0][1]);
    double v2 = mParameter.k[1][0]*(info[0].value(mTimeDimension.sizeN())-mParameter.z[1][0])
            + mParameter.k[1][1]*(info[1].value(mTimeDimension.sizeN())-mParameter.z[1][1]);
    //printf("%f %f\t", v1, v2);

    return intgrl + epsilon*norm() + r*penalty(info);
}

double IFunctional::integral(const DoubleMatrix &u) const
{
    double hx = mSpaceDimensionX.step();
    double hy = mSpaceDimensionY.step();
    unsigned int N1 = mSpaceDimensionX.sizeN();
    unsigned int N2 = mSpaceDimensionY.sizeN();

    double sum = 0.0;

    sum += 0.25*(u[0][0]   - U[0][0])   * (u[0][0]   - U[0][0]);
    sum += 0.25*(u[0][N1]  - U[0][N1])  * (u[0][N1]  - U[0][N1]);
    sum += 0.25*(u[N2][0]  - U[N2][0])  * (u[N2][0]  - U[N2][0]);
    sum += 0.25*(u[N2][N1] - U[N2][N1]) * (u[N2][N1] - U[N2][N1]);

    for (unsigned int n1=1; n1<=N1-1; n1++)
    {
        sum += 0.5*(u[0][n1]  - U[0][n1])*(u[0][n1] - U[0][n1]);
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

double IFunctional::norm() const
{
    double norm = 0.0;

    for (unsigned int i=0; i<mParameter.Lc; i++)
    {
        norm += (mParameter.eta[i].x - mParameter0.eta[i].x)*(mParameter.eta[i].x - mParameter0.eta[i].x);
        norm += (mParameter.eta[i].y - mParameter0.eta[i].y)*(mParameter.eta[i].y - mParameter0.eta[i].y);

        for (unsigned int j=0; j<mParameter.Lo; j++)
        {
            norm += (mParameter.k[i][j] - mParameter0.k[i][j])*(mParameter.k[i][j] - mParameter0.k[i][j]);
            norm += (mParameter.z[i][j] - mParameter0.z[i][j])*(mParameter.z[i][j] - mParameter0.z[i][j]);

            norm += (mParameter.xi[j].x - mParameter0.xi[j].x)*(mParameter.xi[j].x - mParameter0.xi[j].x);
            norm += (mParameter.xi[j].y - mParameter0.xi[j].y)*(mParameter.xi[j].y - mParameter0.xi[j].y);
        }
    }

    return norm;
}

double IFunctional::penalty(vector<ExtendedSpaceNode2D> &info) const
{
    IFunctional* functional = const_cast<IFunctional*>(this);

    double ht = mTimeDimension.step();
    unsigned int L = mTimeDimension.sizeN();

    functional->backward->setPenaltyCoefficient(r);
    functional->backward->info = &info;

    double p_sum = 0.0;

    for (unsigned int i=0; i<mParameter.Lc; i++) p_sum += 0.5*gpi(i, 0, info)*gpi(i, 0, info);
    for (unsigned int l=1; l<=L-1; l++)
    {
        for (unsigned int i=0; i<mParameter.Lc; i++) p_sum += gpi(i, l, info)*gpi(i, l, info);
    }
    for (unsigned int i=0; i<mParameter.Lc; i++) p_sum += 0.5*gpi(i, L, info)*gpi(i, L, info);
    return p_sum*ht;
}

double IFunctional::gpi(unsigned int i, unsigned int layer, const vector<ExtendedSpaceNode2D> &info) const
{
    double p = fabs(g0i(i, layer, info)) - (vmax.at(i) - vmin.at(i))/2.0;
    return p > 0.0 ? p : 0.0;
}

double IFunctional::g0i(unsigned int i, unsigned int layer, const vector<ExtendedSpaceNode2D> &info) const
{
    double vi = 0.0;
    for (unsigned int j=0; j<mParameter.Lo; j++)
    {
        const ExtendedSpaceNode2D node = info.at(j);
        vi += mParameter.k[i][j] * (node.value(layer)-mParameter.z[i][j]);
    }
    return (vmax.at(i) + vmin.at(i))/2.0 - vi;
}

double IFunctional::mu(double x UNUSED_PARAM, double y UNUSED_PARAM) const
{
    return 1.0;
}

void IFunctional::gradient(const DoubleVector &prms, DoubleVector &g)
{
    unsigned int L = mTimeDimension.sizeN();
    double ht = mTimeDimension.step();

    mParameter.fromVector(prms);
    forward->setParameter(mParameter);
    backward->setParameter(mParameter);

    DoubleMatrix u;
    DoubleMatrix p;

    vector<ExtendedSpaceNode2D> u_info;
    forward->calculateMVD(u, u_info, true);

    backward->u = &u;
    backward->U = &U;
    backward->info = &u_info;
    vector<ExtendedSpaceNode2D> p_info;
    backward->calculateMVD(p, p_info, true);

    g.clear();
    g.resize(prms.length(), 0.0);
    unsigned int gi = 0;

    // k
    for (unsigned int i=0; i<mParameter.Lc; i++)
    {
        ExtendedSpaceNode2D &pi = p_info[i];

        for (unsigned int j=0; j<mParameter.Lo; j++)
        {
            ExtendedSpaceNode2D &uj = u_info[j];

            double grad_Kij = 0.0;

            grad_Kij += 0.5 * (pi.value(0)+2.0*r*gpi(i,0,u_info)*sgn(g0i(i,0,u_info))) * (uj.value(0) - mParameter.z.at(i,j));
            for (unsigned int m=1; m<=L-1; m++)
            {
                grad_Kij += (pi.value(m)+2.0*r*gpi(i,m,u_info)*sgn(g0i(i,m,u_info))) * (uj.value(m) - mParameter.z.at(i,j));
            }
            grad_Kij += 0.5 * (pi.value(L)+2.0*r*gpi(i,L,u_info)*sgn(g0i(i,L,u_info))) * (uj.value(L) - mParameter.z.at(i,j));
            grad_Kij *= -ht;
            g[gi++] = grad_Kij + 2.0*epsilon*(mParameter.k.at(i,j) - mParameter0.k.at(i,j));
        }
    }

    // z
    for (unsigned int i=0; i<mParameter.Lc; i++)
    {
        ExtendedSpaceNode2D &pi = p_info[i];

        for (unsigned int j=0; j<mParameter.Lo; j++)
        {
            double grad_Zij = 0.0;

            grad_Zij += 0.5 * (pi.value(0)+2.0*r*gpi(i,0,u_info)*sgn(g0i(i,0,u_info))) * mParameter.k.at(i,j);
            for (unsigned int m=1; m<=L-1; m++)
            {
                grad_Zij += (pi.value(m)+2.0*r*gpi(i,m,u_info)*sgn(g0i(i,m,u_info)))  * mParameter.k.at(i,j);
            }
            grad_Zij += 0.5 * (pi.value(L)+2.0*r*gpi(i,L,u_info)*sgn(g0i(i,L,u_info))) * mParameter.k.at(i,j);
            grad_Zij *= ht;
            g[gi++] = grad_Zij + 2.0*epsilon*(mParameter.z[i][j] - mParameter0.z[i][j]);
        }
    }

    // eta
    for (unsigned int i=0; i<mParameter.Lc; i++)
    {
        ExtendedSpaceNode2D &pi = p_info[i];

        double grad_EtaiX = 0.0;
        double grad_EtaiY = 0.0;
        double vi = 0.0;

        vi = 0.0;
        for (unsigned int j=0; j<mParameter.Lo; j++) vi += mParameter.k.at(i,j) * (u_info[j].value(0) - mParameter.z.at(i,j));
        grad_EtaiX += 0.5 * pi.valueDx(0) * vi;
        grad_EtaiY += 0.5 * pi.valueDy(0) * vi;

        for (unsigned int m=1; m<=L-1; m++)
        {
            vi = 0.0;
            for (unsigned int j=0; j<mParameter.Lo; j++) vi += mParameter.k.at(i,j) * (u_info[j].value(m) - mParameter.z.at(i,j));
            grad_EtaiX += pi.valueDx(m) * vi;
            grad_EtaiY += pi.valueDy(m) * vi;
        }

        vi = 0.0;
        for (unsigned int j=0; j<mParameter.Lo; j++) vi += mParameter.k.at(i,j) * (u_info[j].value(L) - mParameter.z.at(i,j));
        grad_EtaiX += 0.5 * pi.valueDx(L) * vi;
        grad_EtaiY += 0.5 * pi.valueDy(L) * vi;

        grad_EtaiX *= -ht;
        grad_EtaiY *= -ht;

        g[gi++] = grad_EtaiX + 2.0*epsilon*(mParameter.eta[i].x - mParameter0.eta[i].x);
        g[gi++] = grad_EtaiY + 2.0*epsilon*(mParameter.eta[i].y - mParameter0.eta[i].y);
    }

    // xi
    for (unsigned int j=0; j<mParameter.Lo; j++)
    {
        ExtendedSpaceNode2D &uj = u_info[j];

        double gradXijX = 0.0;
        double gradXijY = 0.0;
        double vi = 0.0;

        vi = 0.0;
        for (unsigned int i=0; i<mParameter.Lc; i++) vi += mParameter.k.at(i,j) * (p_info[i].value(0)+2.0*r*gpi(i,0,u_info)*sgn(g0i(i,0,u_info)));
        gradXijX += 0.5 * uj.valueDx(0) * vi;
        gradXijY += 0.5 * uj.valueDy(0) * vi;

        for (unsigned int m=1; m<=L-1; m++)
        {
            vi = 0.0;
            for (unsigned int i=0; i<mParameter.Lc; i++) vi += mParameter.k.at(i,j)*(p_info[i].value(m)+2.0*r*gpi(i,m,u_info)*sgn(g0i(i,m,u_info)));
            gradXijX += uj.valueDx(m) * vi;
            gradXijY += uj.valueDy(m) * vi;
        }

        vi = 0.0;
        for (unsigned int i=0; i<mParameter.Lc; i++) vi += mParameter.k.at(i,j)*(p_info[i].value(L)+2.0*r*gpi(i,L,u_info)*sgn(g0i(i,L,u_info)));
        gradXijX += 0.5 * uj.valueDx(L) * vi;
        gradXijY += 0.5 * uj.valueDy(L) * vi;

        gradXijX *= -ht;
        gradXijY *= -ht;

        g[gi++] = gradXijX + 2.0*epsilon*(mParameter.xi[j].x - mParameter0.xi[j].x);
        g[gi++] = gradXijY + 2.0*epsilon*(mParameter.xi[j].x - mParameter0.xi[j].x);
    }


    u_info.clear();
    p_info.clear();
}

void IFunctional::project(DoubleVector &prm, unsigned int index)
{
    if ( 2*mParameter.Lc*mParameter.Lo <= index && index <= 2*mParameter.Lc*mParameter.Lo + 2*mParameter.Lc + 2*mParameter.Lo - 1 )
    {
        if (prm[index] < 0.05) prm[index] = 0.05;
        if (prm[index] > 0.95) prm[index] = 0.95;
    }
}

void IFunctional::print(unsigned int i, const DoubleVector &x, const DoubleVector &g, double f, GradientMethod::MethodResult result) const
{
    printf("I[%d]: %10.6f %10.6f R: %10.6f\n", i, fx(x), f, r);
    IPrinter::print(x,x.length(),10,4);
    IPrinter::print(g,g.length(),10,4);
}

void IFunctional::setGridParameters(Dimension timeDimension, Dimension spaceDimensionX, Dimension spaceDimensionY)
{
    mTimeDimension = timeDimension;
    mSpaceDimensionX = spaceDimensionX;
    mSpaceDimensionY = spaceDimensionY;

    forward->setTimeDimension(mTimeDimension);
    forward->addSpaceDimension(mSpaceDimensionX);
    forward->addSpaceDimension(mSpaceDimensionY);

    backward->setTimeDimension(mTimeDimension);
    backward->addSpaceDimension(mSpaceDimensionX);
    backward->addSpaceDimension(mSpaceDimensionY);
}

void IFunctional::setGridTimeDimension(Dimension timeDimension)
{
    mTimeDimension = timeDimension;
    forward->setTimeDimension(timeDimension);
    backward->setTimeDimension(timeDimension);
}

void IFunctional::setEquationParameters(double a, double lambda0, double lambda)
{
    forward->setEquationParameters(a, lambda0, lambda);
    backward->setEquationParameters(a, lambda0, lambda);
}

void IFunctional::setIntTemperature(double fi)
{
    forward->setIntTemperature(fi);
    backward->setIntTemperature(fi);
}

void IFunctional::setEnvTemperature(double theta)
{
    forward->setEnvTemperature(theta);
    backward->setEnvTemperature(theta);
}

void IFunctional::setEpsilon(double epsilon)
{
    this->epsilon = epsilon;
}

void IFunctional::setPenaltyCoefficient(double r)
{
    this->r = r;
    backward->setPenaltyCoefficient(r);
}

void IFunctional::setPenaltyLimits(const DoubleVector &vmin, const DoubleVector &vmax)
{
    this->vmin = vmin;
    this->vmax = vmax;
}

void IFunctional::setParameter(const Parameter &parameter)
{
    this->mParameter = parameter;
    forward->setParameter(parameter);
    backward->setParameter(parameter);
}

void IFunctional::setParameter0(const Parameter &parameter0)
{
    this->mParameter0 = parameter0;
}

double IFunctional::sgn(double x) const
{
    if (x<0.0) return -1.0;
    else if (x>0.0) return +1.0;
    return 0.0;
}
