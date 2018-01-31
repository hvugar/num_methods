#include "ifunctional.h"
#include "problem2forward2d.h"
#include "problem2backward2d.h"
#include <math.h>

IFunctional::IFunctional()
{
    forward = new Problem2Forward2D();
    backward = new Problem2Backward2D();
    backward->func = this;
    optimizeK = true;
    optimizeZ = false;
    optimizeC = false;
    optimizeO = false;
}

double IFunctional::fx(const DoubleVector &pv) const
{
    IFunctional* ifunc = const_cast<IFunctional*>(this);

    ifunc->fromVector(pv, ifunc->mParameter);
    //functional->mParameter.fromVector(prms);
    ifunc->forward->setParameter(mParameter);
    ifunc->backward->setParameter(mParameter);

    DoubleMatrix u;
    vector<ExtendedSpaceNode2D> info;
    ifunc->forward->calculateMVD(u, info, true);

    double intgrl = integral(u);
    u.clear();

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

void IFunctional::gradient(const DoubleVector &pv, DoubleVector &g)
{
    unsigned int L = mTimeDimension.sizeN();
    double ht = mTimeDimension.step();

    fromVector(pv, mParameter);
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
    g.resize(pv.length(), 0.0);
    unsigned int gi = 0;

    // k
    if (optimizeK)
    {
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

                g[gi] = grad_Kij + 2.0*epsilon*(mParameter.k.at(i,j) - mParameter0.k.at(i,j)); gi++;
            }
        }
    }
    else
    {
        for (unsigned int i=0; i<mParameter.Lc; i++)
        {
            for (unsigned int j=0; j<mParameter.Lo; j++)
            {
                g[gi] = 0.0; gi++;
            }
        }
    }

    // z
    if (optimizeZ)
    {
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

                g[gi] = grad_Zij + 2.0*epsilon*(mParameter.z[i][j] - mParameter0.z[i][j]); gi++;
            }
        }
    }
    else
    {
        for (unsigned int i=0; i<mParameter.Lc; i++)
        {
            for (unsigned int j=0; j<mParameter.Lo; j++)
            {
                g[gi] = 0.0; gi++;
            }
        }
    }

    // eta
    if (optimizeC)
    {
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

            g[gi] = grad_EtaiX + 2.0*epsilon*(mParameter.eta[i].x - mParameter0.eta[i].x); gi++;
            g[gi] = grad_EtaiY + 2.0*epsilon*(mParameter.eta[i].y - mParameter0.eta[i].y); gi++;
        }
    }
    else
    {
        for (unsigned int i=0; i<mParameter.Lc; i++)
        {
            g[gi] = 0.0; gi++;
            g[gi] = 0.0; gi++;
        }
    }

    // xi
    if (optimizeO)
    {
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

            g[gi] = gradXijX + 2.0*epsilon*(mParameter.xi[j].x - mParameter0.xi[j].x); gi++;
            g[gi] = gradXijY + 2.0*epsilon*(mParameter.xi[j].x - mParameter0.xi[j].x); gi++;
        }
    }
    else
    {
        for (unsigned int j=0; j<mParameter.Lo; j++)
        {
            g[gi] = 0.0; gi++;
            g[gi] = 0.0; gi++;
        }
    }

    //IGradient::Gradient(this, 0.01, prms, g);

    u_info.clear();
    p_info.clear();
}

void IFunctional::project(DoubleVector &pv, unsigned int index)
{
    unsigned int offset = 2*mParameter.Lc*mParameter.Lo;

    // eta
    if ( offset <= index && index <= offset + 2*mParameter.Lc - 1)
    {
        if (pv[index] < 0.05) pv[index] = 0.05;
        if (pv[index] > 0.95) pv[index] = 0.95;
    }

    // xi
    if ( offset + 2*mParameter.Lc <= index && index <= offset + 2*mParameter.Lc + 2*mParameter.Lo - 1)
    {
        if (pv[index] < 0.05) pv[index] = 0.05;
        if (pv[index] > 0.95) pv[index] = 0.95;
    }

    double dx = 0.08;

    if (fabs(pv[8] - pv[12])<dx)
    {
        pv[12] = pv[8] + dx;
        if (pv[12] > 0.95) pv[12] = pv[8] - dx;
    }
    if (fabs(pv[10] - pv[12])<dx)
    {
        pv[12] = pv[10] + dx;
        if (pv[12] > 0.95) pv[12] = pv[10] - dx;
    }

    if (fabs(pv[8] - pv[14])<dx)
    {
        pv[14] = pv[8] + 0.08;
        if (pv[14] > 0.95) pv[14] = pv[8] - dx;
    }
    if (fabs(pv[10] - pv[14])<dx)
    {
        pv[14] = pv[10] + 0.08;
        if (pv[14] > 0.95) pv[14] = pv[10] - dx;
    }

    if (fabs(pv[9] - pv[13])<dx)
    {
        pv[13] = pv[9] + dx;
        if (pv[13] > 0.95) pv[13] = pv[9] - dx;
    }
    if (fabs(pv[11] - pv[13])<dx)
    {
        pv[13] = pv[11] + dx;
        if (pv[13] > 0.95) pv[13] = pv[11] - dx;
    }

    if (fabs(pv[9] - pv[15])<dx)
    {
        pv[15] = pv[9] + dx;
        if (pv[15] > 0.95) pv[15] = pv[9] - dx;
    }
    if (fabs(pv[11] - pv[15])<dx)
    {
        pv[15] = pv[11] + dx;
        if (pv[15] > 0.95) pv[15] = pv[11] - dx;
    }

    //if ( 2*mParameter.Lc*mParameter.Lo <= index && index <= 2*mParameter.Lc*mParameter.Lo + 2*mParameter.Lc + 2*mParameter.Lo - 1 )
    //{
    //    if (prm[index] < 0.05) prm[index] = 0.05;
    //    if (prm[index] > 0.95) prm[index] = 0.95;
    //}
}

void IFunctional::print(unsigned int i, const DoubleVector &x, const DoubleVector &g, double f, GradientMethod::MethodResult result) const
{
    IFunctional* ifunc = const_cast<IFunctional*>(this);
    Parameter prm; prm.Lc = mParameter.Lc; prm.Lo = mParameter.Lo;
    ifunc->fromVector(x, prm);

    ifunc->fromVector(x, ifunc->mParameter);
    ifunc->forward->setParameter(ifunc->mParameter);

    DoubleMatrix u;
    vector<ExtendedSpaceNode2D> info;
    ifunc->forward->calculateMVD(u, info, true);

    printf("I[%3d]: %10.6f %10.6f %14.6f R:%.3f epsilon:%.3f Nt:%d Nx:%d Ny:%d  ", i, integral(u), norm(),  f, r, epsilon, mTimeDimension.sizeN(), mSpaceDimensionX.sizeN(), mSpaceDimensionY.sizeN());
    IPrinter::print(x,x.length(),8,4);
    //IPrinter::print(g,g.length(),10,4);
    //DoubleVector px(x.length());
    //IGradient::Gradient(ifunc, 0.01, x, px);
    //IPrinter::print(px,px.length(),10,4);

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

    ifunc->optimizeK = i%4==3;
    ifunc->optimizeZ = i%4==0;
    ifunc->optimizeC = i%4==1;
    ifunc->optimizeO = i%4==2;


//    if (i%2==0)
//    {
//        ifunc->optimizeK = true;
//        ifunc->optimizeZ = true;
//        ifunc->optimizeC = false;
//        ifunc->optimizeO = false;
//    }
//    else
//    {
//        ifunc->optimizeK = false;
//        ifunc->optimizeZ = false;
//        ifunc->optimizeC = true;
//        ifunc->optimizeO = true;
//    }

//    IPrinter::printSeperatorLine();
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

void IFunctional::toVector(const Parameter &prm, DoubleVector &pv) const
{
    unsigned int Lc = prm.Lc;
    unsigned int Lo = prm.Lo;

    //unsigned int length = 0;
    //if (optimizeK) length += Lc*Lo;
    //if (optimizeZ) length += Lc*Lo;
    //if (optimizeC) length += 2*Lc;
    //if (optimizeO) length += 2*Lo;

    pv.clear();
    pv.resize(2*Lc*Lo+2*Lc+2*Lo);

    //if (optimizeK)
    {
        for (unsigned int i=0; i<Lc; i++)
        {
            for (unsigned int j=0; j<Lo; j++)
            {
                pv[i*Lo + j] = prm.k[i][j];
            }
        }
    }

    //if (optimizeZ)
    {
        for (unsigned int i=0; i<Lc; i++)
        {
            for (unsigned int j=0; j<Lo; j++)
            {
                pv[i*Lo + j + Lc*Lo] = prm.z[i][j];
            }
        }
    }

    //if (optimizeC)
    {
        for (unsigned int i=0; i<Lc; i++)
        {
            pv[2*i + 0 + 2*Lc*Lo] = prm.eta[i].x;
            pv[2*i + 1 + 2*Lc*Lo] = prm.eta[i].y;
        }
    }

    //if (optimizeO)
    {
        for (unsigned int j=0; j<Lo; j++)
        {
            pv[2*j + 0 + 2*Lc + 2*Lc*Lo] = prm.xi[j].x;
            pv[2*j + 1 + 2*Lc + 2*Lc*Lo] = prm.xi[j].y;
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

void IFunctional::fromVector(const DoubleVector &pv, Parameter &prm)
{
    unsigned int Lc = prm.Lc;
    unsigned int Lo = prm.Lo;

    unsigned int index = 0;

    //if (optimizeK)
    {
        prm.k.clear();
        prm.k.resize(Lc, Lo);

        for (unsigned int i=0; i<Lc; i++)
        {
            for (unsigned int j=0; j<Lo; j++)
            {
                prm.k[i][j] = pv[index]; index++;
            }
        }
    }

    //if (optimizeZ)
    {
        prm.z.clear();
        prm.z.resize(Lc, Lo);

        for (unsigned int i=0; i<Lc; i++)
        {
            for (unsigned int j=0; j<Lo; j++)
            {
                prm.z[i][j] = pv[index]; index++;
            }
        }
    }

    //if (optimizeC)
    {
        prm.eta.clear();
        prm.eta.resize(Lc);

        for (unsigned int i=0; i<Lc; i++)
        {
            prm.eta[i].x = pv[index]; index++;
            prm.eta[i].y = pv[index]; index++;
        }
    }

    //if (optimizeO)
    {
        prm.xi.clear();
        prm.xi.resize(Lo);

        for (unsigned int j=0; j<Lo; j++)
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
