#include "abstractproblem22d.h"

AbstactProblem22D::AbstactProblem22D()
{
    forward = new Problem2Forward2D();
    backward = new Problem2Backward2D();
    //backward->func = this;
}

AbstactProblem22D::~AbstactProblem22D()
{
    delete forward;
    delete backward;
}

void AbstactProblem22D::optimization(DoubleVector &prm0)
{
    AbstactProblem22D* p = const_cast<AbstactProblem22D*>(this);

    ConjugateGradient g;
    g.setFunction(this);
    g.setGradient(this);
    g.setPrinter(this);
    g.setProjection(p);
    g.setEpsilon1(0.000);
    g.setEpsilon2(0.000);
    g.setEpsilon3(0.000);
    g.setR1MinimizeEpsilon(1.0, 0.001);
    g.setNormalize(true);
    g.showEndMessage(true);
    g.setResetIteration(false);

    g.calculate(prm0);
}

double AbstactProblem22D::fx(const DoubleVector &prms) const
{
    AbstactProblem22D* ap22d = const_cast<AbstactProblem22D*>(this);

    ap22d->mParameter.fromVector(prms);
    ap22d->forward->setParameter(mParameter);

    DoubleMatrix u;
    vector<ExtendedSpaceNode2D> info;
    ap22d->forward->calculateMVD(u, info, true);

    ap22d->backward->setPenaltyCoefficient(r);
    ap22d->backward->info = &info;

    double intgrl = integral(u);
    u.clear();

    return intgrl + epsilon*norm() + r*penalty(info);
}

double AbstactProblem22D::integral(const DoubleMatrix &u) const
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

double AbstactProblem22D::norm() const
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

double AbstactProblem22D::mu(double x UNUSED_PARAM, double y UNUSED_PARAM) const
{
    return 1.0;
}

double AbstactProblem22D::penalty(const vector<ExtendedSpaceNode2D> &info) const
{
    double ht = mTimeDimension.step();
    unsigned int L = mTimeDimension.sizeN();

    double p_sum = 0.0;

    for (unsigned int i=0; i<mParameter.Lc; i++) p_sum += 0.5*gpi(i, 0, info)*gpi(i, 0, info);
    for (unsigned int l=1; l<=L-1; l++)
    {
        for (unsigned int i=0; i<mParameter.Lc; i++) p_sum += gpi(i, l, info)*gpi(i, l, info);
    }
    for (unsigned int i=0; i<mParameter.Lc; i++) p_sum += 0.5*gpi(i, L, info)*gpi(i, L, info);
    return p_sum*ht;
}

double AbstactProblem22D::gpi(unsigned int i, unsigned int layer, const vector<ExtendedSpaceNode2D> &info) const
{
    double p = fabs(g0i(i, layer, info)) - (vmax.at(i) - vmin.at(i))/2.0;
    return p > 0.0 ? p : 0.0;
}

double AbstactProblem22D::g0i(unsigned int i, unsigned int layer, const vector<ExtendedSpaceNode2D> &info) const
{
    double vi = 0.0;
    for (unsigned int j=0; j<mParameter.Lo; j++)
    {
        const ExtendedSpaceNode2D node = info.at(j);
        vi += mParameter.k[i][j] * (node.value(layer)-mParameter.z[i][j]);
    }
    return (vmax.at(i) + vmin.at(i))/2.0 - vi;
}

void AbstactProblem22D::gradient(const DoubleVector &prms, DoubleVector &g)
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

void AbstactProblem22D::print(unsigned int i, const DoubleVector & x, const DoubleVector & g, double, GradientMethod::MethodResult) const
{
    printf("I[%d]: %10.6f\n", i, fx(x));
    IPrinter::print(x,x.length(),8,4);
    IPrinter::print(g,g.length(),8,4);
    IPrinter::printSeperatorLine();
}

void AbstactProblem22D::project(DoubleVector &prm, unsigned int index)
{
    if ( 2*mParameter.Lc*mParameter.Lo <= index && index <= 2*mParameter.Lc*mParameter.Lo + 2*mParameter.Lc + 2*mParameter.Lo - 1 )
    {
        if (prm[index] < 0.05) prm[index] = 0.05;
        if (prm[index] > 0.95) prm[index] = 0.95;
    }
}


void AbstactProblem22D::setGridParameters(Dimension timeDimension, Dimension spaceDimensionX, Dimension spaceDimensionY)
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

void AbstactProblem22D::setEquationParameters(double a, double lambda0, double lambda)
{
    forward->setEquationParameters(a, lambda0, lambda);
    backward->setEquationParameters(a, lambda0, lambda);
}

void AbstactProblem22D::setIntTemperature(double fi)
{
    forward->setIntTemperature(fi);
    backward->setIntTemperature(fi);
}

void AbstactProblem22D::setEnvTemperature(double theta)
{
    forward->setEnvTemperature(theta);
    backward->setEnvTemperature(theta);
}

void AbstactProblem22D::setEpsilon(double epsilon)
{
    this->epsilon = epsilon;
}

void AbstactProblem22D::setPenaltyCoefficient(double r)
{
    this->r = r;
    backward->setPenaltyCoefficient(r);
}

void AbstactProblem22D::setParameter(const Parameter &parameter)
{
    this->mParameter = parameter;
    forward->setParameter(parameter);
    backward->setParameter(parameter);
}

void AbstactProblem22D::setParameter0(const Parameter &parameter0)
{
    this->mParameter0 = parameter0;
}

void AbstactProblem22D::setPenaltyLimits(const DoubleVector &vmin, const DoubleVector &vmax)
{
    this->vmin = vmin;
    this->vmax = vmax;
}

void AbstactProblem22D::calculateU()
{
    forward->setParameter(mParameter0);
    vector<ExtendedSpaceNode2D> info;
    forward->calculateMVD(U, info, false);
}

double AbstactProblem22D::sgn(double x) const
{
    if (x<0.0) return -1.0;
    else if (x>0.0) return +1.0;
    return 0.0;
}
