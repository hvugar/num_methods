#include "abstractproblem22d.h"

AbstactProblem22D::AbstactProblem22D()
{
    forward = new Problem2Forward2DEx4();
    backward = new Problem2Backward2DEx4();

    forward->setEquationParameters(0.4, 0.01, 1.0, 10.0);
    forward->fi = 0.0;
    backward->setEquationParameters(0.4, 0.01, 1.0, 10.0);
    backward->ap22d = this;

    espilon = 0.001;
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
    g.setR1MinimizeEpsilon(0.1, 0.0001);
    g.setNormalize(true);
    g.showEndMessage(true);
    g.setResetIteration(false);

    g.calculate(prm0);
}

double AbstactProblem22D::fx(const DoubleVector &prms) const
{
    //puts("AbstactProblem22D::fx...");
    AbstactProblem22D* ap22d = const_cast<AbstactProblem22D*>(this);
    ap22d->mParameter.fromVector(prms);
    ap22d->forward->setParamter(mParameter);

    DoubleMatrix u;
    vector<ExtendedSpaceNode2D> info;
    forward->calculateMVD(u, info, false);

    double intgrl = integral(u);
    u.clear();
    //puts("AbstactProblem22D::fx.");

    return intgrl;// + espilon*norm();
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

void AbstactProblem22D::gradient(const DoubleVector &prms, DoubleVector &g)
{
    //puts("AbstactProblem22D::gradient...");
    unsigned int L = mTimeDimension.sizeN();
    double ht = mTimeDimension.step();
    //double hx = mSpaceDimensionX.step();
    //double hy = mSpaceDimensionY.step();

    mParameter.fromVector(prms);
    forward->setParamter(mParameter);
    backward->setParamter(mParameter);

    DoubleMatrix u;
    DoubleMatrix p;

    vector<ExtendedSpaceNode2D> u_info;
    forward->calculateMVD(u, u_info, true);

    backward->u = &u;
    backward->U = &U;
    vector<ExtendedSpaceNode2D> p_info;
    backward->calculateMVD(p, p_info, true);

    g.clear();
    g.resize(prms.length(), 0.0);
    unsigned int gi = 0;

    // k
    for (unsigned int i=0; i<mParameter.Lc; i++)
    {
        //const SpaceNodePDE &eta = mParameter.eta[i];
        ExtendedSpaceNode2D &pi = p_info[i];

        for (unsigned int j=0; j<mParameter.Lo; j++)
        {
            const SpaceNodePDE &xi = mParameter.xi[j];
            ExtendedSpaceNode2D &uj = u_info[j];

            double grad_Kij = 0.0;

            grad_Kij += 0.5 * pi.value(0) * (uj.value(0) - mParameter.z.at(i,j));
            for (unsigned int m=1; m<=L-1; m++)
            {
                grad_Kij += pi.value(m) * (uj.value(m) - mParameter.z.at(i,j));
            }
            grad_Kij += 0.5 * pi.value(L) * (uj.value(L) - mParameter.z.at(i,j));
            grad_Kij *= -ht;
            g[gi++] = grad_Kij;// + 2.0*espilon*(mParameter.k.at(i,j) - mParameter0.k.at(i,j));
        }
    }

    // z
    for (unsigned int i=0; i<mParameter.Lc; i++)
    {
        //const SpaceNodePDE &eta = mParameter.eta[i];
        ExtendedSpaceNode2D &pi = p_info[i];

        for (unsigned int j=0; j<mParameter.Lo; j++)
        {
            double grad_Zij = 0.0;

            grad_Zij += 0.5 * pi.value(0) * mParameter.k.at(i,j);
            for (unsigned int m=1; m<=L-1; m++)
            {
                grad_Zij += pi.value(m)  * mParameter.k.at(i,j);
            }
            grad_Zij += 0.5 * pi.value(L) * mParameter.k.at(i,j);
            grad_Zij *= ht;
            g[gi++] = grad_Zij;// + 2.0*espilon*(mParameter.z[i][j] - mParameter0.z[i][j]);
        }
    }

    // eta
    for (unsigned int i=0; i<mParameter.Lc; i++)
    {
        //const SpaceNodePDE &eta = mParameter.eta[i];
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

        g[gi++] = grad_EtaiX;// + 2.0*espilon*(mParameter.eta[i].x - mParameter0.eta[i].x);
        g[gi++] = grad_EtaiY;// + 2.0*espilon*(mParameter.eta[i].y - mParameter0.eta[i].y);
    }

    // xi
    for (unsigned int j=0; j<mParameter.Lo; j++)
    {
        //const SpaceNodePDE &xi = mParameter.xi[j];
        ExtendedSpaceNode2D &uj = u_info[j];

        double gradXijX = 0.0;
        double gradXijY = 0.0;
        double vi = 0.0;

        vi = 0.0;
        for (unsigned int i=0; i<mParameter.Lc; i++) vi += mParameter.k.at(i,j) * p_info[i].value(0);
        gradXijX += 0.5 * uj.valueDx(0) * vi;
        gradXijY += 0.5 * uj.valueDy(0) * vi;

        for (unsigned int m=1; m<=L-1; m++)
        {
            vi = 0.0;
            for (unsigned int i=0; i<mParameter.Lc; i++) vi += mParameter.k.at(i,j)*p_info[i].value(m);
            gradXijX += uj.valueDx(m) * vi;
            gradXijY += uj.valueDy(m) * vi;
        }

        vi = 0.0;
        for (unsigned int i=0; i<mParameter.Lc; i++) vi += mParameter.k.at(i,j)*p_info[i].value(L);
        gradXijX += 0.5 * uj.valueDx(L) * vi;
        gradXijY += 0.5 * uj.valueDy(L) * vi;

        gradXijX *= -ht;
        gradXijY *= -ht;

        g[gi++] = gradXijX;// + 2.0*espilon*(mParameter.xi[j].x - mParameter0.xi[j].x);;
        g[gi++] = gradXijY;// + 2.0*espilon*(mParameter.xi[j].x - mParameter0.xi[j].x);;
    }


    u_info.clear();
    p_info.clear();
    //puts("AbstactProblem22D::gradient.");
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

const Parameter &AbstactProblem22D::parameter() const
{
    return mParameter;
}

void AbstactProblem22D::setParameter(const Parameter &parameter)
{
    mParameter = parameter;
    forward->setParamter(parameter);
    backward->setParamter(parameter);
}

void AbstactProblem22D::print(unsigned int i, const DoubleVector & x, const DoubleVector & g, double, GradientMethod::MethodResult) const
{
    printf("I[%d]: %10.6f\n", i, fx(x));
    IPrinter::print(x);
    IPrinter::print(g);
    IPrinter::printSeperatorLine();
}

void AbstactProblem22D::project(DoubleVector &prm, unsigned int index)
{
    if ( 2*mParameter.Lc*mParameter.Lo < index && index < 2*mParameter.Lc*mParameter.Lo + 2*mParameter.Lc + 2*mParameter.Lo - 1 )
    {
        if (prm[index] < 0.0) prm[index] = 0.0;
        if (prm[index] > 1.0) prm[index] = 1.0;
    }
}

void AbstactProblem22D::calculateU(const Parameter &prm0)
{
    mParameter0 = prm0;
    forward->setParamter(mParameter0);
    vector<ExtendedSpaceNode2D> info;
    forward->calculateMVD(U, info, false);
}
