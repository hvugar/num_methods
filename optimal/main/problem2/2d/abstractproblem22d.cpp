#include "abstractproblem22d.h"

AbstactProblem22D::AbstactProblem22D()
{
    forward = new Problem2Forward2DEx4();
    backward = new Problem2Backward2DEx4();

    forward->setEquationParameters(0.4, 0.01, 1.0, 10.0);
    forward->fi = 0.0;
    backward->setEquationParameters(0.4, 0.01, 1.0, 10.0);
    backward->ap22d = this;
}

AbstactProblem22D::~AbstactProblem22D()
{
    delete forward;
    delete backward;
}

double AbstactProblem22D::fx(const DoubleVector &prms) const
{
    AbstactProblem22D* ap22d = const_cast<AbstactProblem22D*>(this);
    ap22d->mParameter.fromVector(prms);
    ap22d->forward->setParamter(mParameter);

    DoubleMatrix u;
    vector<ExtendedSpaceNode2D> info;
    forward->calculateMVD(u, info);

    double intgrl = integral(u);
    u.clear();

    return intgrl;
}

double AbstactProblem22D::integral(const DoubleMatrix &u) const
{
    double hx = mSpaceDimensionX.step();
    double hy = mSpaceDimensionY.step();
    unsigned int N1 = mSpaceDimensionX.sizeN();
    unsigned int N2 = mSpaceDimensionY.sizeN();

    double sum = 0.0;

    sum += 0.25*(u[0][0]  - U[0][0])    * (u[0][0]   - U[0][0]);
    sum += 0.25*(u[0][N1]  - U[0][N1])  * (u[0][N1]  - U[0][N1]);
    sum += 0.25*(u[N2][0]  - U[N2][0])  * (u[N2][0]  - U[N2][0]);
    sum += 0.25*(u[N2][N1] - U[N2][N1]) * (u[N2][N1] - U[N2][N1]);

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

double AbstactProblem22D::mu(double x UNUSED_PARAM, double y UNUSED_PARAM) const
{
    return 1.0;
}

void AbstactProblem22D::gradient(const DoubleVector &prms, DoubleVector &g)
{
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
    forward->calculateMVD(u, u_info);

    backward->u = &u;
    backward->U = &U;
    vector<ExtendedSpaceNode2D> p_info;
    backward->calculateMVD(p, p_info);

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

            grad_Kij += 0.5 * pi.value(0) * (uj.value(0) - mParameter.z[i][j]);
            for (unsigned int m=1; m<=L-1; m++)
            {
                grad_Kij += pi.value(m) * (uj.value(m) - mParameter.z[i][j]);
            }
            grad_Kij += 0.5 * pi.value(L) * (uj.value(L) - mParameter.z[i][j]);
            grad_Kij *= -ht;
            g[gi++] = grad_Kij;
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

            grad_Zij += 0.5 * pi.value(0) * mParameter.k[i][j];
            for (unsigned int m=1; m<=L-1; m++)
            {
                grad_Zij += pi.value(m)  * mParameter.k[i][j];
            }
            grad_Zij += 0.5 * pi.value(L) * mParameter.k[i][j];
            grad_Zij *= ht;
            g[gi++] = grad_Zij;
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
        for (unsigned int j=0; j<mParameter.Lo; j++) vi += mParameter.k[i][j]*(u_info[j].value(0) - mParameter.z[i][j]);
        grad_EtaiX += 0.5 * /*pi.valueDxN(0,hx)*/pi.valueDx(0) * vi;
        grad_EtaiY += 0.5 * /*pi.valueDyN(0,hy)*/pi.valueDy(0) * vi;

        for (unsigned int m=1; m<=L-1; m++)
        {
            vi = 0.0;
            for (unsigned int j=0; j<mParameter.Lo; j++) vi += mParameter.k[i][j]*(u_info[j].value(m) - mParameter.z[i][j]);
            grad_EtaiX += /*pi.valueDxN(m,hx)*/pi.valueDx(m) * vi;
            grad_EtaiY += /*pi.valueDyN(m,hy)*/pi.valueDy(m) * vi;
        }

        vi = 0.0;
        for (unsigned int j=0; j<mParameter.Lo; j++) vi += mParameter.k[i][j]*(u_info[j].value(L) - mParameter.z[i][j]);
        grad_EtaiX += 0.5 * /*pi.valueDxN(L,hx)*/pi.valueDx(L) * vi;
        grad_EtaiY += 0.5 * /*pi.valueDyN(L,hy)*/pi.valueDy(L) * vi;

        grad_EtaiX *= -ht;
        grad_EtaiY *= -ht;

        g[gi++] = grad_EtaiX;
        g[gi++] = grad_EtaiY;
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
        for (unsigned int i=0; i<mParameter.Lc; i++) vi += mParameter.k[i][j] * p_info[i].value(0);
        gradXijX += 0.5 * /*uj.valueDxN(0,hx)*/uj.valueDx(0) * vi;
        gradXijY += 0.5 * /*uj.valueDyN(0,hy)*/uj.valueDy(0) * vi;

        for (unsigned int m=1; m<=L-1; m++)
        {
            vi = 0.0;
            for (unsigned int i=0; i<mParameter.Lc; i++) vi += mParameter.k[i][j]*p_info[i].value(m);
            gradXijX += /*uj.valueDxN(m,hx)*/uj.valueDx(m) * vi;
            gradXijY += /*uj.valueDyN(m,hy)*/uj.valueDy(m) * vi;
        }

        vi = 0.0;
        for (unsigned int i=0; i<mParameter.Lc; i++) vi += mParameter.k[i][j]*p_info[i].value(L);
        gradXijX += 0.5 * /*uj.valueDxN(L,hx)*/uj.valueDx(L) * vi;
        gradXijY += 0.5 * /*uj.valueDyN(L,hy)*/uj.valueDy(L) * vi;

        gradXijX *= -ht;
        gradXijY *= -ht;

        g[gi++] = gradXijX;
        g[gi++] = gradXijY;
    }

    u_info.clear();
    p_info.clear();
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

//-------------------------------------------------------------------------------------------------------//

double Problem2Forward2DEx4::initial(const SpaceNodePDE &) const
{
    return fi;
}

double Problem2Forward2DEx4::boundary(const SpaceNodePDE &, const TimeNodePDE &, BoundaryType) const { return NAN; }
double Problem2Forward2DEx4::f(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }
double Problem2Forward2DEx4::g1(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }
double Problem2Forward2DEx4::g2(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }
double Problem2Forward2DEx4::g3(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }
double Problem2Forward2DEx4::g4(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }

void Problem2Forward2DEx4::layerInfo(const DoubleMatrix &u, unsigned int n) const
{
    //    QPixmap pixmap;
    //    visualizeMatrixHeat(u, u.min(), u.max(), pixmap, u.cols(), u.rows());
    //    pixmap.save(QString("image%1.png").arg(n), "PNG");
}

//-------------------------------------------------------------------------------------------------------//

double Problem2Backward2DEx4::initial(const SpaceNodePDE & sn) const
{
    return -2.0 * ap22d->mu(sn.x, sn.y) * ((*u)[sn.j][sn.i] - (*U)[sn.j][sn.i]);
}

double Problem2Backward2DEx4::boundary(const SpaceNodePDE &, const TimeNodePDE &, BoundaryType) const { return NAN; }
double Problem2Backward2DEx4::f(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }
double Problem2Backward2DEx4::h(const SpaceNodePDE &) const { return 0.0; }
double Problem2Backward2DEx4::g1(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }
double Problem2Backward2DEx4::g2(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }
double Problem2Backward2DEx4::g3(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }
double Problem2Backward2DEx4::g4(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }

void Problem2Backward2DEx4::layerInfo(const DoubleMatrix &p, unsigned int ln) const
{
    if (ln==timeDimension().sizeN())
    {
        IPrinter::printSeperatorLine();
        IPrinter::printMatrix(p);
        IPrinter::printSeperatorLine();
    }
}

//-------------------------------------------------------------------------------------------------------//
