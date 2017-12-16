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

double AbstactProblem22D::fx(const DoubleVector &prms) const
{
    AbstactProblem22D* ap22d = const_cast<AbstactProblem22D*>(this);
    ap22d->msetting.fromVector(prms);
    ap22d->forward->setSetting(msetting);

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

    sum += 0.25*(u[0][0]  - U[0][0])   * (u[0][0]   - U[0][0]);
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

    msetting.fromVector(prms);
    forward->setSetting(msetting);
    backward->setSetting(msetting);

    DoubleMatrix u;
    DoubleMatrix p;

    vector<ExtendedSpaceNode2D> u_info;
    forward->calculateMVD(u, u_info);

    //    IPrinter::printSeperatorLine();
    //    IPrinter::printMatrix(U,10,10);
    //    IPrinter::printSeperatorLine();
    //    IPrinter::printMatrix(u ,10,10);
    //    IPrinter::printSeperatorLine();


    backward->u = &u;
    backward->U = &U;
    vector<ExtendedSpaceNode2D> p_info;
    backward->calculateMVD(p, p_info);

    g.clear();
    g.resize(prms.length(), 0.0);
    unsigned int gi = 0;

    // k
    for (unsigned int i=0; i<msetting.Lc; i++)
    {
        const SpaceNodePDE &eta = msetting.eta[i];
        ExtendedSpaceNode2D &pi = p_info[i];

        for (unsigned int j=0; j<msetting.Lo; j++)
        {
            const SpaceNodePDE &xi = msetting.xi[j];
            ExtendedSpaceNode2D &uj = u_info[j];

            double gradKij = 0.0;

            gradKij += 0.5 * pi.value(0) * (uj.value(0) - msetting.z[i][j]);
            for (unsigned int m=1; m<=L-1; m++)
            {
                gradKij += pi.value(m) * (uj.value(m) - msetting.z[i][j]);
            }
            gradKij += 0.5 * pi.value(L) * (uj.value(L) - msetting.z[i][j]);
            gradKij *= -ht;
            g[gi++] = gradKij;
        }
    }

    // z
    for (unsigned int i=0; i<msetting.Lc; i++)
    {
        const SpaceNodePDE &eta = msetting.eta[i];
        ExtendedSpaceNode2D &pi = p_info[i];

        for (unsigned int j=0; j<msetting.Lo; j++)
        {
            double gradZij = 0.0;

            gradZij += 0.5 * pi.value(0) * msetting.k[i][j];
            for (unsigned int m=1; m<=L-1; m++)
            {
                gradZij += pi.value(m)  * msetting.k[i][j];
            }
            gradZij += 0.5 * pi.value(L) * msetting.k[i][j];
            gradZij *= ht;
            g[gi++] = gradZij;
        }
    }

    // eta
    for (unsigned int i=0; i<msetting.Lc; i++)
    {
        const SpaceNodePDE &eta = msetting.eta[i];
        ExtendedSpaceNode2D &pi = p_info[i];

        double gradEtaiX = 0.0;
        double gradEtaiY = 0.0;
        double vi = 0.0;
        double h = 0.01;

        vi = 0.0;
        for (unsigned int j=0; j<msetting.Lo; j++) vi += msetting.k[i][j]*(u_info[j].value(0) - msetting.z[i][j]);
        gradEtaiX += 0.5 * pi.valueDxN(0,h) * vi;
        gradEtaiY += 0.5 * pi.valueDyN(0,h) * vi;

        for (unsigned int m=1; m<=L-1; m++)
        {
            vi = 0.0;
            for (unsigned int j=0; j<msetting.Lo; j++) vi += msetting.k[i][j]*(u_info[j].value(m) - msetting.z[i][j]);
            gradEtaiX += pi.valueDxN(m,h) * vi;
            gradEtaiY += pi.valueDyN(m,h) * vi;
        }

        vi = 0.0;
        for (unsigned int j=0; j<msetting.Lo; j++) vi += msetting.k[i][j]*(u_info[j].value(L) - msetting.z[i][j]);
        gradEtaiX += 0.5 * pi.valueDxN(L,h) * vi;
        gradEtaiY += 0.5 * pi.valueDyN(L,h) * vi;

        gradEtaiX *= -ht;
        gradEtaiY *= -ht;

        g[gi++] = gradEtaiX;
        g[gi++] = gradEtaiY;
    }

    // xi
    for (unsigned int j=0; j<msetting.Lo; j++)
    {
        const SpaceNodePDE &xi = msetting.xi[j];
        ExtendedSpaceNode2D &uj = u_info[j];

        double gradXijX = 0.0;
        double gradXijY = 0.0;
        double vi = 0.0;
        double h = 0.01;

        vi = 0.0;
        for (unsigned int i=0; i<msetting.Lc; i++) vi += msetting.k[i][j] * p_info[i].value(0);
        gradXijX += 0.5 * uj.valueDxN(0,h) * vi;
        gradXijY += 0.5 * uj.valueDyN(0,h) * vi;

        for (unsigned int m=1; m<=L-1; m++)
        {
            vi = 0.0;
            for (unsigned int i=0; i<msetting.Lc; i++) vi += msetting.k[i][j]*p_info[i].value(m);
            gradXijX += 0.5 * uj.valueDxN(m,h) * vi;
            gradXijY += 0.5 * uj.valueDyN(m,h) * vi;
        }

        vi = 0.0;
        for (unsigned int i=0; i<msetting.Lc; i++) vi += msetting.k[i][j]*p_info[i].value(L);
        gradXijX += 0.5 * uj.valueDxN(L,h) * vi;
        gradXijY += 0.5 * uj.valueDyN(L,h) * vi;

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

const P2Setting &AbstactProblem22D::setting() const
{
    return msetting;
}

void AbstactProblem22D::setP2Setting(const P2Setting &setting)
{
    msetting = setting;
    forward->setSetting(setting);
    backward->setSetting(setting);
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

void Problem2Forward2DEx4::layerInfo(const DoubleMatrix & u, unsigned int n) const
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

void Problem2Backward2DEx4::layerInfo(const DoubleMatrix &, unsigned int) const {}

//-------------------------------------------------------------------------------------------------------//
