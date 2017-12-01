#include "abstractproblem22d.h"

void AbstactProblem22D::setGridParameters(Dimension timeDimension, Dimension spaceDimensionX, Dimension spaceDimensionY)
{
    mTimeDimension = timeDimension;
    mSpaceDimensionX = spaceDimensionX;
    mSpaceDimensionY = spaceDimensionY;
}

double AbstactProblem22D::fx(const DoubleVector &prms) const
{
    const_cast<AbstactProblem22D*>(this)->msetting.fromVector(prms);

    forward->setTimeDimension(mTimeDimension);
    forward->addSpaceDimension(mSpaceDimensionX);
    forward->addSpaceDimension(mSpaceDimensionY);
    forward->setSettings(msetting);
    std::vector<DoubleMatrix> u;
    forward->calculateMVD(u);

    DoubleMatrix uT = u[u.size()-1];

    double intgrl = integral(uT);

    for (unsigned int i=0; i<u.size(); i++)
    {
        u[i].clear();
    }
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
    unsigned int N = mSpaceDimensionX.sizeN();
    unsigned int M = mSpaceDimensionY.sizeN();
    unsigned int L = mTimeDimension.sizeN();

    double hx = mSpaceDimensionX.step();
    double hy = mSpaceDimensionY.step();
    double ht = mTimeDimension.step();

    msetting.fromVector(prms);

    forward->setTimeDimension(mTimeDimension);
    forward->addSpaceDimension(mSpaceDimensionX);
    forward->addSpaceDimension(mSpaceDimensionY);
    forward->setSettings(msetting);
    forward->calculateMVD(u);

    DoubleMatrix uT = u[u.size()-1];

    backward->setTimeDimension(mTimeDimension);
    backward->addSpaceDimension(mSpaceDimensionX);
    backward->addSpaceDimension(mSpaceDimensionY);
    backward->setSettings(msetting);

    //backward->U = this->U;
    //backward->uT = u[u.size()-1];

    //backward->mu.clear();
    //backward->mu.resize(M+1, N+1, 1.0);

    std::vector<DoubleMatrix> p;
    backward->calculateMVD(p);

    g.resize(prms.length(), 0.0);
    unsigned int gi = 0;

    // k
    for (unsigned int i=0; i<msetting.Lc; i++)
    {
        const SpaceNodePDE &eta = msetting.eta[i];
        unsigned int etaX = (unsigned int)(round(eta.x*N));
        unsigned int etaY = (unsigned int)(round(eta.y*M));

        for (unsigned int j=0; j<msetting.Lo; j++)
        {
            const SpaceNodePDE &xi = msetting.xi[j];
            unsigned int xiX = (unsigned int)(round(xi.x*N));
            unsigned int xiY = (unsigned int)(round(xi.y*M));

            double gradKij = 0.0;
            gradKij += 0.5 * p[0][etaY][etaX] * (u[0][xiY][xiX] - msetting.z[i][j]);
            for (unsigned int m=1; m<=L-1; m++)
            {
                gradKij += p[m][etaY][etaX] * (u[m][xiY][xiX] - msetting.z[i][j]);
            }
            gradKij += 0.5 * p[L][etaY][etaX] * (u[L][xiY][xiX] - msetting.z[i][j]);
            gradKij *= -ht;
            g[gi++] = gradKij;
        }
    }

    // z
    for (unsigned int i=0; i<msetting.Lc; i++)
    {
        const SpaceNodePDE &eta = msetting.eta[i];
        unsigned int etaX = (unsigned int)(round(eta.x*N));
        unsigned int etaY = (unsigned int)(round(eta.y*M));

        for (unsigned int j=0; j<msetting.Lo; j++)
        {
            double gradZij = 0.0;
            gradZij += 0.5 * p[0][etaY][etaX] * msetting.k[i][j];
            for (unsigned int m=1; m<=L-1; m++)
            {
                gradZij += p[m][etaY][etaX]  * msetting.k[i][j];
            }
            gradZij += 0.5 * p[L][etaY][etaX] * msetting.k[i][j];
            gradZij *= ht;
            g[gi++] = gradZij;
        }
    }

    // xi
    for (unsigned int j=0; j<msetting.Lo; j++)
    {
        g[gi++] = 0.0;
    }

    // eta
    for (unsigned int j=0; j<msetting.Lc; j++)
    {
        g[gi++] = 0.0;
    }
}

void AbstactProblem22D::setForward(IProblem2Forward2D *f)
{
    forward = f;
}

void AbstactProblem22D::setBackward(IProblem2Backward2D *b)
{
    backward = b;
}

const P2Setting &AbstactProblem22D::setting() const
{
    return msetting;
}

void AbstactProblem22D::setP2Setting(const P2Setting &setting)
{
    msetting = setting;
}
