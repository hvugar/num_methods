#include "problem22d.h"

void Problem22D::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    P2Setting setting;
    setting.a = 1.0;
    setting.lambda = 0.01;
    setting.lambda0 = 1.0;
    setting.theta = 10.0;
    setting.Lc = 2;
    setting.Lo = 1;

    setting.k.resize(setting.Lc, setting.Lo);
    setting.z.resize(setting.Lc, setting.Lo);
    for (unsigned int i=0; i<setting.Lc; i++)
    {
        for (unsigned int j=0; j<setting.Lo; j++)
        {
            setting.k[i][j] = +5.0;
            setting.z[i][j] = +4.0;
        }
    }

    setting.eta.resize(setting.Lc);
    setting.eta[0].x = 0.25; setting.eta[0].y = 0.25;
    setting.eta[1].x = 0.75; setting.eta[1].y = 0.75;

    setting.xi.resize(setting.Lo);
    setting.xi[0].x = 0.65;  setting.xi[0].y = 0.65;
//    setting.xi[1].x = 0.50;  setting.xi[1].y = 0.50;
//    setting.xi[2].x = 0.65;  setting.xi[2].y = 0.65;


//    setting.eta.resize(setting.Lc);
//    setting.eta[0].x = 0.33; setting.eta[0].y = 0.33;
//    setting.eta[1].x = 0.33; setting.eta[1].y = 0.66;
//    setting.eta[2].x = 0.66; setting.eta[2].y = 0.66;
//    setting.eta[3].x = 0.66; setting.eta[3].y = 0.33;

//    setting.xi.resize(setting.Lo);
//    setting.xi[0].x = 0.25;  setting.xi[0].y = 0.25;
//    setting.xi[1].x = 0.25;  setting.xi[1].y = 0.75;
//    setting.xi[2].x = 0.75;  setting.xi[2].y = 0.75;
//    setting.xi[3].x = 0.75;  setting.xi[3].y = 0.25;
//    setting.xi[4].x = 0.50;  setting.xi[4].y = 0.50;


    Problem22D p22d;
    p22d.setGridParameters(Dimension(0.01, 0, 100), Dimension(0.01, 0, 100), Dimension(0.01, 0, 100));
    unsigned int N1 = p22d.mSpaceDimensionX.sizeN();
    unsigned int N2 = p22d.mSpaceDimensionY.sizeN();
    p22d.U.resize(N2+1, N1+1, 10.0);

//    std::vector<DoubleMatrix> u;
//    p22d.forward.calculateMVD(u);
//    IPrinter::printMatrix(u[u.size()-1]);
//    IPrinter::printSeperatorLine();

//    p22d.backward.U = p22d.U;
//    p22d.backward.uT = u[u.size()-1];

//    std::vector<DoubleMatrix> p;
//    p22d.backward.calculateMVD(p);
//    IPrinter::printMatrix(p[0]);
//    IPrinter::printSeperatorLine();

    DoubleVector prm;
    DoubleVector agrd;
    DoubleVector ngrd;

    setting.toVector(prm);

    agrd.resize(prm.length(), 0.0);
    ngrd.resize(prm.length(), 0.0);

    p22d.gradient(prm, agrd);

    for (unsigned int i=0; i<setting.Lc; i++)
    {
        for (unsigned int j=0; j<setting.Lo; j++)
        {
            unsigned int index = i*setting.Lo + j;
            DoubleVector x2 = prm; x2[index] += 0.01; double f2 = p22d.fx(x2);
            DoubleVector x1 = prm; x1[index] -= 0.01; double f1 = p22d.fx(x1);
            ngrd[index] = (f2 - f1)/0.02;
        }
    }

    for (unsigned int i=0; i<p22d.setting.Lc; i++)
    {
        for (unsigned int j=0; j<p22d.setting.Lo; j++)
        {
            unsigned int index = p22d.setting.Lc*p22d.setting.Lo + i*p22d.setting.Lo + j;
            DoubleVector x2 = prm; x2[index] += 0.01; double f2 = p22d.fx(x2);
            DoubleVector x1 = prm; x1[index] -= 0.01; double f1 = p22d.fx(x1);
            ngrd[index] = (f2 - f1)/0.02;
        }
    }

    IPrinter::print(prm);
    agrd.L2Normalize();
    IPrinter::print(agrd);
    ngrd.L2Normalize();
    IPrinter::print(ngrd);
}

void Problem22D::setGridParameters(Dimension timeDimension, Dimension spaceDimensionX, Dimension spaceDimensionY)
{
    mTimeDimension = timeDimension;
    mSpaceDimensionX = spaceDimensionX;
    mSpaceDimensionY = spaceDimensionY;
}

double Problem22D::fx(const DoubleVector &prms) const
{
    P2Setting setting1 = setting;
    setting1.fromVector(prms);

    Problem2Forward2D forward;
    forward.setTimeDimension(mTimeDimension);
    forward.addSpaceDimension(mSpaceDimensionX);
    forward.addSpaceDimension(mSpaceDimensionY);
    forward.setSettings(setting1);
    std::vector<DoubleMatrix> u;
    forward.calculateMVD(u);

    double intgrl = integral(u[u.size()-1]);

    return intgrl;
}

void Problem22D::gradient(const DoubleVector &prms UNUSED_PARAM, DoubleVector &g UNUSED_PARAM)
{
    unsigned int Nx = mSpaceDimensionX.sizeN();
    unsigned int Ny = mSpaceDimensionY.sizeN();
    unsigned int M = mTimeDimension.sizeN();

    double hx = mSpaceDimensionX.step();
    double hy = mSpaceDimensionY.step();
    double ht = mTimeDimension.step();

    P2Setting setting1 = setting;
    setting1.fromVector(prms);


    Problem2Forward2D forward;
    forward.setTimeDimension(mTimeDimension);
    forward.addSpaceDimension(mSpaceDimensionX);
    forward.addSpaceDimension(mSpaceDimensionY);
    forward.setSettings(setting1);
    std::vector<DoubleMatrix> u;
    forward.calculateMVD(u);

    Problem2Backward2D backward;
    backward.setTimeDimension(mTimeDimension);
    backward.addSpaceDimension(mSpaceDimensionX);
    backward.addSpaceDimension(mSpaceDimensionY);
    backward.setSettings(setting1);
    std::vector<DoubleMatrix> p;
    backward.calculateMVD(p);

    g.resize(prms.length(), 0.0);
    unsigned int gi = 0;

    // k
    for (unsigned int i=0; i<setting1.Lc; i++)
    {
        const SpaceNodePDE &eta = setting1.eta[i];
        unsigned int etaX = (unsigned int)(round(eta.x*Nx));
        unsigned int etaY = (unsigned int)(round(eta.y*Ny));

        for (unsigned int j=0; j<setting1.Lo; j++)
        {
            const SpaceNodePDE &xi = setting1.xi[j];
            unsigned int xiX = (unsigned int)(round(xi.x*Nx));
            unsigned int xiY = (unsigned int)(round(xi.y*Ny));

            double gradKij = 0.0;
            gradKij += 0.5 * p[0][etaY][etaX] * (u[0][xiY][xiX] - setting1.z[i][j]);
            for (unsigned int m=1; m<=M-1; m++)
            {
                gradKij += p[m][etaY][etaX] * (u[m][xiY][xiX] - setting1.z[i][j]);
            }
            gradKij += 0.5 * p[M][etaY][etaX] * (u[M][xiY][xiX] - setting1.z[i][j]);
            gradKij *= -ht;
            g[gi++] = gradKij;
        }
    }

    // z
    for (unsigned int i=0; i<setting1.Lc; i++)
    {
        const SpaceNodePDE &eta = setting1.eta[i];
        unsigned int etaX = (unsigned int)(round(eta.x*Nx));
        unsigned int etaY = (unsigned int)(round(eta.y*Ny));

        for (unsigned int j=0; j<setting1.Lo; j++)
        {
            double gradZij = 0.0;
            gradZij += 0.5 * p[0][etaY][etaX] * setting1.k[i][j];
            for (unsigned int m=1; m<=M-1; m++)
            {
                gradZij += p[m][etaY][etaX]  * setting1.k[i][j];
            }
            gradZij += 0.5 * p[M][etaY][etaX] * setting1.k[i][j];
            gradZij *= ht;
            g[gi++] = gradZij;
        }
    }

    // xi
    for (unsigned int j=0; j<setting1.Lo; j++)
    {
        g[gi++] = 0.0;
    }

    // eta
    for (unsigned int j=0; j<setting1.Lc; j++)
    {
        g[gi++] = 0.0;
    }
}

double Problem22D::integral(const DoubleMatrix &u) const
{
    double hx = mSpaceDimensionX.step();
    double hy = mSpaceDimensionY.step();
    unsigned int N1 = mSpaceDimensionX.sizeN();
    unsigned int N2 = mSpaceDimensionY.sizeN();

    double sum = 0.0;

    sum += 0.25* (u[0][0]  - U[0][0])   * (u[0][0]   - U[0][0]);
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

double Problem22D::mu(double x UNUSED_PARAM, double y UNUSED_PARAM) const
{
    return 1.0;
}


