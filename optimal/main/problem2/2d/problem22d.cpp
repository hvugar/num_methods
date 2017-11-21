#include "problem22d.h"

void Problem22D::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    Problem22D p22d(1.0, 0.1, 1.0, 10.0, 4, 5);
    p22d.setGridParameters(Dimension(0.01, 0, 100), Dimension(0.01, 0, 100), Dimension(0.01, 0, 100));

    std::vector<DoubleMatrix> u;
    p22d.forward.calculateMVD(u);
    IPrinter::printMatrix(u[u.size()-1]);
    IPrinter::printSeperatorLine();

    p22d.backward.U = p22d.U;
    p22d.backward.uT = u[u.size()-1];

    DoubleMatrix p;
    p22d.backward.calculateMVD(p);
    IPrinter::printMatrix(p);
    IPrinter::printSeperatorLine();
}

Problem22D::Problem22D(double a, double lambda0, double lambda, double theta, double Lc, double Lo)
{
    setting.a = a;
    setting.lambda = lambda;
    setting.lambda0 = lambda0;
    setting.theta = theta;
    setting.Lc = Lc;
    setting.Lo = Lo;

    setting.k.resize(setting.Lc, setting.Lo);
    setting.z.resize(setting.Lc, setting.Lo);
    for (unsigned int i=0; i<setting.Lc; i++)
    {
        for (unsigned int j=0; j<setting.Lo; j++)
        {
            setting.k[i][j] = -2.0;
            setting.z[i][j] = 10.0;
        }
    }

    setting.eta.resize(setting.Lc);
    setting.eta[0].x = 0.33; setting.eta[0].y = 0.33;
    setting.eta[1].x = 0.33; setting.eta[1].y = 0.66;
    setting.eta[2].x = 0.66; setting.eta[2].y = 0.66;
    setting.eta[3].x = 0.66; setting.eta[3].y = 0.33;

    setting.xi.resize(setting.Lo);
    setting.xi[0].x = 0.25;  setting.xi[0].y = 0.25;
    setting.xi[1].x = 0.25;  setting.xi[1].y = 0.75;
    setting.xi[2].x = 0.75;  setting.xi[2].y = 0.75;
    setting.xi[3].x = 0.75;  setting.xi[3].y = 0.25;
    setting.xi[4].x = 0.50;  setting.xi[4].y = 0.50;

    forward.setSettings(setting);
    backward.setSettings(setting);
}

void Problem22D::setGridParameters(Dimension timeDimension, Dimension spaceDimensionX, Dimension spaceDimensionY)
{
    mTimeDimension = timeDimension;
    mSpaceDimensionX = spaceDimensionX;
    mSpaceDimensionY = spaceDimensionY;

    forward.setTimeDimension(mTimeDimension);
    forward.addSpaceDimension(mSpaceDimensionX);
    forward.addSpaceDimension(mSpaceDimensionY);

    backward.setTimeDimension(mTimeDimension);
    backward.addSpaceDimension(mSpaceDimensionX);
    backward.addSpaceDimension(mSpaceDimensionY);

    unsigned int N1 = mSpaceDimensionX.sizeN();
    unsigned int N2 = mSpaceDimensionY.sizeN();

    U.resize(N2+1, N1+1, 10.0);
}

double Problem22D::fx(const DoubleVector &prms) const
{
    P2Setting setting;
    array2Parameters(prms, setting.k, setting.z, setting.xi, setting.eta);

    Problem2Forward2D forward;
    forward.setTimeDimension(mTimeDimension);
    forward.addSpaceDimension(mSpaceDimensionX);
    forward.addSpaceDimension(mSpaceDimensionY);
    forward.setSettings(setting);
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

    P2Setting setting;
    array2Parameters(prms, setting.k, setting.z, setting.xi, setting.eta);

    Problem2Forward2D forward;
    forward.setTimeDimension(mTimeDimension);
    forward.addSpaceDimension(mSpaceDimensionX);
    forward.addSpaceDimension(mSpaceDimensionY);
    forward.setSettings(setting);
    std::vector<DoubleMatrix> u;
    forward.calculateMVD(u);

    Problem2Backward2D backward;
    backward.setTimeDimension(mTimeDimension);
    backward.addSpaceDimension(mSpaceDimensionX);
    backward.addSpaceDimension(mSpaceDimensionY);
    backward.setSettings(setting);
    DoubleMatrix p;
    backward.calculateMVD(p);

    g.resize(prms.length(), 0.0);

    for (unsigned int i=0; i<setting.Lc; i++)
    {
        const SpaceNodePDE &eta = setting.eta[i];
        unsigned int etaX = (unsigned int)(round(eta.x*Nx));
        unsigned int etaY = (unsigned int)(round(eta.y*Ny));

        for (unsigned int j=0; j<setting.Lo; j++)
        {
            const SpaceNodePDE &xi = setting.xi[j];
            unsigned int xiX = (unsigned int)(round(xi.x*Nx));
            unsigned int xiY = (unsigned int)(round(xi.y*Ny));

//            double sumKij = 0.0;
//            sumKij += 0.5*p[0][etai]*(u[0][xij] - z[i][j]);
//            for (unsigned int m=1; m<=M-1; m++)
//            {
//                sumKij += p[m][etai]*(u[m][xij] - z[i][j]);
//            }
//            sumKij += 0.5*p[M][etai]*(u[M][xij] - z[i][j]);
//            sumKij *= -ht;

//            double sumZij = 0.0;
//            sumZij += 0.5*p[0][etai];
//            for (unsigned int m=1; m<=M-1; m++)
//            {
//                sumZij += p[m][etai];
//            }
//            sumZij += 0.5*p[M][etai];
//            sumZij *= ht*k[i][j];

//            g[i*Lo+j+0*Lc*Lo] = sumKij;
//            g[i*Lo+j+1*Lc*Lo] = sumZij;
        }
    }

//    for (unsigned int j=0; j<Lo; j++)
//    {
//        double SUM = 0.0;
//        for (unsigned int i=0; i<Lc; i++)
//        {
//            unsigned int etai = (unsigned int)(eta[i]*N);

//            double sumXij = 0.0;
//            sumXij += 0.5*p[0][etai]*du(xi[j],N,hx,u,0);
//            for (unsigned int m=1; m<=M-1; m++)
//            {
//                sumXij += p[m][etai]*du(xi[j],N,hx,u,m);
//            }
//            sumXij += 0.5*p[M][etai]*du(xi[j],N,hx,u,M);
//            sumXij *= -ht*k[i][j];

//            SUM += sumXij;
//        }
//        g[j+2*Lc*Lo] = SUM;
//    }

//    for (unsigned int i=0; i<Lc; i++)
//    {
//        //unsigned int etai = (unsigned int)(eta[i]*N);

//        double SUM = 0.0;
//        for (unsigned int j=0; j<Lc; j++)
//        {
//            unsigned int xij = (unsigned int)(xi[j]*N);

//            double sumEtai = 0.0;
//            sumEtai += 0.5*dp(eta[i],N,hx,p,0)*(u[0][xij] - z[i][j]);
//            for (unsigned int m=1; m<=M-1; m++)
//            {
//                sumEtai += dp(eta[i],N,hx,p,m)*(u[m][xij] - z[i][j]);
//            }
//            sumEtai += 0.5*dp(eta[i],N,hx,p,M)*(u[M][xij] - z[i][j]);
//            sumEtai *= -ht*k[i][j];

//            SUM += sumEtai;
//        }

//        g[i+Lo+2*Lc*Lo] = SUM;
//    }
}

double Problem22D::integral(const DoubleMatrix &u) const
{
    double hx = mSpaceDimensionX.step();
    double hy = mSpaceDimensionY.step();
    unsigned int N1 = mSpaceDimensionX.sizeN();
    unsigned int N2 = mSpaceDimensionY.sizeN();

    double sum = 0.0;

    sum += 0.25* (u[0][0] - U[0][0])*(u[0][0] - U[0][0]);
    sum += 0.25*(u[0][N1] - U[0][N1])*(u[0][N1] - U[0][N1]);
    sum += 0.25*(u[N2][0] - U[N2][0])*(u[N2][0] - U[N2][0]);
    sum += 0.25*(u[N2][N1] - U[N2][N1])*(u[N2][N1] - U[N2][N1]);

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

void Problem22D::array2Parameters(const DoubleVector &prms, DoubleMatrix &k, DoubleMatrix &z, std::vector<SpaceNodePDE> &xi, std::vector<SpaceNodePDE> &eta) const
{
    k.clear();
    k.resize(setting.Lc, setting.Lo);

    z.clear();
    z.resize(setting.Lc, setting.Lo);

    xi.clear();
    xi.resize(setting.Lo);

    eta.clear();
    eta.resize(setting.Lc);

    for (unsigned int i=0; i<setting.Lc; i++)
    {
        for (unsigned int j=0; j<setting.Lo; j++)
        {
            k[i][j] = prms[i*setting.Lo+j+0*setting.Lc*setting.Lo];
            z[i][j] = prms[i*setting.Lo+j+1*setting.Lc*setting.Lo];
        }
    }

    for (unsigned int j=0; j<setting.Lo; j++)
    {
        xi[j].x = prms[2*j + 0 + 2*setting.Lc*setting.Lo];
        xi[j].y = prms[2*j + 1 + 2*setting.Lc*setting.Lo];
    }

    for (unsigned int i=0; i<setting.Lc; i++)
    {
        eta[i].x = prms[2*i + 0 + 2*setting.Lo + 2*setting.Lc*setting.Lo];
        eta[i].y = prms[2*i + 1 + 2*setting.Lo + 2*setting.Lc*setting.Lo];
    }
}

void Problem22D::paremeters2Array(const DoubleMatrix &k, const DoubleMatrix &z, const std::vector<SpaceNodePDE> &xi, const std::vector<SpaceNodePDE> &eta, DoubleVector &prms) const
{
    prms.clear();
    prms.resize(setting.Lc*setting.Lo + setting.Lc*setting.Lo + 2*setting.Lo + 2*setting.Lc);

    for (unsigned int i=0; i<setting.Lc; i++)
    {
        for (unsigned int j=0; j<setting.Lo; j++)
        {
            prms[i*setting.Lo + j + 0*setting.Lc*setting.Lo] = k[i][j];
            prms[i*setting.Lo + j + 1*setting.Lc*setting.Lo] = z[i][j];
        }
    }

    for (unsigned int j=0; j<setting.Lo; j++)
    {
        prms[2*j + 0 + 2*setting.Lc*setting.Lo] = xi[j].x;
        prms[2*j + 1 + 2*setting.Lc*setting.Lo] = xi[j].y;
    }

    for (unsigned int i=0; i<setting.Lc; i++)
    {
        prms[2*i + 0 + 2*setting.Lo + 2*setting.Lc*setting.Lo] = eta[i].x;
        prms[2*i + 1 + 2*setting.Lo + 2*setting.Lc*setting.Lo] = eta[i].y;
    }
}

double Problem22D::mu(double x UNUSED_PARAM, double y UNUSED_PARAM) const
{
    return 1.0;
}


