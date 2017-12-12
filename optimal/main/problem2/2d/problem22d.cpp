#include "problem22d.h"

void Problem22D::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    P2Setting setting1;
    //setting1.a = 1.0;
    //setting1.lambda = 0.01;
    //setting1.lambda0 = 0.1;
    //setting1.theta = 10.0;
    setting1.Lc = 2;
    setting1.Lo = 1;

    setting1.k.resize(setting1.Lc, setting1.Lo);
    setting1.z.resize(setting1.Lc, setting1.Lo);
    for (unsigned int i=0; i<setting1.Lc; i++)
    {
        for (unsigned int j=0; j<setting1.Lo; j++)
        {
            setting1.k[i][j] = -1.0;
            setting1.z[i][j] = +10.0;
        }
    }

    setting1.eta.resize(setting1.Lc);
    setting1.eta[0].x = 0.50; setting1.eta[0].y = 0.50;
    setting1.eta[1].x = 0.30; setting1.eta[1].y = 0.80;

    setting1.xi.resize(setting1.Lo);
    setting1.xi[0].x = 0.504;  setting1.xi[0].y = 0.502;
    //setting.xi[1].x = 0.50;  setting.xi[1].y = 0.50;
    //setting.xi[2].x = 0.65;  setting.xi[2].y = 0.65;

//    setting1.eta.resize(setting1.Lc);
//    setting1.eta[0].x = 0.33; setting1.eta[0].y = 0.33;
//    setting1.eta[1].x = 0.33; setting1.eta[1].y = 0.66;
//    setting1.eta[2].x = 0.66; setting1.eta[2].y = 0.66;
//    setting1.eta[3].x = 0.66; setting1.eta[3].y = 0.33;

//    setting1.xi.resize(setting1.Lo);
//    setting1.xi[0].x = 0.25;  setting1.xi[0].y = 0.25;
//    setting1.xi[1].x = 0.25;  setting1.xi[1].y = 0.75;
//    setting1.xi[2].x = 0.75;  setting1.xi[2].y = 0.75;
//    setting1.xi[3].x = 0.75;  setting1.xi[3].y = 0.25;
//    setting1.xi[4].x = 0.50;  setting1.xi[4].y = 0.50;

    Problem22D p22d;
    p22d.setting = setting1;
    p22d.setGridParameters(Dimension(0.01, 0, 100), Dimension(0.01, 0, 100), Dimension(0.01, 0, 100));
    //p22d.setGridParameters(Dimension(0.01, 0, 100), Dimension(0.010, 0, 100), Dimension(0.010, 0, 100));
    //p22d.setGridParameters(Dimension(0.01, 0, 100), Dimension(0.005, 0, 200), Dimension(0.010, 0, 100));
    //p22d.setGridParameters(Dimension(0.01, 0, 100), Dimension(0.010, 0, 100), Dimension(0.005, 0, 200));

    //clock_t t1 = clock();
    p22d.testForwardEquation(p22d.setting);
    //clock_t t2 = clock();
    //printf ("It took me %d clicks (%f seconds).\n",t2-t1,((float)(t2-t1))/CLOCKS_PER_SEC);
    //p22d.testBackwardEquation(p22d.setting);
    //clock_t t3 = clock();
    //printf ("It took me %d clicks (%f seconds).\n",t3-t2,((float)(t3-t2))/CLOCKS_PER_SEC);
    return;

    unsigned int N1 = p22d.mSpaceDimensionX.sizeN();
    unsigned int N2 = p22d.mSpaceDimensionY.sizeN();

    p22d.U.resize(N2+1, N1+1, 10.0);

    //p22d.backward.U = p22d.U;
    //p22d.backward.uT = u[u.size()-1];

    //std::vector<DoubleMatrix> p;
    //p22d.backward.calculateMVD(p);
    //IPrinter::printMatrix(p[0]);
    //IPrinter::printSeperatorLine();

    DoubleVector prm;
    DoubleVector agrd;
    DoubleVector ngrd;

    p22d.setting.toVector(prm);

    agrd.resize(prm.length(), 0.0);
    ngrd.resize(prm.length(), 0.0);

    p22d.gradient(prm, agrd);

    for (unsigned int i=0; i<p22d.setting.Lc; i++)
    {
        for (unsigned int j=0; j<p22d.setting.Lo; j++)
        {
            unsigned int index = i*p22d.setting.Lo + j;
            DoubleVector x2 = prm; x2[index] += 0.01; double f2 = p22d.fx(x2);
            DoubleVector x1 = prm; x1[index] -= 0.01; double f1 = p22d.fx(x1);
            ngrd[index] = (f2 - f1)/0.02;
        }
    }

    for (unsigned int i=0; i<setting1.Lc; i++)
    {
        for (unsigned int j=0; j<setting1.Lo; j++)
        {
            unsigned int index = p22d.setting.Lc*p22d.setting.Lo + i*p22d.setting.Lo + j;
            DoubleVector x2 = prm; x2[index] += 0.01; double f2 = p22d.fx(x2);
            DoubleVector x1 = prm; x1[index] -= 0.01; double f1 = p22d.fx(x1);
            ngrd[index] = (f2 - f1)/0.02;
        }
    }

    FILE *file = fopen("gradients.txt", "w");
    IPrinter::print(prm,prm.length(),14, 10, file);
    agrd.L2Normalize();
    IPrinter::print(agrd,agrd.length(),14, 10, file);
    ngrd.L2Normalize();
    IPrinter::print(ngrd,agrd.length(),14, 10, file);
    fclose(file);
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
    DoubleMatrix u;
    vector<ExtendedSpaceNode2D> info;
    forward.calculateMVD(u, info);

    double intgrl = integral(u);

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

    setting.fromVector(prms);

//    Problem2Forward2D forward;
//    forward.setTimeDimension(mTimeDimension);
//    forward.addSpaceDimension(mSpaceDimensionX);
//    forward.addSpaceDimension(mSpaceDimensionY);
//    forward.setSettings(setting);
//    std::vector<DoubleMatrix> u;
//    forward.calculateMVD(u);

//    Problem2Backward2D backward;
//    backward.setTimeDimension(mTimeDimension);
//    backward.addSpaceDimension(mSpaceDimensionX);
//    backward.addSpaceDimension(mSpaceDimensionY);
//    backward.setSettings(setting);

//    backward.U = U;
//    backward.uT = u[u.size()-1];
//    backward.mu.resize(Ny+1, Nx+1, 1.0);

//    std::vector<DoubleMatrix> p;
//    backward.calculateMVD(p);

//    g.resize(prms.length(), 0.0);
//    unsigned int gi = 0;

//    // k
//    for (unsigned int i=0; i<setting.Lc; i++)
//    {
//        const SpaceNodePDE &eta = setting.eta[i];
//        unsigned int etaX = (unsigned int)(round(eta.x*Nx));
//        unsigned int etaY = (unsigned int)(round(eta.y*Ny));

//        for (unsigned int j=0; j<setting.Lo; j++)
//        {
//            const SpaceNodePDE &xi = setting.xi[j];
//            unsigned int xiX = (unsigned int)(round(xi.x*Nx));
//            unsigned int xiY = (unsigned int)(round(xi.y*Ny));

//            double gradKij = 0.0;
//            gradKij += 0.5 * p[0][etaY][etaX] * (u[0][xiY][xiX] - setting.z[i][j]);
//            for (unsigned int m=1; m<=M-1; m++)
//            {
//                gradKij += p[m][etaY][etaX] * (u[m][xiY][xiX] - setting.z[i][j]);
//            }
//            gradKij += 0.5 * p[M][etaY][etaX] * (u[M][xiY][xiX] - setting.z[i][j]);
//            gradKij *= -ht;
//            g[gi++] = gradKij;
//        }
//    }

//    // z
//    for (unsigned int i=0; i<setting.Lc; i++)
//    {
//        const SpaceNodePDE &eta = setting.eta[i];
//        unsigned int etaX = (unsigned int)(round(eta.x*Nx));
//        unsigned int etaY = (unsigned int)(round(eta.y*Ny));

//        for (unsigned int j=0; j<setting.Lo; j++)
//        {
//            double gradZij = 0.0;
//            gradZij += 0.5 * p[0][etaY][etaX] * setting.k[i][j];
//            for (unsigned int m=1; m<=M-1; m++)
//            {
//                gradZij += p[m][etaY][etaX]  * setting.k[i][j];
//            }
//            gradZij += 0.5 * p[M][etaY][etaX] * setting.k[i][j];
//            gradZij *= ht;
//            g[gi++] = gradZij;
//        }
//    }

//    // xi
//    for (unsigned int j=0; j<setting.Lo; j++)
//    {
//        g[gi++] = 0.0;
//    }

//    // eta
//    for (unsigned int j=0; j<setting.Lc; j++)
//    {
//        g[gi++] = 0.0;
//    }
}

double Problem22D::integral(const DoubleMatrix &u) const
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

double Problem22D::mu(double x UNUSED_PARAM, double y UNUSED_PARAM) const
{
    return 1.0;
}

void Problem22D::testForwardEquation(const P2Setting &setting) const
{
    CProblem2Forward2D forward;
    forward.a = 1.0;
    forward.lambda = 0.01;
    forward.lambda0 = 0.1;
    forward.theta = 10.0;

    forward.setTimeDimension(mTimeDimension);
    forward.addSpaceDimension(mSpaceDimensionX);
    forward.addSpaceDimension(mSpaceDimensionY);
    forward.setSettings(setting);

    DoubleMatrix u;
    vector<ExtendedSpaceNode2D> info;
    forward.calculateMVD(u, info);
    IPrinter::printMatrix(u);
    IPrinter::printSeperatorLine();
    u.clear();

    for (unsigned int i=0; i<info.size(); i++)
    {
        const ExtendedSpaceNode2D &esn = info[i];
        //for (unsigned int ln=0; ln<esn.layerNumber; ln++)
        unsigned int ln = esn.layerNumber;
        {
            printf("%12.8f %12.8f %12.8f %12.8f\n%12.8f %12.8f %12.8f %12.8f\n%12.8f %12.8f %12.8f %12.8f\n%12.8f %12.8f %12.8f %12.8f\n---\n",
                    esn.wi[3][0].u[ln],esn.wi[3][1].u[ln],esn.wi[3][2].u[ln],esn.wi[3][3].u[ln],
                    esn.wi[2][0].u[ln],esn.wi[2][1].u[ln],esn.wi[2][2].u[ln],esn.wi[2][3].u[ln],
                    esn.wi[1][0].u[ln],esn.wi[1][1].u[ln],esn.wi[1][2].u[ln],esn.wi[1][3].u[ln],
                    esn.wi[0][0].u[ln],esn.wi[0][1].u[ln],esn.wi[0][2].u[ln],esn.wi[0][3].u[ln]);
        }

        printf("%12.8f %12.8f %12.8f %12.8f\n%12.8f %12.8f %12.8f %12.8f\n%12.8f %12.8f %12.8f %12.8f\n%12.8f %12.8f %12.8f %12.8f\n---\n",
                esn.wi[3][0].w,esn.wi[3][1].w,esn.wi[3][2].w,esn.wi[3][3].w,
                esn.wi[2][0].w,esn.wi[2][1].w,esn.wi[2][2].w,esn.wi[2][3].w,
                esn.wi[1][0].w,esn.wi[1][1].w,esn.wi[1][2].w,esn.wi[1][3].w,
                esn.wi[0][0].w,esn.wi[0][1].w,esn.wi[0][2].w,esn.wi[0][3].w);

        printf("%f\n", esn.value(esn.x, esn.y, ln));
    }
}

void Problem22D::testBackwardEquation(const P2Setting &setting) const
{
    CProblem2Backward2D backward;
    backward.a = 1.0;
    backward.lambda = 0.01;
    backward.lambda0 = 0.1;
    backward.theta = 10.0;

    backward.setTimeDimension(mTimeDimension);
    backward.addSpaceDimension(mSpaceDimensionX);
    backward.addSpaceDimension(mSpaceDimensionY);
    backward.setSettings(setting);

    backward.U.resize(mSpaceDimensionY.sizeN()+1, mSpaceDimensionY.sizeN()+1, 10.0);
    backward.uT.resize(mSpaceDimensionY.sizeN()+1, mSpaceDimensionY.sizeN()+1, 9.0);
    backward.mu.resize(mSpaceDimensionY.sizeN()+1, mSpaceDimensionY.sizeN()+1, 1.0);

    DoubleMatrix p;
    backward.calculateMVD(p);
    IPrinter::printMatrix(p);
    IPrinter::printSeperatorLine();
    p.clear();
}

void Problem22D::testExample1()
{
    Problem22D p22d;
}



