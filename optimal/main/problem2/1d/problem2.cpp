#include "problem2.h"

void Problem2::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
        Dimension timeDimension(0.1, 0, 10);
        Dimension spaceDimensionX(0.01, 0, 100);

        IProblem2Forward ipf;
        ipf.setTimeDimension(timeDimension);
        ipf.addSpaceDimension(spaceDimensionX);

        DoubleMatrix u;
        ipf.gridMethod(u);

        IPrinter::printMatrix(14, 10, u);
        IPrinter::printSeperatorLine();

    //    IProblem2Backward ipb(u.row(10));
    //    ipb.setTimeDimension(timeDimension);
    //    ipb.addSpaceDimension(spaceDimensionX);

    //    DoubleMatrix p;
    //    ipb.gridMethod(p);

    //    IPrinter::printMatrix(14, 10, p);

//    Problem2 p;
//    p.Lo = 3;
//    p.Lc = 2;
//    p.k.resize(p.Lc, p.Lo);
//    p.z.resize(p.Lc, p.Lo);

//    p.xi.resize(p.Lo);
//    p.xi[0] = 0.25;
//    p.xi[1] = 0.50;
//    p.xi[2] = 0.75;

//    p.eta.resize(p.Lc);
//    p.eta[0] = 0.33;
//    p.eta[1] = 0.66;

//    for (unsigned int i=0; i<p.Lc; i++)
//    {
//        for (unsigned int j=0; j<p.Lo; j++)
//        {
//            p.k[i][j] = 5.0;//*(i+1)*(j-1);
//            p.z[i][j] = 8.0;//*(i+2)*(j-21);
//        }
//    }

//    p.a = 1.0;
//    p.lambda0 = 1.0;
//    p.lambda1 = 2.0;
//    p.lambda2 = 1.5;
//    p.theta = 20.0;

//    DoubleVector prms;
//    p.paremeters2Array(p.k, p.z, p.xi, p.eta, prms);

//    DoubleVector g1;
//    p.gradient(prms, g1);
//    //for (unsigned int k=0; k<15; k++) g1[k] = 0.0;
//    for (unsigned int k=15; k<g1.length(); k++) g1[k] = 0.0;
//    g1.L2Normalize();

//    IPrinter::print(g1,g1.length(),10,6);

//    DoubleVector g2(prms.length());
//    IGradient::Gradient(&p, 0.01, prms, g2);
//    //for (unsigned int k=0; k<15; k++) g2[k] = 0.0;
//    for (unsigned int k=15; k<g2.length(); k++) g2[k] = 0.0;
//    //g2.L2Normalize();
//    IPrinter::print(g2,g2.length(),10,6);
}

Problem2::Problem2()
{
    mTimeDimension = Dimension(0.01, 0, 100);
    mSpaceDimension = Dimension(0.01, 0, 100);

    forward.setTimeDimension(mTimeDimension);
    forward.addSpaceDimension(mSpaceDimension);

    backward.setTimeDimension(mTimeDimension);
    backward.addSpaceDimension(mSpaceDimension);

    unsigned int N = mSpaceDimension.sizeN();

    U.clear();
    U.resize(N+1);

    backward.U.clear();
    backward.U.resize(N+1);
    for (unsigned int n=0; n<=N; n++)
    {
        U[n] = 10.0;
        backward.U[n] = U[n];
    }

    backward.mu.clear();
    backward.mu.resize(N+1);
    for (unsigned int n=0; n<=N; n++) backward.mu[n] = 1.0;

    alpha0 = 1.0;
}

double Problem2::fx(const DoubleVector &prms) const
{
    DoubleMatrix k;
    DoubleMatrix z;
    DoubleVector xi;
    DoubleVector eta;

    array2Parameters(prms, k, z, xi, eta);

    DoubleMatrix u;
    Problem2Forward forward;
    forward.setTimeDimension(mTimeDimension);
    forward.addSpaceDimension(mSpaceDimension);
    forward.a = a; forward.lambda0 = lambda0; forward.lambda1 = lambda1; forward.lambda2 = lambda2; forward.theta = theta;
    forward.k = k; forward.z = z; forward.xi = xi; forward.eta = eta; forward.Lo = Lo; forward.Lc = Lc;
    forward.gridMethod(u);

    return alpha0*integral(u);
}

double Problem2::integral(const DoubleMatrix &u) const
{
    unsigned int N = mSpaceDimension.sizeN();
    unsigned int M = mTimeDimension.sizeN();
    double hx = mSpaceDimension.step();

    double sum = 0.0;
    sum += 0.5*mu(0)*(u[M][0]-U[0])*(u[M][0]-U[0]);
    for (unsigned int n=1; n<=N-1; n++)
    {
        sum += mu(n)*(u[M][n]-U[n])*(u[M][n]-U[n]);
    }
    sum += 0.5*mu(N)*(u[M][N]-U[N])*(u[M][N]-U[N]);

    sum *= hx;
    return sum;
}

double du(double x, unsigned int N, double hx, const DoubleMatrix &u, unsigned int m)
{
    unsigned int i = (unsigned int)(x*N);
    return (u[m][i+1] - u[m][i-1])/(2.0*hx);

//    double xm1 = (i-1)*hx;
//    double xp0 = (i+0)*hx;
//    double xp1 = (i+1)*hx;
//    double xp2 = (i+2)*hx;

//    return -1.0/(6.0*hx*hx*hx)*((x-xp1)*(x-xp2)+(x-xp0)*(x-xp2)+(x-xp0)*(x-xp1))*u[m][i-1]
//           +1.0/(2.0*hx*hx*hx)*((x-xp1)*(x-xp2)+(x-xm1)*(x-xp2)+(x-xm1)*(x-xp1))*u[m][i+0]
//           -1.0/(2.0*hx*hx*hx)*((x-xp0)*(x-xp2)+(x-xm1)*(x-xp2)+(x-xm1)*(x-xp0))*u[m][i+1]
//           +1.0/(6.0*hx*hx*hx)*((x-xp0)*(x-xp1)+(x-xm1)*(x-xp1)+(x-xm1)*(x-xp0))*u[m][i+2];
}

double dp(double x, unsigned int N, double hx, const DoubleMatrix &p, unsigned int m)
{
    unsigned int i = (unsigned int)(x*N);
    return (p[m][i+1] - p[m][i-1])/(2.0*hx);

//    double xm1 = (i-1)*hx;
//    double xp0 = (i+0)*hx;
//    double xp1 = (i+1)*hx;
//    double xp2 = (i+2)*hx;

//    return -1.0/(6.0*hx*hx*hx)*((x-xp1)*(x-xp2)+(x-xp0)*(x-xp2)+(x-xp0)*(x-xp1))*p[m][i-1]
//           +1.0/(2.0*hx*hx*hx)*((x-xp1)*(x-xp2)+(x-xm1)*(x-xp2)+(x-xm1)*(x-xp1))*p[m][i+0]
//           -1.0/(2.0*hx*hx*hx)*((x-xp0)*(x-xp2)+(x-xm1)*(x-xp2)+(x-xm1)*(x-xp0))*p[m][i+1]
//           +1.0/(6.0*hx*hx*hx)*((x-xp0)*(x-xp1)+(x-xm1)*(x-xp1)+(x-xm1)*(x-xp0))*p[m][i+2];
}

void Problem2::gradient(const DoubleVector &prms, DoubleVector &g)
{
    unsigned int N = mSpaceDimension.sizeN();
    unsigned int M = mTimeDimension.sizeN();
    double hx = mSpaceDimension.step();
    double ht = mTimeDimension.step();

    DoubleMatrix k;
    DoubleMatrix z;
    DoubleVector xi;
    DoubleVector eta;

    array2Parameters(prms, k, z, xi, eta);

    DoubleMatrix u;
    forward.a = a; forward.lambda0 = lambda0; forward.lambda1 = lambda1; forward.lambda2 = lambda2; forward.theta = theta;
    forward.k = k; forward.z = z; forward.xi = xi; forward.eta = eta; forward.Lo = Lo; forward.Lc = Lc;
    forward.gridMethod(u);
    //IPrinter::printMatrix(u);

    DoubleMatrix p;
    backward.uT = u.row(M);
    backward.a = a; backward.lambda0 = lambda0; backward.lambda1 = lambda1; backward.lambda2 = lambda2; backward.theta = theta;
    backward.k = k; backward.z = z; backward.xi = xi; backward.eta = eta; backward.Lo = Lo; backward.Lc = Lc;
    backward.gridMethod(p);

    g.resize(prms.length(), 0.0);

    for (unsigned int i=0; i<Lc; i++)
    {
        unsigned int etai = (unsigned int)(eta[i]*N);

        for (unsigned int j=0; j<Lo; j++)
        {
            unsigned int xij = (unsigned int)(xi[j]*N);

            double sumKij = 0.0;
            sumKij += 0.5*p[0][etai]*(u[0][xij] - z[i][j]);
            for (unsigned int m=1; m<=M-1; m++)
            {
                sumKij += p[m][etai]*(u[m][xij] - z[i][j]);
            }
            sumKij += 0.5*p[M][etai]*(u[M][xij] - z[i][j]);
            sumKij *= -ht;

            double sumZij = 0.0;
            sumZij += 0.5*p[0][etai];
            for (unsigned int m=1; m<=M-1; m++)
            {
                sumZij += p[m][etai];
            }
            sumZij += 0.5*p[M][etai];
            sumZij *= ht*k[i][j];

            g[i*Lo+j+0*Lc*Lo] = sumKij;
            g[i*Lo+j+1*Lc*Lo] = sumZij;
        }
    }

    for (unsigned int j=0; j<Lo; j++)
    {
        double SUM = 0.0;
        for (unsigned int i=0; i<Lc; i++)
        {
            unsigned int etai = (unsigned int)(eta[i]*N);

            double sumXij = 0.0;
            sumXij += 0.5*p[0][etai]*du(xi[j],N,hx,u,0);
            for (unsigned int m=1; m<=M-1; m++)
            {
                sumXij += p[m][etai]*du(xi[j],N,hx,u,m);
            }
            sumXij += 0.5*p[M][etai]*du(xi[j],N,hx,u,M);
            sumXij *= -ht*k[i][j];

            SUM += sumXij;
        }
        g[j+2*Lc*Lo] = SUM;
    }

    for (unsigned int i=0; i<Lc; i++)
    {
        //unsigned int etai = (unsigned int)(eta[i]*N);

        double SUM = 0.0;
        for (unsigned int j=0; j<Lc; j++)
        {
            unsigned int xij = (unsigned int)(xi[j]*N);

            double sumEtai = 0.0;
            sumEtai += 0.5*dp(eta[i],N,hx,p,0)*(u[0][xij] - z[i][j]);
            for (unsigned int m=1; m<=M-1; m++)
            {
                sumEtai += dp(eta[i],N,hx,p,m)*(u[m][xij] - z[i][j]);
            }
            sumEtai += 0.5*dp(eta[i],N,hx,p,M)*(u[M][xij] - z[i][j]);
            sumEtai *= -ht*k[i][j];

            SUM += sumEtai;
        }

        g[i+Lo+2*Lc*Lo] = SUM;
    }
}

void Problem2::array2Parameters(const DoubleVector &prms, DoubleMatrix &k, DoubleMatrix &z, DoubleVector &xi, DoubleVector &eta) const
{
    k.clear();
    k.resize(Lc, Lo);

    z.clear();
    z.resize(Lc, Lo);

    xi.clear();
    xi.resize(Lo);

    eta.clear();
    eta.resize(Lc);

    for (unsigned int i=0; i<Lc; i++)
    {
        for (unsigned int j=0; j<Lo; j++)
        {
            k[i][j] = prms[i*Lo+j+0*Lc*Lo];
            z[i][j] = prms[i*Lo+j+1*Lc*Lo];
        }
    }

    for (unsigned int j=0; j<Lo; j++)
    {
        xi[j] = prms[j+2*Lc*Lo];
    }

    for (unsigned int i=0; i<Lc; i++)
    {
        eta[i] = prms[i+Lo+2*Lc*Lo];
    }
}

void Problem2::paremeters2Array(const DoubleMatrix &k, const DoubleMatrix &z, const DoubleVector &xi, const DoubleVector &eta, DoubleVector &prms) const
{
    prms.clear();
    prms.resize(Lc*Lo + Lc*Lo + Lo + Lc);

    for (unsigned int i=0; i<Lc; i++)
    {
        for (unsigned int j=0; j<Lo; j++)
        {
            prms[i*Lo+j+0*Lc*Lo] = k[i][j];
            prms[i*Lo+j+1*Lc*Lo] = z[i][j];
        }
    }

    for (unsigned int j=0; j<Lo; j++)
    {
        prms[j+2*Lc*Lo] = xi[j];
    }

    for (unsigned int i=0; i<Lc; i++)
    {
        prms[i+Lo+2*Lc*Lo] = eta[i];
    }
}

double Problem2::mu(unsigned int n) const
{
    return 1.0;
}
