#include "iproblem2forward2d.h"

IProblem2Forward2D::IProblem2Forward2D(double a, double lambda0, double lambda, double theta, unsigned int Lc1, unsigned int Lo1)
{
    this->a = a;
    this->lambda0 = lambda0;
    this->lambda = lambda;
    this->theta = theta;

    this->Lc = Lc1;
    this->Lo = Lo1;

    this->Lo = 5;
    this->Lc = 4;

    k.resize(Lc, Lo);
    z.resize(Lc, Lo);

    for (unsigned int i=0; i<this->Lc; i++)
    {
        for (unsigned int j=0; j<this->Lo; j++)
        {
            k[i][j] = 5.0;
            z[i][j] = 8.0;
        }
    }

    xi.resize(this->Lo);
    xi[0].x = 0.2543;  xi[0].y = 0.2554; //xi[0].i = 25; xi[0].j = 25;
    xi[1].x = 0.2566;  xi[1].y = 0.7556; //xi[1].i = 25; xi[1].j = 75;
    xi[2].x = 0.7532;  xi[2].y = 0.7576; //xi[2].i = 75; xi[2].j = 75;
    xi[3].x = 0.7511;  xi[3].y = 0.2577; //xi[3].i = 75; xi[3].j = 25;
    xi[4].x = 0.5098;  xi[4].y = 0.5086; //xi[4].i = 50; xi[4].j = 50;

    eta.resize(this->Lc);
    eta[0].x = 0.3345; eta[0].y = 0.3318; //eta[0].i = 33; eta[0].j = 33;
    eta[1].x = 0.3322; eta[1].y = 0.6675; //eta[1].i = 33; eta[1].j = 66;
    eta[2].x = 0.6634; eta[2].y = 0.6665; //eta[2].i = 66; eta[2].j = 66;
    eta[3].x = 0.6652; eta[3].y = 0.3356; //eta[3].i = 66; eta[3].j = 33;
}

void IProblem2Forward2D::calculateMVD(DoubleMatrix &u)
{
    Dimension xd = spaceDimension(Dimension::DimensionX);
    Dimension yd = spaceDimension(Dimension::DimensionY);
    Dimension td = timeDimension();

    unsigned int N = xd.sizeN();
    unsigned int M = yd.sizeN();
    unsigned int L = td.sizeN();
    double hx = xd.step();
    double hy = yd.step();
    double ht = td.step();

    u.clear();
    u.resize(M+1, N+1);

    DoubleMatrix uh(M+1, N+1);

    SpaceNodePDE sn;
    for (unsigned int m=0; m<=M; m++)
    {
        sn.j = m; sn.y = m*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = n; sn.x = n*hx;
            u[m][n] = initial(sn);
        }
    }

    IPrinter::printMatrix(18, 14, u);
    IPrinter::printSeperatorLine();

    TimeNodePDE tn;
    for (unsigned int l=1; l<=1; l++)
    {
        /////////////////////////////////////////////////////////////////////////////////////////

        tn.i = l;
        tn.t = l*ht - 0.5*ht;

        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = n; sn.x = n*hx;

            DoubleMatrix w1(M+1, M+1, 0.0);
            DoubleVector d1(M+1, 0.0);

            for (unsigned int m=0; m<=M; m++)
            {
                sn.j = m; sn.y = m*hy;

                d1[m] = 2.0*u[m][n] + lambda0*theta*ht + ht*f(sn, tn);

                if (n==0)       d1[m] += ((a*a*ht)/(hx*hx))*(u[m][0]   - 2.0*u[m][1]   + u[m][2]);
                if (n>0 && n<N) d1[m] += ((a*a*ht)/(hx*hx))*(u[m][n-1] - 2.0*u[m][n]   + u[m][n+1]);
                if (n==N)       d1[m] += ((a*a*ht)/(hx*hx))*(u[m][N-2] - 2.0*u[m][N-1] + u[m][N]);

                if (m == 0)
                {
                    w1[0][0] = 2.0 + (2.0*a*a*ht)/(hy*hy) + lambda0*ht + (2.0*a*a*lambda*ht)/(hy);
                    w1[0][1] = -(2.0*a*a*ht)/(hy*hy);

                    d1[0] += (2.0*a*a*lambda*theta*ht)/(hy) + ((2.0*a*a*ht)/(hy))*g3(sn, tn);
                }
                else if (m == M)
                {
                    w1[M][M-1] = -(2.0*a*a*ht)/(hy*hy);
                    w1[M][M-0] = 2.0 + (2.0*a*a*ht)/(hy*hy) + lambda0*ht + (2.0*a*a*lambda*ht)/(hy);

                    d1[M] += (2.0*a*a*lambda*theta*ht)/(hy) + ((2.0*a*a*ht)/(hy))*g4(sn, tn);
                }
                else
                {
                    w1[m][m-1] = -(a*a*ht)/(hy*hy);
                    w1[m][m+0] = 2.0 + (2.0*a*a*ht)/(hy*hy) + lambda0*ht;
                    w1[m][m+1] = -(a*a*ht)/(hy*hy);
                }

                // Adding delta part ***********************************************************************************

                for (unsigned int i=0; i<Lc; i++)
                {
                    double _delta = delta(sn, i);

                    for (unsigned int j=0; j<Lo; j++)
                    {
                        d1[m] += -ht*k[i][j]*z[i][j] * _delta;

                        //unsigned int jinx = (unsigned int)(xi[j].x*N);
                        unsigned int jiny = (unsigned int)(xi[j].y*M);

                        //w1[m][jiny] += -ht*k[i][j];

                        double hx3 = hx*hx*hx;
                        double hx32 = (1.0/(2.0*hx3));
                        double hx36 = (1.0/(6.0*hx3));

                        double hy3 = hy*hy*hy;
                        double hy32 = (1.0/(2.0*hy3));
                        double hy36 = (1.0/(6.0*hy3));

                        for (unsigned int m1=0; m1<=M; m1++)
                        {

                            double dh = fabs(m1*hy - xi[j].y);

                            if (dh <= hy)
                            {
                                w1[m][m1] += -ht*k[i][j] * ((2.0*hy-dh)*(hy-dh)*(hy+dh)) * hy32 * _delta;
                            }

                            if (hy < dh && dh <= 2.0*hy)
                            {
                                w1[m][m1] += -ht*k[i][j] * ((2.0*hy-dh)*(hy-dh)*(3.0*hy-dh)) * hy36 * _delta;
                            }
                        }
                    }
                }

                //*****************************************************************************************************
            }

            DoubleVector x1(M+1);
            LinearEquation::GaussianElimination(w1,d1,x1);
            //IPrinter::printVector(x);

            w1.clear();
            d1.clear();

            for (unsigned int m=0; m<=M; m++) uh[m][n] = x1[m];
            x1.clear();
        }

        IPrinter::printMatrix(18, 14, uh);
        IPrinter::printSeperatorLine();

        tn.i = l;
        tn.t = l*ht;

        for (unsigned int m=0; m<=M; m++)
        {
            sn.j = m; sn.y = m*hy;

            DoubleMatrix w2(N+1, N+1, 0.0);
            DoubleVector d2(N+1, 0.0);

            for (unsigned int n=0; n<=N; n++)
            {
                sn.i = n; sn.x = n*hx;

                d2[n] = 2.0*uh[m][n] + lambda0*theta*ht + ht*f(sn, tn);

                if (m==0)       d2[n] += ((a*a*ht)/(hy*hy))*(uh[0][n]   - 2.0*uh[1][n]   + uh[2][n]);
                if (m>0 && m<M) d2[n] += ((a*a*ht)/(hy*hy))*(uh[m-1][n] - 2.0*uh[m][n]   + uh[m+1][n]);
                if (m==M)       d2[n] += ((a*a*ht)/(hy*hy))*(uh[M-2][n] - 2.0*uh[M-1][n] + uh[M][n]);

                if (n == 0)
                {
                    w2[0][0] = 2.0 + (2.0*a*a*ht)/(hx*hx) + lambda0*ht + (2.0*a*a*lambda*ht)/(hx);
                    w2[0][1] = -(2.0*a*a*ht)/(hx*hx);

                    d2[0] += (2.0*a*a*lambda*theta*ht)/(hx) + ((2.0*a*a*ht)/(hx))*g1(sn, tn);
                }
                else if (n == N)
                {
                    w2[N][N-1] = -(2.0*a*a*ht)/(hx*hx);
                    w2[N][N-0] = 2.0 + (2.0*a*a*ht)/(hx*hx) + lambda0*ht + (2.0*a*a*lambda*ht)/(hx);

                    d2[N] += (2.0*a*a*lambda*theta*ht)/(hx) + ((2.0*a*a*ht)/(hx))*g2(sn, tn);
                }
                else
                {
                    w2[n][n-1] = -(a*a*ht)/(hx*hx);
                    w2[n][n+0] = 2.0 + (2.0*a*a*ht)/(hx*hx) + lambda0*ht;
                    w2[n][n+1] = -(a*a*ht)/(hx*hx);
                }

                // Adding delta part ***********************************************************************************

                for (unsigned int i=0; i<Lc; i++)
                {
                    double _delta = delta(sn, i);

                    for (unsigned int j=0; j<Lo; j++)
                    {
                        d2[n] += -ht*k[i][j]*z[i][j] * _delta;

                        unsigned int jinx = (unsigned int)(xi[j].x*N);
                        //unsigned int jiny = (unsigned int)(xi[j].y*M);

                        //w2[n][jinx] += -ht*k[i][j]*_delta;

                        double hx3 = hx*hx*hx;
                        double hx32 = (1.0/(2.0*hx3));
                        double hx36 = (1.0/(6.0*hx3));

                        double hy3 = hy*hy*hy;
                        double hy32 = (1.0/(2.0*hy3));
                        double hy36 = (1.0/(6.0*hy3));

                        for (unsigned int n1=0; n1<=N; n1++)
                        {

                            double dh = fabs(n1*hx - xi[j].x);

                            if (dh <= hx)
                            {
                                w2[n][n1] += -ht*k[i][j] * ((2.0*hx-dh)*(hx-dh)*(hx+dh)) * hx32 * _delta;
                            }

                            if (hx < dh && dh <= 2.0*hx)
                            {
                                w2[n][n1] += -ht*k[i][j] * ((2.0*hx-dh)*(hx-dh)*(3.0*hx-dh)) * hx36 * _delta;
                            }
                        }
                    }
                }

                //*****************************************************************************************************
            }

            DoubleVector x2(N+1);
            LinearEquation::GaussianElimination(w2,d2,x2);
            //IPrinter::printVector(x);

            w2.clear();
            d2.clear();

            for (unsigned int n=0; n<=N; n++) u[m][n] = x2[n];
            x2.clear();
        }

        IPrinter::printMatrix(18, 14, u);
        IPrinter::printSeperatorLine();
    }
}

double IProblem2Forward2D::initial(const SpaceNodePDE &sn) const
{
    double x = sn.x;
    double y = sn.y;
    return x*x + y*y;
}

double IProblem2Forward2D::boundary(const SpaceNodePDE &sn UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM, BoundaryType boundary UNUSED_PARAM) const
{
    return NAN;
}

double IProblem2Forward2D::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    double x = sn.x;
    double y = sn.y;
    double t = tn.t;

    double res = 1.0 - 4.0*a*a + lambda0*(x*x + y*y + t - theta);

    double W = 0.0;
    for (unsigned int i=0; i<Lc; i++)
    {
        double vi = 0.0;
        for (unsigned int j=0; j<Lo; j++)
        {
            SpaceNodePDE xij = xi.at(j);
            vi += k[i][j] * ( U(xij.x, xij.y, t) - z[i][j] );
        }
        W += vi*delta(sn, i);
    }
    res -= W;

    return res;
}

double IProblem2Forward2D::delta(const SpaceNodePDE &sn UNUSED_PARAM, unsigned int i UNUSED_PARAM) const
{
    Dimension xd = spaceDimension(Dimension::DimensionX);
    Dimension yd = spaceDimension(Dimension::DimensionY);
    double hx = xd.step();
    double hy = yd.step();

    // Approximation delta function using normal distribution formula
    static double sigmaX = 3.0*hx;
    static double sigmaY = 3.0*hy;
    //double res = (1.0/(2.0*M_PI*sigmaX*sigmaY)) * exp( - ((sn.x-eta[i].x)*(sn.x-eta[i].x))/(2.0*sigmaX*sigmaX) - ((sn.y-eta[i].y)*(sn.y-eta[i].y))/(2.0*sigmaY*sigmaY) );
    double res = (1.0/(2.0*M_PI*sigmaX*sigmaY)) *
            exp(-0.5*(((sn.x-eta[i].x)*(sn.x-eta[i].x))/(sigmaX*sigmaX)+((sn.y-eta[i].y)*(sn.y-eta[i].y))/(sigmaY*sigmaY)));
    return res;

    // Approximation delta function using L4 Lagrange interpolation
//    double h3 = hx*hx*hx;
//    double h32 = (1.0/(2.0*h3));
//    double h36 = (1.0/(6.0*h3));
//    double dh = fabs(sn.x - eta[i]);
//    if (dh <= hx)                return ((2.0*hx-dh)*(hx-dh)*(hx+dh)) * h32;
//    if (hx < dh && dh <= 2.0*hx) return ((2.0*hx-dh)*(hx-dh)*(3.0*hx-dh)) * h36;
//    return 0.0;

    //double res = 0.0;
    //if ( sn.x == eta[i].x && sn.y == eta[i].y ) res = 1.0/(hx*hy);
    //return res;
}

double IProblem2Forward2D::g1(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    double y = sn.y;
    double t = tn.t;
    return lambda*(y*y + t - theta);
}

double IProblem2Forward2D::g2(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    double y = sn.y;
    double t = tn.t;
    return 2.0 + lambda*(1.0 + y*y + t - theta);
}

double IProblem2Forward2D::g3(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    double x = sn.x;
    double t = tn.t;
    return lambda*(x*x + t - theta);
}

double IProblem2Forward2D::g4(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    double x = sn.x;
    double t = tn.t;
    return 2.0 + lambda*(1.0 + x*x + t - theta);
}

double IProblem2Forward2D::U(double x, double y, double t) const
{
    return x*x + y*y + t;
}

