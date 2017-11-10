#include "iproblem2backward2d.h"

IProblem2Backward2D::IProblem2Backward2D() {}

void IProblem2Backward2D::setSettings(double a, double lambda0, double lambda, double theta, unsigned int Lc, unsigned int Lo)
{
    this->a = a;
    this->lambda0 = lambda0;
    this->lambda = lambda;
    this->theta = theta;
    this->Lc = Lc;
    this->Lo = Lo;

    k.resize(this->Lc, this->Lo);
    z.resize(this->Lc, this->Lo);

    for (unsigned int i=0; i<this->Lc; i++)
    {
        for (unsigned int j=0; j<this->Lo; j++)
        {
            k[i][j] = 5.0;
            z[i][j] = 8.0;
        }
    }

    eta.resize(this->Lc);
    eta[0].x = 0.30; eta[0].y = 0.20; eta[0].i = 3; eta[0].j = 2;
    eta[1].x = 0.20; eta[1].y = 0.80; eta[1].i = 2; eta[1].j = 8;

    xi.resize(this->Lo);
    xi[0].x = 0.10;  xi[0].y = 0.80; xi[0].i = 1; xi[0].j = 8;
    xi[1].x = 0.80;  xi[1].y = 0.50; xi[1].i = 8; xi[1].j = 5;
    xi[2].x = 0.40;  xi[2].y = 0.70; xi[2].i = 4; xi[2].j = 7;
}

void IProblem2Backward2D::calculateMVD(DoubleMatrix &p)
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

    mu.clear();
    mu.resize(M+1, N+1, 1.0);

    p.clear();
    p.resize(M+1, N+1);

    DoubleMatrix ph(M+1, N+1);

    //************************************* initial conditions *************************************//
    for (unsigned int m=0; m<=M; m++)
    {
        SpaceNodePDE sn;
        sn.j = m; sn.y = m*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = n; sn.x = n*hx;
            p[m][n] = initial(sn);
        }
    }
    //************************************* initial conditions *************************************//

    //IPrinter::printMatrix(10, 6, p);
    //IPrinter::printSeperatorLine();

    TimeNodePDE tn;
    for (unsigned int l=M-1; l!=UINT32_MAX; l--)
    {
        //************************************* approximatin to y direction conditions *************************************//
        tn.i = l;
        tn.t = l*ht + 0.5*ht;

        DoubleMatrix w1((M+1)*(N+1), (M+1)*(N+1), 0.0);
        DoubleVector d1((M+1)*(N+1), 0.0);

        for (unsigned int n=0; n<=N; n++)
        {
            SpaceNodePDE sn1;
            sn1.i = n; sn1.x = n*hx;

            for (unsigned int m=0; m<=M; m++)
            {
                sn1.j = m; sn1.y = m*hy;
                unsigned int offset = n*(M+1);

                d1[offset+m] = 2.0*p[m][n] - ht*f(sn1, tn);

                if (n==0)       d1[offset+m] += ((a*a*ht)/(hx*hx))*(p[m][0]   - 2.0*p[m][1]   + p[m][2]);
                if (n>0 && n<N) d1[offset+m] += ((a*a*ht)/(hx*hx))*(p[m][n-1] - 2.0*p[m][n]   + p[m][n+1]);
                if (n==N)       d1[offset+m] += ((a*a*ht)/(hx*hx))*(p[m][N-2] - 2.0*p[m][N-1] + p[m][N]);

                if (m == 0)
                {
                    w1[offset+0][offset+0] += 2.0 + (2.0*a*a*ht)/(hy*hy) + lambda0*ht - (2.0*a*a*lambda*ht)/(hy);
                    w1[offset+0][offset+1] += -(2.0*a*a*ht)/(hy*hy);

                    d1[offset+0] += ((2.0*a*a*ht)/(hy))*g3(sn1, tn);
                }
                else if (m == M)
                {
                    w1[offset+M][offset+M-1] += -(2.0*a*a*ht)/(hy*hy);
                    w1[offset+M][offset+M-0] += 2.0 + (2.0*a*a*ht)/(hy*hy) + lambda0*ht - (2.0*a*a*lambda*ht)/(hy);

                    d1[offset+M] += ((2.0*a*a*ht)/(hy))*g4(sn1, tn);
                }
                else
                {
                    w1[offset+m][offset+m-1] += -(a*a*ht)/(hy*hy);
                    w1[offset+m][offset+m+0] += 2.0 + (2.0*a*a*ht)/(hy*hy) + lambda0*ht;
                    w1[offset+m][offset+m+1] += -(a*a*ht)/(hy*hy);
                }

                //************************************* Adding delta part *************************************//
                for (unsigned int j=0; j<Lo; j++)
                {
                    double _delta = delta(sn1, j, 2);

                    for (unsigned int i=0; i<Lc; i++)
                    {
                        //d1[offset+m] += -ht*k[i][j]*z[i][j] * _delta;

                        unsigned int jinx = eta[i].i;
                        unsigned int jiny = eta[i].j;
                        w1[offset+m][jinx*(M+1)+jiny] += -ht*k[i][j] * _delta;
                    }
                }
                //************************************* Adding delta part *************************************//
            }
        }

        DoubleVector x1((M+1)*(N+1));
        LinearEquation::GaussianElimination(w1,d1,x1);

        w1.clear();
        d1.clear();

        for (unsigned int m=0; m<=M; m++)
        {
            for (unsigned int n=0; n<=N; n++)
            {
                ph[m][n] = x1[n*(M+1)+m];
            }
        }

        x1.clear();

        //IPrinter::printMatrix(10, 6, ph);
        //IPrinter::printSeperatorLine();
        //************************************* approximatin to y direction conditions *************************************//

        //************************************* approximatin to x direction conditions *************************************//
        tn.i = l;
        tn.t = l*ht;

        DoubleMatrix w2((N+1)*(M+1), (N+1)*(M+1), 0.0);
        DoubleVector d2((N+1)*(M+1), 0.0);

        for (unsigned int m=0; m<=M; m++)
        {
            SpaceNodePDE sn2;
            sn2.j = m; sn2.y = m*hy;

            for (unsigned int n=0; n<=N; n++)
            {
                sn2.i = n; sn2.x = n*hx;

                unsigned int offset = m*(N+1);

                d2[offset+n] = 2.0*ph[m][n] - ht*f(sn2, tn);

                if (m==0)       d2[offset+n] += ((a*a*ht)/(hy*hy))*(ph[0][n]   - 2.0*ph[1][n]   + ph[2][n]);
                if (m>0 && m<M) d2[offset+n] += ((a*a*ht)/(hy*hy))*(ph[m-1][n] - 2.0*ph[m][n]   + ph[m+1][n]);
                if (m==M)       d2[offset+n] += ((a*a*ht)/(hy*hy))*(ph[M-2][n] - 2.0*ph[M-1][n] + ph[M][n]);

                if (n == 0)
                {
                    w2[offset+0][offset+0] += 2.0 + (2.0*a*a*ht)/(hx*hx) + lambda0*ht - (2.0*a*a*lambda*ht)/(hx);
                    w2[offset+0][offset+1] += -(2.0*a*a*ht)/(hx*hx);

                    d2[offset+0] += ((2.0*a*a*ht)/hx)*g1(sn2, tn);
                }
                else if (n == N)
                {
                    w2[offset+N][offset+N-1] += -(2.0*a*a*ht)/(hx*hx);
                    w2[offset+N][offset+N-0] += 2.0 + (2.0*a*a*ht)/(hx*hx) + lambda0*ht - (2.0*a*a*lambda*ht)/(hx);

                    d2[offset+N] += ((2.0*a*a*ht)/(hx))*g2(sn2, tn);
                }
                else
                {
                    w2[offset+n][offset+n-1] += -(a*a*ht)/(hx*hx);
                    w2[offset+n][offset+n+0] += 2.0 + (2.0*a*a*ht)/(hx*hx) + lambda0*ht;
                    w2[offset+n][offset+n+1] += -(a*a*ht)/(hx*hx);
                }

                //************************************* Adding delta part *************************************//
                for (unsigned int j=0; j<Lo; j++)
                {
                    double _delta = delta(sn2, j, 2);

                    for (unsigned int i=0; i<Lc; i++)
                    {
                        //d2[offset+n] += -ht*k[i][j]*z[i][j] * _delta;

                        unsigned int jinx = eta[i].i;
                        unsigned int jiny = eta[i].j;
                        w2[offset+n][jiny*(N+1)+jinx] += -ht*k[i][j]*_delta;
                    }
                }
                //************************************* Adding delta part *************************************//
            }
        }

        DoubleVector x2((N+1)*(M+1));
        LinearEquation::GaussianElimination(w2,d2,x2);

        w2.clear();
        d2.clear();

        for (unsigned int m=0; m<=M; m++)
        {
            for (unsigned int n=0; n<=N; n++)
            {
                p[m][n] = x2[m*(N+1)+n];
            }
        }

        x2.clear();

        //IPrinter::printMatrix(10, 6, p);
        //IPrinter::printSeperatorLine();
        //************************************* approximatin to x direction conditions *************************************//
    }
}

double IProblem2Backward2D::initial(const SpaceNodePDE &sn) const
{
    unsigned int i = sn.i;
    unsigned int j = sn.j;
    return -2.0 * mu[j][i] * (uT[j][i]-U[j][i]) + h(sn);
}

double IProblem2Backward2D::boundary(const SpaceNodePDE &sn UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM, BoundaryType boundary UNUSED_PARAM) const
{
    return NAN;
}

double IProblem2Backward2D::f(const SpaceNodePDE &sn UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM) const
{
    double x = sn.x;
    double y = sn.y;
    double t = tn.t;

    double res = 1.0 + 4.0*a*a - lambda0*(x*x + y*y + t);

    double W = 0.0;
    for (unsigned int j=0; j<Lo; j++)
    {
        double _delta = delta(sn, j, 0);
        double vi = 0.0;
        for (unsigned int i=0; i<Lc; i++)
        {
            vi += k[i][j] * P(eta[i].x, eta[i].y, t);
        }
        W += vi * _delta;
    }
    res += W;

    return res;
}

double IProblem2Backward2D::delta(const SpaceNodePDE &sn UNUSED_PARAM, unsigned int j UNUSED_PARAM, unsigned int source UNUSED_PARAM) const
{
    double res = 0.0;

    Dimension xd = spaceDimension(Dimension::DimensionX);
    Dimension yd = spaceDimension(Dimension::DimensionY);
    double hx = xd.step();
    double hy = yd.step();

    // Approximation delta function using normal distribution formula
    //static double sigmaX = hx;
    //static double sigmaY = hy;
    //res = (1.0/(2.0*M_PI*sigmaX*sigmaY)) *
    //        exp(-0.5*(((sn.x-xi[j].x)*(sn.x-xi[j].x))/(sigmaX*sigmaX)+((sn.y-xi[j].y)*(sn.y-xi[j].y))/(sigmaY*sigmaY)));

    //if ( sn.x == eta[i].x && sn.y == eta[i].y ) res = 1.0/(hx*hy);
    if ( sn.i == xi[j].i && sn.j == xi[j].j ) res = 1.0/(hx*hy);

    return res;
}

double IProblem2Backward2D::g1(const SpaceNodePDE &sn UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM) const
{
    return -lambda*(sn.y*sn.y+tn.t);
}

double IProblem2Backward2D::g2(const SpaceNodePDE &sn UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM) const
{
    return 2.0 - lambda*(1.0+sn.y*sn.y+tn.t);
}

double IProblem2Backward2D::g3(const SpaceNodePDE &sn UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM) const
{
    return -lambda*(sn.x*sn.x+tn.t);
}

double IProblem2Backward2D::g4(const SpaceNodePDE &sn UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM) const
{
    return 2.0 - lambda*(sn.x*sn.x+1.0+tn.t);
}

double IProblem2Backward2D::h(const SpaceNodePDE &sn) const
{
    double x = sn.x; unsigned int i = sn.i;
    double y = sn.y; unsigned int j = sn.j;

    return (x*x + y*y + 1.0) + 2.0 * mu[j][i] * (uT[j][i] - U[j][i]);
}

double IProblem2Backward2D::P(double x, double y, double t) const
{
    return x*x + y*y + t;
}
