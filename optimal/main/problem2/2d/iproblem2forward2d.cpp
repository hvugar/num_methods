#include "iproblem2forward2d.h"

IProblem2Forward2D::IProblem2Forward2D() {}

void IProblem2Forward2D::setSettings(double a, double lambda0, double lambda, double theta, unsigned int Lc, unsigned int Lo)
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
    eta[0].x = 0.30; eta[0].y = 0.20; eta[0].i = 30; eta[0].j = 20;
    eta[1].x = 0.20; eta[1].y = 0.80; eta[1].i = 20; eta[1].j = 80;

    xi.resize(this->Lo);
    xi[0].x = 0.10;  xi[0].y = 0.80; xi[0].i = 10; xi[0].j = 80;
    xi[1].x = 0.80;  xi[1].y = 0.50; xi[1].i = 80; xi[1].j = 50;
    xi[2].x = 0.40;  xi[2].y = 0.70; xi[2].i = 40; xi[2].j = 70;

    //xi.resize(this->Lo);
    //xi[0].x = 0.25;  xi[0].y = 0.25; //xi[0].i = 25; xi[0].j = 25;
    //xi[1].x = 0.25;  xi[1].y = 0.75; //xi[1].i = 25; xi[1].j = 75;
    //xi[2].x = 0.75;  xi[2].y = 0.75; //xi[2].i = 75; xi[2].j = 75;
    //xi[3].x = 0.75;  xi[3].y = 0.25; //xi[3].i = 75; xi[3].j = 25;
    //xi[4].x = 0.50;  xi[4].y = 0.50; //xi[4].i = 50; xi[4].j = 50;

    //eta.resize(this->Lc);
    //eta[0].x = 0.33; eta[0].y = 0.33; //eta[0].i = 33; eta[0].j = 33;
    //eta[1].x = 0.33; eta[1].y = 0.66; //eta[1].i = 33; eta[1].j = 66;
    //eta[2].x = 0.66; eta[2].y = 0.66; //eta[2].i = 66; eta[2].j = 66;
    //eta[3].x = 0.66; eta[3].y = 0.33; //eta[3].i = 66; eta[3].j = 33;
}

void IProblem2Forward2D::calculateMVD(DoubleMatrix &u) const
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

    //************************************* initial conditions *************************************//
    for (unsigned int m=0; m<=M; m++)
    {
        SpaceNodePDE sn;
        sn.j = m; sn.y = m*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = n; sn.x = n*hx;
            u[m][n] = initial(sn);
        }
    }
    //************************************* initial conditions *************************************//

    //IPrinter::printMatrix(10, 6, u);
    //IPrinter::printSeperatorLine();

    TimeNodePDE tn;
    for (unsigned int l=1; l<=M; l++)
    {
        //************************************* approximatin to y direction conditions *************************************//
        tn.i = l;
        tn.t = l*ht - 0.5*ht;

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

                d1[offset+m] = 2.0*u[m][n] + lambda0*theta*ht + ht*f(sn1, tn);

                if (n==0)       d1[offset+m] += ((a*a*ht)/(hx*hx))*(u[m][0]   - 2.0*u[m][1]   + u[m][2]);
                if (n>0 && n<N) d1[offset+m] += ((a*a*ht)/(hx*hx))*(u[m][n-1] - 2.0*u[m][n]   + u[m][n+1]);
                if (n==N)       d1[offset+m] += ((a*a*ht)/(hx*hx))*(u[m][N-2] - 2.0*u[m][N-1] + u[m][N]);

                if (m == 0)
                {
                    w1[offset+0][offset+0] += 2.0 + (2.0*a*a*ht)/(hy*hy) + lambda0*ht + (2.0*a*a*lambda*ht)/(hy);
                    w1[offset+0][offset+1] += -(2.0*a*a*ht)/(hy*hy);

                    d1[offset+0] += (2.0*a*a*lambda*theta*ht)/(hy) + ((2.0*a*a*ht)/(hy))*g3(sn1, tn);
                }
                else if (m == M)
                {
                    w1[offset+M][offset+M-1] += -(2.0*a*a*ht)/(hy*hy);
                    w1[offset+M][offset+M-0] += 2.0 + (2.0*a*a*ht)/(hy*hy) + lambda0*ht + (2.0*a*a*lambda*ht)/(hy);

                    d1[offset+M] += (2.0*a*a*lambda*theta*ht)/(hy) + ((2.0*a*a*ht)/(hy))*g4(sn1, tn);
                }
                else
                {
                    w1[offset+m][offset+m-1] += -(a*a*ht)/(hy*hy);
                    w1[offset+m][offset+m+0] += 2.0 + (2.0*a*a*ht)/(hy*hy) + lambda0*ht;
                    w1[offset+m][offset+m+1] += -(a*a*ht)/(hy*hy);
                }

                //************************************* Adding delta part *************************************//
                for (unsigned int i=0; i<Lc; i++)
                {
                    double _delta = delta(sn1, i, 2);

                    for (unsigned int j=0; j<Lo; j++)
                    {
                        d1[offset+m] += -ht*k[i][j]*z[i][j] * _delta;

                        unsigned int jinx = xi[j].i;
                        unsigned int jiny = xi[j].j;
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
                uh[m][n] = x1[n*(M+1)+m];
            }
        }

        x1.clear();

        //IPrinter::printMatrix(10, 6, uh);
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

                d2[offset+n] = 2.0*uh[m][n] + lambda0*theta*ht + ht*f(sn2, tn);

                if (m==0)       d2[offset+n] += ((a*a*ht)/(hy*hy))*(uh[0][n]   - 2.0*uh[1][n]   + uh[2][n]);
                if (m>0 && m<M) d2[offset+n] += ((a*a*ht)/(hy*hy))*(uh[m-1][n] - 2.0*uh[m][n]   + uh[m+1][n]);
                if (m==M)       d2[offset+n] += ((a*a*ht)/(hy*hy))*(uh[M-2][n] - 2.0*uh[M-1][n] + uh[M][n]);

                if (n == 0)
                {
                    w2[offset+0][offset+0] += 2.0 + (2.0*a*a*ht)/(hx*hx) + lambda0*ht + (2.0*a*a*lambda*ht)/(hx);
                    w2[offset+0][offset+1] += -(2.0*a*a*ht)/(hx*hx);

                    d2[offset+0] += (2.0*a*a*lambda*theta*ht)/(hx) + ((2.0*a*a*ht)/hx)*g1(sn2, tn);
                }
                else if (n == N)
                {
                    w2[offset+N][offset+N-1] += -(2.0*a*a*ht)/(hx*hx);
                    w2[offset+N][offset+N-0] += 2.0 + (2.0*a*a*ht)/(hx*hx) + lambda0*ht + (2.0*a*a*lambda*ht)/(hx);

                    d2[offset+N] += (2.0*a*a*lambda*theta*ht)/(hx) + ((2.0*a*a*ht)/(hx))*g2(sn2, tn);
                }
                else
                {
                    w2[offset+n][offset+n-1] += -(a*a*ht)/(hx*hx);
                    w2[offset+n][offset+n+0] += 2.0 + (2.0*a*a*ht)/(hx*hx) + lambda0*ht;
                    w2[offset+n][offset+n+1] += -(a*a*ht)/(hx*hx);
                }

                //************************************* Adding delta part *************************************//
                for (unsigned int i=0; i<Lc; i++)
                {
                    double _delta = delta(sn2, i, 2);

                    for (unsigned int j=0; j<Lo; j++)
                    {
                        d2[offset+n] += -ht*k[i][j]*z[i][j] * _delta;

                        unsigned int jinx = xi[j].i;
                        unsigned int jiny = xi[j].j;
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
                u[m][n] = x2[m*(N+1)+n];
            }
        }

        x2.clear();

        //IPrinter::printMatrix(10, 6, u);
        //IPrinter::printSeperatorLine();
        //************************************* approximatin to x direction conditions *************************************//
    }
}

void IProblem2Forward2D::calculateMVD1(DoubleMatrix &u) const
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

    //************************************* initial conditions *************************************//
    for (unsigned int m=0; m<=M; m++)
    {
        SpaceNodePDE sn;
        sn.j = m; sn.y = m*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = n; sn.x = n*hx;
            u[m][n] = initial(sn);
        }
    }
    //************************************* initial conditions *************************************//

    IPrinter::printMatrix(10, 6, u);
    IPrinter::printSeperatorLine();

    TimeNodePDE tn;
    for (unsigned int l=1; l<=1; l++)
    {
        //************************************* approximatin to y direction conditions *************************************//
        tn.i = l;
        tn.t = l*ht - 0.5*ht;

        for (unsigned int n=0; n<=N; n++)
        {
            SpaceNodePDE sn1;
            sn1.i = n; sn1.x = n*hx;

            DoubleVector a1(M+1, 0.0);
            DoubleVector b1(M+1, 0.0);
            DoubleVector c1(M+1, 0.0);
            DoubleVector d1(M+1, 0.0);
            DoubleVector x1(M+1, 0.0);

            for (unsigned int m=0; m<=M; m++)
            {
                sn1.j = m; sn1.y = m*hy;

                d1[m] = 2.0*u[m][n] + lambda0*theta*ht + ht*f(sn1, tn);

                if (n==0)       d1[m] += ((a*a*ht)/(hx*hx))*(u[m][0]   - 2.0*u[m][1]   + u[m][2]);
                if (n>0 && n<N) d1[m] += ((a*a*ht)/(hx*hx))*(u[m][n-1] - 2.0*u[m][n]   + u[m][n+1]);
                if (n==N)       d1[m] += ((a*a*ht)/(hx*hx))*(u[m][N-2] - 2.0*u[m][N-1] + u[m][N]);

                if (m == 0)
                {
                    d1[0] += (2.0*a*a*lambda*theta*ht)/(hy) + ((2.0*a*a*ht)/(hy))*g3(sn1, tn);

                    b1[0] = 2.0 + (2.0*a*a*ht)/(hy*hy) + lambda0*ht + (2.0*a*a*lambda*ht)/(hy);
                    c1[0] = -(2.0*a*a*ht)/(hy*hy);

                    //w1[0][0] += 2.0 + (2.0*a*a*ht)/(hy*hy) + lambda0*ht + (2.0*a*a*lambda*ht)/(hy);
                    //w1[0][ 1] += -(2.0*a*a*ht)/(hy*hy);

                }
                else if (m == M)
                {
                    d1[M] += (2.0*a*a*lambda*theta*ht)/(hy) + ((2.0*a*a*ht)/(hy))*g4(sn1, tn);

                    a1[M] = -(2.0*a*a*ht)/(hy*hy);
                    b1[M] = 2.0 + (2.0*a*a*ht)/(hy*hy) + lambda0*ht + (2.0*a*a*lambda*ht)/(hy);

                    //w1[offset+M][offset+M-1] += -(2.0*a*a*ht)/(hy*hy);
                    //w1[offset+M][offset+M-0] += 2.0 + (2.0*a*a*ht)/(hy*hy) + lambda0*ht + (2.0*a*a*lambda*ht)/(hy);

                }
                else
                {
                    a1[m] = -(a*a*ht)/(hy*hy);
                    b1[m] = 2.0 + (2.0*a*a*ht)/(hy*hy) + lambda0*ht;
                    c1[m] = -(a*a*ht)/(hy*hy);

                    //w1[offset+m][offset+m-1] += -(a*a*ht)/(hy*hy);
                    //w1[offset+m][offset+m+0] += 2.0 + (2.0*a*a*ht)/(hy*hy) + lambda0*ht;
                    //w1[offset+m][offset+m+1] += -(a*a*ht)/(hy*hy);
                }

                //************************************* Adding delta part *************************************//
                //for (unsigned int i=0; i<Lc; i++)
                //{
                //    double _delta = delta(sn1, i, 2);

                //    for (unsigned int j=0; j<Lo; j++)
                //    {
                //        d1[offset+m] += -ht*k[i][j]*z[i][j] * _delta;

                //        unsigned int jinx = xi[j].i;
                //        unsigned int jiny = xi[j].j;
                //        w1[offset+m][jinx*(M+1)+jiny] += -ht*k[i][j] * _delta;
                //    }
                //}
                //************************************* Adding delta part *************************************//
            }

            tomasAlgorithm(a1.data(), b1.data(), c1.data(), d1.data(), x1.data(), M+1);

            for (unsigned int m=0; m<=M; m++) uh[m][n] = x1[m];

            a1.clear();
            b1.clear();
            c1.clear();
            d1.clear();
            x1.clear();

        }

        //DoubleVector x1((M+1)*(N+1));
        //LinearEquation::GaussianElimination(w1,d1,x1);

        //w1.clear();
        //d1.clear();

        for (unsigned int m=0; m<=M; m++)
        {
            for (unsigned int n=0; n<=N; n++)
            {
                uh[m][n] = U(n*hx, m*hy, tn.t);/*x1[n*(M+1)+m]*/;
            }
        }

        //x1.clear();

        IPrinter::printMatrix(10, 6, uh);
        IPrinter::printSeperatorLine();
        //************************************* approximatin to y direction conditions *************************************//

        //************************************* approximatin to x direction conditions *************************************//
        tn.i = l;
        tn.t = l*ht;

        //DoubleMatrix w2((N+1)*(M+1), (N+1)*(M+1), 0.0);
        //DoubleVector d2((N+1)*(M+1), 0.0);

        for (unsigned int m=0; m<=M; m++)
        {
            SpaceNodePDE sn2;
            sn2.j = m; sn2.y = m*hy;

            DoubleVector a2(N+1, 0.0);
            DoubleVector b2(N+1, 0.0);
            DoubleVector c2(N+1, 0.0);
            DoubleVector d2(N+1, 0.0);
            DoubleVector x2(N+1, 0.0);

            bool ok = false;
            for (unsigned int i=0; i<Lc; i++)
            {
                if (sn2.j == eta[i].j) ok = true;
            }

            if (ok == false)
            {
                for (unsigned int n=0; n<=N; n++)
                {
                    sn2.i = n; sn2.x = n*hx;

                    d2[n] = 2.0*uh[m][n] + lambda0*theta*ht + ht*f(sn2, tn);

                    if (m==0)       d2[n] += ((a*a*ht)/(hy*hy))*(uh[0][n]   - 2.0*uh[1][n]   + uh[2][n]);
                    if (m>0 && m<M) d2[n] += ((a*a*ht)/(hy*hy))*(uh[m-1][n] - 2.0*uh[m][n]   + uh[m+1][n]);
                    if (m==M)       d2[n] += ((a*a*ht)/(hy*hy))*(uh[M-2][n] - 2.0*uh[M-1][n] + uh[M][n]);

                    if (n == 0)
                    {
                        b2[0] += 2.0 + (2.0*a*a*ht)/(hx*hx) + lambda0*ht + (2.0*a*a*lambda*ht)/(hx);
                        c2[0] += -(2.0*a*a*ht)/(hx*hx);

                        d2[0] += (2.0*a*a*lambda*theta*ht)/(hx) + ((2.0*a*a*ht)/hx)*g1(sn2, tn);
                    }
                    else if (n == N)
                    {
                        a2[N] += -(2.0*a*a*ht)/(hx*hx);
                        b2[N] += 2.0 + (2.0*a*a*ht)/(hx*hx) + lambda0*ht + (2.0*a*a*lambda*ht)/(hx);

                        d2[N] += (2.0*a*a*lambda*theta*ht)/(hx) + ((2.0*a*a*ht)/(hx))*g2(sn2, tn);
                    }
                    else
                    {
                        a2[n] += -(a*a*ht)/(hx*hx);
                        b2[n] += 2.0 + (2.0*a*a*ht)/(hx*hx) + lambda0*ht;
                        c2[n] += -(a*a*ht)/(hx*hx);
                    }

                    //************************************* Adding delta part *************************************//
                    for (unsigned int i=0; i<Lc; i++)
                    {
                        double _delta = delta(sn2, i, 2);

                        for (unsigned int j=0; j<Lo; j++)
                        {
                            d2[n] += -ht*k[i][j]*z[i][j] * _delta;

                            //unsigned int jinx = xi[j].i;
                            //unsigned int jiny = xi[j].j;
                            //w2[offset+n][jiny*(N+1)+jinx] += -ht*k[i][j]*_delta;
                        }
                    }
                    //************************************* Adding delta part *************************************//
                }

                tomasAlgorithm(a2.data(), b2.data(), c2.data(), d2.data(), x2.data(), N+1);

                for (unsigned int n=0; n<=N; n++) u[m][n] = x2[m];
            }

            a2.clear();
            b2.clear();
            c2.clear();
            d2.clear();
            x2.clear();
        }

        //DoubleVector x2((N+1)*(M+1));
        //LinearEquation::GaussianElimination(w2,d2,x2);

        //w2.clear();
        //d2.clear();

        //        for (unsigned int m=0; m<=M; m++)
        //        {
        //            for (unsigned int n=0; n<=N; n++)
        //            {
        //                u[m][n] = x2[m*(N+1)+n];
        //            }
        //        }

        //        x2.clear();

        IPrinter::printMatrix(10, 6, u);
        IPrinter::printSeperatorLine();
        //************************************* approximatin to x direction conditions *************************************//
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

    //printf("t: %f\n", t);

    double res = 1.0 - 4.0*a*a + lambda0*(U(x,y,t) - theta);

    double W = 0.0;
    for (unsigned int i=0; i<Lc; i++)
    {
        double _delta = delta(sn, i, 0);
        double vi = 0.0;
        for (unsigned int j=0; j<Lo; j++)
        {
            vi += k[i][j] * ( U(xi[j].x, xi[j].y, t)  - z[i][j]);
        }
        W += vi * _delta;
    }
    res -= W;

    return res;
}

double IProblem2Forward2D::delta(const SpaceNodePDE &sn UNUSED_PARAM, unsigned int i UNUSED_PARAM, unsigned int source UNUSED_PARAM) const
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
    //        exp(-0.5*(((sn.x-eta[i].x)*(sn.x-eta[i].x))/(sigmaX*sigmaX)+((sn.y-eta[i].y)*(sn.y-eta[i].y))/(sigmaY*sigmaY)));

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
    if ( sn.i == eta[i].i && sn.j == eta[i].j ) res = 1.0/(hx*hy);

    return res;
}

double IProblem2Forward2D::g1(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    double y = sn.y;
    double t = tn.t;
    return lambda*(y*y+t - theta);
}

double IProblem2Forward2D::g2(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    double y = sn.y;
    double t = tn.t;
    //return 2.0 + lambda*(U(1.0, y, t) - theta);
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

