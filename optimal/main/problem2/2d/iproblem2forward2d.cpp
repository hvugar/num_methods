#include "iproblem2forward2d.h"

IProblem2Forward2D::IProblem2Forward2D(double a, double lambda0, double lambda, double theta, unsigned int Lc, unsigned int Lo)
{
    this->a = a;
    this->lambda0 = lambda0;
    this->lambda = lambda;
    this->theta = theta;

    this->Lc = Lc;
    this->Lo = Lo;

    k.resize(Lc, Lo);
    z.resize(Lc, Lo);

    for (unsigned int i=0; i<Lc; i++)
    {
        for (unsigned int j=0; j<Lo; j++)
        {
            k[i][j] = 5.0;
            z[i][j] = 8.0;
        }
    }

    Lo = 5;
    xi.resize(Lo);
    xi[0].x = 0.25; xi[0].i = 25; xi[0].y = 0.25; xi[0].j = 25;
    xi[1].x = 0.25; xi[1].i = 25; xi[1].y = 0.75; xi[1].j = 75;
    xi[2].x = 0.75; xi[2].i = 75; xi[2].y = 0.75; xi[2].j = 75;
    xi[3].x = 0.75; xi[3].i = 75; xi[3].y = 0.25; xi[3].j = 25;
    xi[4].x = 0.50; xi[4].i = 50; xi[4].y = 0.50; xi[4].j = 50;

    Lc = 4;
    eta.resize(Lc);
    eta[0].x = 0.33; eta[0].i = 33; eta[0].y = 0.33; eta[0].y = 33;
    eta[1].x = 0.33; eta[1].i = 33; eta[1].y = 0.66; eta[1].y = 66;
    eta[2].x = 0.66; eta[2].i = 66; eta[2].y = 0.66; eta[2].y = 66;
    eta[3].x = 0.66; eta[3].i = 66; eta[3].y = 0.33; eta[3].y = 33;
}

void IProblem2Forward2D::calculateMVD(DoubleMatrix &u)
{
    Dimension xd = spaceDimension(Dimension::DimensionX);
    Dimension yd = spaceDimension(Dimension::DimensionY);
    Dimension td = timeDimension();

    unsigned int N = xd.sizeN();
    unsigned int M = yd.sizeN();
    unsigned int K = td.sizeN();
    double hx = xd.step();
    double hy = yd.step();
    double ht = td.step();

    u.clear();
    u.resize(N+1, M+1);

    DoubleMatrix uh(N+1, M+1);

    SpaceNodePDE sn;
    for (unsigned int m=0; m<=N; m++)
    {
        sn.j = m; sn.y = m*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = n; sn.x = i*hx;
            u[m][n] = initial(sn);
        }
    }

    TimeNodePDE tn;
    for (unsigned int k=1; k<=M; k++)
    {
        tn.i = k;
        tn.t = k*ht - 0.5*ht;

        DoubleMatrix w(N+1, N+1, 0.0);
        DoubleVector d(N+1, 0.0);

        for (unsigned int m=0; m<=M; m++)
        {
            sn.i = m; sn.x = m*hx;

            d[m] = u[k-1][m] + lambda0*ht*theta + ht*f(sn, tn);

            if (m == 0)
            {
                w[0][0] = 1.0 + (2.0*a*a*ht)/(hx*hx) + lambda0*ht + (2.0*a*a*ht*lambda1)/(hx);
                w[0][1] = -(2.0*a*a*ht)/(hx*hx);

                d[0] += (2.0*a*a*ht*lambda1*theta)/(hx) - ((2.0*a*a*ht)/(hx))*g0(tn);
            }
            else if (m == N)
            {
                w[N][N-1] = -(2.0*a*a*ht)/(hx*hx);
                w[N][N-0] = 1.0 + (2.0*a*a*ht)/(hx*hx) + lambda0*ht + (2.0*a*a*ht*lambda2)/(hx);

                d[N] += (2.0*a*a*ht*lambda2*theta)/(hx) + ((2.0*a*a*ht)/(hx))*g1(tn);
            }
            else
            {
                w[m][m-1] = -(a*a*ht)/(hx*hx);
                w[m][m+0] = 1.0 + (2.0*a*a*ht)/(hx*hx) + lambda0*ht;
                w[m][m+1] = -(a*a*ht)/(hx*hx);
            }


            for (unsigned int i=0; i<Lc; i++)
            {
                double _delta = delta(sn, i);
                //if ( dt > 0.0 )
                //if (n*hx == eta[i])
                {
                    for (unsigned int j=0; j<Lo; j++)
                    {
                        d[m] += -ht*k[i][j]*z[i][j] * _delta;

                        unsigned int jinx = (unsigned int)(xi[j]*N);
                        //w[n][jinx] += -ht*k[i][j];

                        double h3 = hx*hx*hx;
                        double h32 = (1.0/(2.0*h3));
                        double h36 = (1.0/(6.0*h3));

                        for (unsigned int n1=jinx-2; n1<=jinx+2; n1++)
                        {

                            double dh = fabs(n1*hx - xi[j]);

                            if (dh <= hx)
                            {
                                w[m][n1] += -ht*k[i][j] * ((2.0*hx-dh)*(hx-dh)*(hx+dh)) * h32 * _delta;
                            }

                            if (hx < dh && dh <= 2.0*hx)
                            {
                                w[m][n1] += -ht*k[i][j] * ((2.0*hx-dh)*(hx-dh)*(3.0*hx-dh)) * h36 * _delta;
                            }
                        }
                    }

                }
            }
        }

        DoubleVector x(N+1);
        LinearEquation::GaussianElimination(w,d,x);
        //IPrinter::printVector(x);

        w.clear();
        d.clear();

        for (unsigned int n=0; n<=N; n++) u[k][n] = x[n];


        //        FILE* file = fopen("D:\\data_m_1.txt", "w");
        //        IPrinter::printMatrix(10, 4, w, M, N, NULL, file);
        //        fclose(file);

        //        FILE* file1 = fopen("D:\\data_d_1.txt", "w");
        //        IPrinter::printVector(d, NULL, N, 0,0, file1);
        //        fclose(file1);
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
            vi += k.at(i,j) * ( U(xij.x, xij.y, t) - z.at(i,j) );
        }
        W += vi*delta(sn, i);
    }
    res -= W;

    return res;
}

double IProblem2Forward2D::delta(const SpaceNodePDE &sn UNUSED_PARAM, unsigned int i UNUSED_PARAM) const
{
    return 1.0;
}

double IProblem2Forward2D::g1(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    double y = sn.y;
    double t = tn.t;
    return -lambda*(y*y + t - theta);
}

double IProblem2Forward2D::g2(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    double y = sn.y;
    double t = tn.t;
    return 2.0 - lambda*(1.0 + y*y + t - theta);
}

double IProblem2Forward2D::g3(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    double x = sn.x;
    double t = tn.t;
    return -lambda*(x*x + t - theta);
}

double IProblem2Forward2D::g4(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    double x = sn.x;
    double t = tn.t;
    return 2.0 - lambda*(1.0 + x*x + t - theta);
}

double IProblem2Forward2D::U(double x, double y, double t) const
{
    return x*x + y*y + t;
}

