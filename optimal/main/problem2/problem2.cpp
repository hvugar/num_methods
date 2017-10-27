#include "problem2.h"

void IProblem2::Main(int argc, char *argv[])
{
    IProblem2 ip2;
    ip2.setTimeDimension(Dimension(0.1, 0, 10));
    ip2.addSpaceDimension(Dimension(0.01, 0, 100));

    DoubleMatrix u;
    ip2.gridMethod(u);

    //IPrinter::printMatrix(10, 6, u);
}

IProblem2::IProblem2()
{
    a = 1.0;
    lambda0 = 1.0;
    lambda1 = 2.0;
    lambda2 = 1.5;
    theta = 5.0;

    Lc = 2;
    Lo = 3;
    k.resize(Lc, Lo);
    z.resize(Lc, Lo);

    eta.resize(Lc);
    eta[0] = 0.33;
    eta[1] = 0.66;

    xi.resize(Lo);
    xi[0] = 0.2514;
    xi[1] = 0.5045;
    xi[2] = 0.7515;

    for (unsigned int i=0; i<Lc; i++)
    {
        for (unsigned int j=0; j<Lo; j++)
        {
            k[i][j] = 5.0;//(i+1)*(j+2);
            z[i][j] = 8.0;
        }
    }
}

void IProblem2::gridMethod(DoubleMatrix &u) const
{
    Dimension sd = spaceDimension(Dimension::DimensionX);
    Dimension td = timeDimension();

    unsigned int N = sd.sizeN();
    double hx = sd.step();
    unsigned int M = td.sizeN();
    double ht = td.step();

    u.clear();
    u.resize(M+1, N+1);

    SpaceNodePDE sn;
    for (unsigned int n=0; n<=N; n++)
    {
        sn.i = n; sn.x = n*hx;
        u[0][n] = initial(sn);
    }

    TimeNodePDE tn;
    for (unsigned int m=1; m<=1; m++)
    {
        tn.i = m;
        tn.t = m*ht;

        DoubleMatrix w(N+1, N+1, 0.0);
        DoubleVector d(N+1, 0.0);

        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = n; sn.x = n*hx;
            d[n] = u[m-1][n] + lambda0*ht*theta + ht*f(sn, tn);

            if (n == 0)
            {
                w[0][0] = 1.0 + (2.0*a*a*ht)/(hx*hx) + lambda0*ht + (2.0*a*a*ht*lambda1)/(hx);
                w[0][1] = -(2.0*a*a*ht)/(hx*hx);

                d[0] += (2.0*a*a*ht*lambda1*theta)/(hx) - ((2.0*a*a*ht)/(hx))*g0(tn);
            }
            else if (n == N)
            {
                w[N][N-1] = -(2.0*a*a*ht)/(hx*hx);
                w[N][N-0] = 1.0 + (2.0*a*a*ht)/(hx*hx) + lambda0*ht + (2.0*a*a*ht*lambda2)/(hx);

                d[N] += (2.0*a*a*ht*lambda2*theta)/(hx) + ((2.0*a*a*ht)/(hx))*g1(tn);
            }
            else
            {
                w[n][n-1] = -(a*a*ht)/(hx*hx);
                w[n][n+0] = 1.0 + (2.0*a*a*ht)/(hx*hx) + lambda0*ht;
                w[n][n+1] = -(a*a*ht)/(hx*hx);
            }


            for (unsigned int i=0; i<Lc; i++)
            {
                double dt = delta(sn, i);
                //if ( dt > 0.0 )
                //if (n*hx == eta[i])
                {
                    for (unsigned int j=0; j<Lo; j++)
                    {
                        d[n] += -ht*k[i][j]*z[i][j] * dt;

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
                                w[n][n1] += -ht*k[i][j] * ((2.0*hx-dh)*(hx-dh)*(hx+dh)) * h32 * dt;
                            }

                            if (hx < dh && dh <= 2.0*hx)
                            {
                                w[n][n1] += -ht*k[i][j] * ((2.0*hx-dh)*(hx-dh)*(3.0*hx-dh)) * h36 * dt;
                            }
                        }
                    }

                }
            }
        }

        DoubleVector x(N+1);
        LinearEquation::GaussianElimination(w,d,x);
        IPrinter::printVector(x);

        w.clear();
        d.clear();

        for (unsigned int n=0; n<=N; n++) u[m][n] = x[n];


        //        FILE* file = fopen("D:\\data_m_1.txt", "w");
        //        IPrinter::printMatrix(10, 4, w, M, N, NULL, file);
        //        fclose(file);

        //        FILE* file1 = fopen("D:\\data_d_1.txt", "w");
        //        IPrinter::printVector(d, NULL, N, 0,0, file1);
        //        fclose(file1);
    }
}

double IProblem2::initial(const SpaceNodePDE &sn) const
{
    return 0.0;
}

double IProblem2::boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryType boundary) const
{
    return NAN;
}

double IProblem2::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    double t = tn.t;
    double x = sn.x;

    double res = x*x - 2.0 * a*a * t + lambda0 * (U(x,t)- theta);

    double W = 0.0;
    for (unsigned int i=0; i<Lc; i++)
    {
        if ( x == eta[i] )
        {
            double vi = 0.0;
            for (unsigned int j=0; j<Lo; j++)
            {
                vi += k.at(i,j) * (U(xi[j], t) - z.at(i,j));
            }
            W += vi;
        }
    }
    res -= W;

    return res;
}

double IProblem2::g0(const TimeNodePDE &tn) const
{
    return lambda1*theta;
}

double IProblem2::g1(const TimeNodePDE &tn) const
{
    double t = tn.t;

    return 2.0*t + lambda2*(t - theta);
}

double pp = 0;
double IProblem2::delta(const SpaceNodePDE &sn, unsigned int i) const
{
    Dimension sd = spaceDimension(Dimension::DimensionX);
    double hx = sd.step();
    double sigma = 3.0*hx;
    double aa = hx*(1.0/(sqrt(2.0*M_PI)*sigma)) * exp( -((sn.x-eta[i])*(sn.x-eta[i])) / (2.0*sigma*sigma) );
    if (i==0) pp+=aa;
    if (i==0) printf("%d %d %f %f\n", sn.i, i, aa, pp);
    return aa;

    if ( sn.x == eta[i] ) return 1.0;
    return 0.0;
}

double IProblem2::U(double x, double t) const
{
    return x*x*t;
}





















